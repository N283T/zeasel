// Threading infrastructure for parallel pipelines.
//
// Provides a start barrier (master/worker synchronization), a thread pool
// (master/worker model), a work queue (producer-consumer mailbox pattern),
// and a dual work queue (reader/worker buffer-recycling pattern) for
// HMMER-style parallel search.

const std = @import("std");
const Allocator = std.mem.Allocator;

/// Master/worker start synchronization barrier, modeled after Easel's
/// ESL_THREADS start protocol (esl_threads_WaitForStart / esl_threads_Started).
///
/// Workers call `started()` to signal readiness, then block until the master
/// calls `waitForStart()`. This guarantees all workers are initialized before
/// any work begins.
pub const StartBarrier = struct {
    mutex: std.Thread.Mutex,
    cond: std.Thread.Condition,
    started_count: usize,
    n_workers: usize,
    go_signal: bool,

    pub fn init(n_workers: usize) StartBarrier {
        return .{
            .mutex = .{},
            .cond = .{},
            .started_count = 0,
            .n_workers = n_workers,
            .go_signal = false,
        };
    }

    /// Block the master thread until all workers have called `started()`,
    /// then release them all to begin work.
    pub fn waitForStart(self: *StartBarrier) void {
        self.mutex.lock();
        defer self.mutex.unlock();

        while (self.started_count < self.n_workers) {
            self.cond.wait(&self.mutex);
        }

        // All workers are ready; release them.
        self.go_signal = true;
        self.cond.broadcast();
    }

    /// Called by a worker thread to signal that it has initialized and is
    /// ready to begin work. Blocks until the master calls `waitForStart()`.
    pub fn started(self: *StartBarrier) void {
        self.mutex.lock();
        defer self.mutex.unlock();

        self.started_count += 1;
        self.cond.broadcast();

        // Wait for the master's go signal.
        while (!self.go_signal) {
            self.cond.wait(&self.mutex);
        }
    }
};

/// Per-worker argument bundle passed to the worker function.
/// Contains a unique worker index, the shared context, and a pointer
/// to the start barrier for synchronization.
pub fn WorkerInfo(comptime Context: type) type {
    return struct {
        worker_index: usize,
        context: Context,
        barrier: *StartBarrier,
    };
}

/// A thread pool that manages a fixed number of worker threads with
/// master/worker start synchronization, modeled after Easel's ESL_THREADS.
///
/// Workers receive a `WorkerInfo` containing a unique index (0..n_workers-1),
/// a shared context, and a pointer to a `StartBarrier`. Workers should call
/// `info.barrier.started()` early in their execution; the master calls
/// `pool.waitForStart()` to block until all workers are ready, then
/// releases them all simultaneously.
pub fn ThreadPool(comptime Context: type, comptime worker_fn: fn (*WorkerInfo(Context)) void) type {
    const Info = WorkerInfo(Context);

    return struct {
        const Self = @This();

        threads: []std.Thread,
        worker_infos: []Info,
        barrier: *StartBarrier,
        allocator: Allocator,

        /// Create a thread pool with `n_workers` threads.
        ///
        /// Each worker receives a `WorkerInfo` struct containing its unique index
        /// (0..n_workers-1), the shared context, and a pointer to the barrier.
        ///
        /// Workers should call `info.barrier.started()` early in their execution.
        /// The master then calls `pool.waitForStart()` to block until all workers
        /// are ready, then releases them all simultaneously.
        pub fn init(
            allocator: Allocator,
            n_workers: usize,
            context: Context,
        ) !Self {
            const threads = try allocator.alloc(std.Thread, n_workers);
            errdefer allocator.free(threads);
            const worker_infos = try allocator.alloc(Info, n_workers);
            errdefer allocator.free(worker_infos);
            const barrier = try allocator.create(StartBarrier);
            errdefer allocator.destroy(barrier);
            barrier.* = StartBarrier.init(n_workers);

            var spawned: usize = 0;
            errdefer {
                // If spawn fails partway, signal already-started workers to proceed
                // and join them so they can exit cleanly.
                barrier.mutex.lock();
                barrier.go_signal = true;
                barrier.cond.broadcast();
                barrier.mutex.unlock();
                for (0..spawned) |i| threads[i].join();
            }

            for (0..n_workers) |i| {
                worker_infos[i] = Info{
                    .worker_index = i,
                    .context = context,
                    .barrier = barrier,
                };
                threads[i] = try std.Thread.spawn(.{}, worker_fn, .{&worker_infos[i]});
                spawned += 1;
            }

            return Self{
                .threads = threads,
                .worker_infos = worker_infos,
                .barrier = barrier,
                .allocator = allocator,
            };
        }

        /// Block the master thread until all workers have called
        /// `barrier.started()`, then release them all to begin work.
        pub fn waitForStart(self: *Self) void {
            self.barrier.waitForStart();
        }

        /// Return the number of worker threads in the pool.
        pub fn workerCount(self: *const Self) usize {
            return self.threads.len;
        }

        /// Wait for all worker threads to complete and free resources.
        pub fn deinit(self: *Self) void {
            for (self.threads) |t| t.join();
            self.allocator.destroy(self.barrier);
            self.allocator.free(self.worker_infos);
            self.allocator.free(self.threads);
        }
    };
}

/// A thread-safe work queue for producer-consumer patterns.
/// Items are sent by the producer and received by consumers.
/// Uses a circular buffer for O(1) enqueue and dequeue.
pub fn WorkQueue(comptime T: type) type {
    return struct {
        const Self = @This();

        buf: []T,
        head: usize,
        count: usize,
        mutex: std.Thread.Mutex,
        not_empty: std.Thread.Condition,
        closed: bool,
        allocator: Allocator,

        pub fn init(allocator: Allocator, capacity: usize) !Self {
            const buf = try allocator.alloc(T, capacity);
            return Self{
                .buf = buf,
                .head = 0,
                .count = 0,
                .mutex = .{},
                .not_empty = .{},
                .closed = false,
                .allocator = allocator,
            };
        }

        /// Send an item to the queue. Wakes one waiting consumer.
        /// Returns error.Overflow if the queue is full.
        pub fn send(self: *Self, item: T) !void {
            self.mutex.lock();
            defer self.mutex.unlock();
            if (self.count >= self.buf.len) return error.Overflow;
            const tail = (self.head + self.count) % self.buf.len;
            self.buf[tail] = item;
            self.count += 1;
            self.not_empty.signal();
        }

        /// Receive an item from the queue. Blocks until an item is available
        /// or the queue is closed. Returns null if the queue is closed and empty.
        pub fn receive(self: *Self) ?T {
            self.mutex.lock();
            defer self.mutex.unlock();
            while (self.count == 0) {
                if (self.closed) return null;
                self.not_empty.wait(&self.mutex);
            }
            const item = self.buf[self.head];
            self.head = (self.head + 1) % self.buf.len;
            self.count -= 1;
            return item;
        }

        /// Close the queue. All waiting consumers will be woken and
        /// will receive null after draining remaining items.
        pub fn close(self: *Self) void {
            self.mutex.lock();
            defer self.mutex.unlock();
            self.closed = true;
            self.not_empty.broadcast();
        }

        /// Return the number of items currently in the queue.
        pub fn len(self: *Self) usize {
            self.mutex.lock();
            defer self.mutex.unlock();
            return self.count;
        }

        pub fn deinit(self: *Self) void {
            self.allocator.free(self.buf);
        }
    };
}

/// A dual work queue for zero-copy buffer recycling between a reader and workers.
///
/// Modeled after Easel's ESL_WORK_QUEUE. Two internal circular queues:
///   - reader queue: holds processed/empty items ready for the reader to fill
///   - worker queue: holds filled items ready for workers to process
///
/// The reader calls `readerUpdate(old, &new)` to return a processed item
/// and get an empty one. Workers call `workerUpdate(old, &new)` to return
/// a filled item and get the next one to process.
pub fn DualWorkQueue(comptime T: type) type {
    return struct {
        const Self = @This();

        reader_buf: []T,
        reader_head: usize,
        reader_count: usize,

        worker_buf: []T,
        worker_head: usize,
        worker_count: usize,

        capacity: usize,
        pending_workers: usize,

        mutex: std.Thread.Mutex,
        reader_cond: std.Thread.Condition,
        worker_cond: std.Thread.Condition,

        allocator: Allocator,

        pub fn init(allocator: Allocator, capacity: usize) !Self {
            const reader_buf = try allocator.alloc(T, capacity);
            errdefer allocator.free(reader_buf);
            const worker_buf = try allocator.alloc(T, capacity);

            return Self{
                .reader_buf = reader_buf,
                .reader_head = 0,
                .reader_count = 0,
                .worker_buf = worker_buf,
                .worker_head = 0,
                .worker_count = 0,
                .capacity = capacity,
                .pending_workers = 0,
                .mutex = .{},
                .reader_cond = .{},
                .worker_cond = .{},
                .allocator = allocator,
            };
        }

        /// Seed the reader queue with an initial item (call before starting workers).
        /// This is analogous to Easel's esl_workqueue_Init.
        pub fn addItem(self: *Self, item: T) !void {
            self.mutex.lock();
            defer self.mutex.unlock();
            if (self.reader_count >= self.capacity) return error.Overflow;
            const tail = (self.reader_head + self.reader_count) % self.capacity;
            self.reader_buf[tail] = item;
            self.reader_count += 1;
            self.reader_cond.signal();
        }

        /// Reader (producer) update: atomically enqueue a filled item onto
        /// the worker queue and dequeue the next empty item from the reader queue.
        ///
        /// Pass `null` for `in_item` on the first call to just get an empty item.
        /// Pass `null` for `out` if you don't need a return item (e.g. final flush).
        pub fn readerUpdate(self: *Self, in_item: ?T, out: ?*T) !void {
            self.mutex.lock();
            defer self.mutex.unlock();

            // Place the filled item on the worker queue.
            if (in_item) |item| {
                if (self.worker_count >= self.capacity) return error.Overflow;
                const tail = (self.worker_head + self.worker_count) % self.capacity;
                self.worker_buf[tail] = item;
                self.worker_count += 1;

                // Wake any workers waiting for work.
                if (self.pending_workers != 0) {
                    self.worker_cond.broadcast();
                }
            }

            // Wait for an empty item on the reader queue.
            if (out) |ptr| {
                while (self.reader_count == 0) {
                    self.reader_cond.wait(&self.mutex);
                }
                ptr.* = self.reader_buf[self.reader_head];
                self.reader_head = (self.reader_head + 1) % self.capacity;
                self.reader_count -= 1;
            }
        }

        /// Worker (consumer) update: atomically return a processed item to the
        /// reader queue and get the next filled item from the worker queue.
        ///
        /// Pass `null` for `in_item` on the first call to just get an item.
        pub fn workerUpdate(self: *Self, in_item: ?T, out: ?*T) void {
            self.mutex.lock();
            defer self.mutex.unlock();

            // Return the processed item to the reader queue.
            if (in_item) |item| {
                const tail = (self.reader_head + self.reader_count) % self.capacity;
                self.reader_buf[tail] = item;
                const was_empty = self.reader_count == 0;
                self.reader_count += 1;
                if (was_empty) {
                    self.reader_cond.signal();
                }
            }

            // Get the next filled item from the worker queue.
            if (out) |ptr| {
                if (self.worker_count == 0) {
                    self.pending_workers += 1;
                    while (self.worker_count == 0) {
                        self.worker_cond.wait(&self.mutex);
                    }
                    self.pending_workers -= 1;
                }

                ptr.* = self.worker_buf[self.worker_head];
                self.worker_head = (self.worker_head + 1) % self.capacity;
                self.worker_count -= 1;
            }
        }

        /// Signal all waiting workers to wake up. Call after enqueueing
        /// sentinel/termination items so workers can detect shutdown.
        pub fn complete(self: *Self) void {
            self.mutex.lock();
            defer self.mutex.unlock();
            if (self.pending_workers != 0) {
                self.worker_cond.broadcast();
            }
        }

        /// Reset the queue for reuse by moving all worker-queue items
        /// back to the reader queue.
        pub fn reset(self: *Self) void {
            self.mutex.lock();
            defer self.mutex.unlock();

            while (self.worker_count > 0) {
                const tail = (self.reader_head + self.reader_count) % self.capacity;
                self.reader_buf[tail] = self.worker_buf[self.worker_head];
                self.reader_count += 1;
                self.worker_head = (self.worker_head + 1) % self.capacity;
                self.worker_count -= 1;
            }

            self.pending_workers = 0;
        }

        pub fn deinit(self: *Self) void {
            self.allocator.free(self.reader_buf);
            self.allocator.free(self.worker_buf);
        }
    };
}

/// Return the number of available hardware threads.
pub fn cpuCount() usize {
    return std.Thread.getCpuCount() catch 1;
}

// --- Tests ---

test "WorkQueue: send and receive" {
    var q = try WorkQueue(u32).init(std.testing.allocator, 16);
    defer q.deinit();

    try q.send(42);
    try q.send(99);

    try std.testing.expectEqual(@as(?u32, 42), q.receive());
    try std.testing.expectEqual(@as(?u32, 99), q.receive());
}

test "WorkQueue: close returns null" {
    var q = try WorkQueue(u32).init(std.testing.allocator, 16);
    defer q.deinit();

    try q.send(1);
    q.close();

    try std.testing.expectEqual(@as(?u32, 1), q.receive());
    try std.testing.expectEqual(@as(?u32, null), q.receive());
}

test "WorkQueue: overflow returns error" {
    var q = try WorkQueue(u32).init(std.testing.allocator, 2);
    defer q.deinit();

    try q.send(1);
    try q.send(2);
    try std.testing.expectError(error.Overflow, q.send(3));
}

test "WorkQueue: threaded producer-consumer" {
    var q = try WorkQueue(u32).init(std.testing.allocator, 16);
    defer q.deinit();

    const producer = try std.Thread.spawn(.{}, struct {
        fn run(queue: *WorkQueue(u32)) void {
            for (0..10) |i| {
                queue.send(@intCast(i)) catch {};
            }
            queue.close();
        }
    }.run, .{&q});

    var sum: u32 = 0;
    while (q.receive()) |v| {
        sum += v;
    }
    producer.join();

    try std.testing.expectEqual(@as(u32, 45), sum); // 0+1+...+9
}

test "WorkQueue: FIFO ordering with wraparound" {
    var q = try WorkQueue(u32).init(std.testing.allocator, 4);
    defer q.deinit();

    // Fill and drain to advance head past start
    try q.send(10);
    try q.send(20);
    try std.testing.expectEqual(@as(?u32, 10), q.receive());
    try std.testing.expectEqual(@as(?u32, 20), q.receive());

    // Now head is at index 2; fill all 4 slots to force wraparound
    try q.send(1);
    try q.send(2);
    try q.send(3);
    try q.send(4);

    try std.testing.expectEqual(@as(?u32, 1), q.receive());
    try std.testing.expectEqual(@as(?u32, 2), q.receive());
    try std.testing.expectEqual(@as(?u32, 3), q.receive());
    try std.testing.expectEqual(@as(?u32, 4), q.receive());
}

test "DualWorkQueue: single-threaded round-trip" {
    var q = try DualWorkQueue(u32).init(std.testing.allocator, 4);
    defer q.deinit();

    // Seed two empty buffers (values don't matter yet)
    try q.addItem(0);
    try q.addItem(0);

    // Reader gets an empty buffer, fills it, sends it
    var buf: u32 = undefined;
    try q.readerUpdate(null, &buf); // get empty buffer
    buf = 42;
    try q.readerUpdate(buf, &buf); // send filled, get next empty
    buf = 99;
    try q.readerUpdate(buf, null); // send filled, don't wait for empty

    // Worker gets filled items
    var work: u32 = undefined;
    q.workerUpdate(null, &work); // get first filled item
    try std.testing.expectEqual(@as(u32, 42), work);
    q.workerUpdate(work, &work); // return processed, get next
    try std.testing.expectEqual(@as(u32, 99), work);
    q.workerUpdate(work, null); // return last, don't wait
}

test "DualWorkQueue: threaded reader-worker pipeline" {
    const ncpu = 2;
    const n_items = 20;

    var q = try DualWorkQueue(u32).init(std.testing.allocator, ncpu * 2);
    defer q.deinit();

    // Seed empty buffers
    for (0..ncpu * 2) |_| {
        try q.addItem(0);
    }

    var result_sum = std.atomic.Value(u32).init(0);

    // Spawn worker threads
    var workers: [ncpu]std.Thread = undefined;
    for (0..ncpu) |i| {
        workers[i] = try std.Thread.spawn(.{}, struct {
            fn run(queue: *DualWorkQueue(u32), sum: *std.atomic.Value(u32)) void {
                var item: u32 = undefined;
                queue.workerUpdate(null, &item);
                while (item > 0) {
                    _ = sum.fetchAdd(item, .monotonic);
                    queue.workerUpdate(item, &item);
                }
                // Return the sentinel
                queue.workerUpdate(item, null);
            }
        }.run, .{ &q, &result_sum });
    }

    // Reader: fill items 1..n_items, then send ncpu sentinels (0)
    var buf: u32 = undefined;
    try q.readerUpdate(null, &buf);
    for (1..n_items + 1) |i| {
        buf = @intCast(i);
        try q.readerUpdate(buf, &buf);
    }
    // Send sentinel for each worker
    for (0..ncpu) |_| {
        buf = 0;
        try q.readerUpdate(buf, &buf);
    }

    for (0..ncpu) |i| workers[i].join();

    // Sum of 1..20 = 210
    try std.testing.expectEqual(@as(u32, 210), result_sum.load(.monotonic));
}

test "DualWorkQueue: reset moves items back to reader queue" {
    var q = try DualWorkQueue(u32).init(std.testing.allocator, 4);
    defer q.deinit();

    try q.addItem(1);
    try q.addItem(2);

    // Reader sends both to worker queue
    var buf: u32 = undefined;
    try q.readerUpdate(null, &buf);
    buf = 10;
    try q.readerUpdate(buf, &buf);
    buf = 20;
    try q.readerUpdate(buf, null);

    // Reset moves worker items back to reader queue
    q.reset();

    // Reader should be able to get items back
    try q.readerUpdate(null, &buf);
    try std.testing.expect(buf == 10 or buf == 20);
}

test "ThreadPool: basic spawn and join lifecycle" {
    const Pool = ThreadPool(void, struct {
        fn worker(info: *WorkerInfo(void)) void {
            info.barrier.started();
        }
    }.worker);

    var pool = try Pool.init(std.testing.allocator, 4, {});
    pool.waitForStart();
    pool.deinit();
}

test "ThreadPool: all workers execute with unique indices" {
    const n_workers = 8;
    var executed = [_]std.atomic.Value(u32){std.atomic.Value(u32).init(0)} ** n_workers;

    const Ctx = struct {
        flags: *[n_workers]std.atomic.Value(u32),
    };
    const Pool = ThreadPool(Ctx, struct {
        fn worker(info: *WorkerInfo(Ctx)) void {
            info.barrier.started();
            info.context.flags[info.worker_index].store(1, .release);
        }
    }.worker);

    var pool = try Pool.init(std.testing.allocator, n_workers, Ctx{ .flags = &executed });
    pool.waitForStart();
    pool.deinit();

    // Verify all workers executed exactly once with unique indices.
    for (0..n_workers) |i| {
        try std.testing.expectEqual(@as(u32, 1), executed[i].load(.acquire));
    }
}

test "ThreadPool: workers receive correct context and indices" {
    const n_workers = 4;
    var index_sum = std.atomic.Value(u64).init(0);
    var context_sum = std.atomic.Value(u64).init(0);

    const Ctx = struct {
        index_sum: *std.atomic.Value(u64),
        context_sum: *std.atomic.Value(u64),
        magic: u64,
    };
    const Pool = ThreadPool(Ctx, struct {
        fn worker(info: *WorkerInfo(Ctx)) void {
            info.barrier.started();
            _ = info.context.index_sum.fetchAdd(@intCast(info.worker_index), .monotonic);
            _ = info.context.context_sum.fetchAdd(info.context.magic, .monotonic);
        }
    }.worker);

    var pool = try Pool.init(std.testing.allocator, n_workers, Ctx{
        .index_sum = &index_sum,
        .context_sum = &context_sum,
        .magic = 42,
    });
    pool.waitForStart();
    pool.deinit();

    // Indices 0+1+2+3 = 6
    try std.testing.expectEqual(@as(u64, 6), index_sum.load(.monotonic));
    // Each of 4 workers adds 42 => 168
    try std.testing.expectEqual(@as(u64, 168), context_sum.load(.monotonic));
}

test "ThreadPool: work distribution via WorkQueue" {
    const n_workers = 4;
    const n_items = 100;

    var queue = try WorkQueue(u32).init(std.testing.allocator, n_items + n_workers);
    defer queue.deinit();

    var result_sum = std.atomic.Value(u64).init(0);

    const Ctx = struct {
        queue: *WorkQueue(u32),
        result_sum: *std.atomic.Value(u64),
    };
    const Pool = ThreadPool(Ctx, struct {
        fn worker(info: *WorkerInfo(Ctx)) void {
            info.barrier.started();
            while (info.context.queue.receive()) |item| {
                _ = info.context.result_sum.fetchAdd(item, .monotonic);
            }
        }
    }.worker);

    var pool = try Pool.init(std.testing.allocator, n_workers, Ctx{
        .queue = &queue,
        .result_sum = &result_sum,
    });
    pool.waitForStart();

    // Produce work items after all workers are ready.
    for (1..n_items + 1) |i| {
        queue.send(@intCast(i)) catch unreachable;
    }
    queue.close();

    pool.deinit();

    // Sum of 1..100 = 5050
    try std.testing.expectEqual(@as(u64, 5050), result_sum.load(.monotonic));
}

test "ThreadPool: start barrier ensures workers wait for go signal" {
    const n_workers = 4;
    var pre_barrier = std.atomic.Value(u32).init(0);
    var post_barrier = std.atomic.Value(u32).init(0);

    const Ctx = struct {
        pre_barrier: *std.atomic.Value(u32),
        post_barrier: *std.atomic.Value(u32),
    };
    const Pool = ThreadPool(Ctx, struct {
        fn worker(info: *WorkerInfo(Ctx)) void {
            _ = info.context.pre_barrier.fetchAdd(1, .monotonic);
            info.barrier.started();
            _ = info.context.post_barrier.fetchAdd(1, .monotonic);
        }
    }.worker);

    var pool = try Pool.init(std.testing.allocator, n_workers, Ctx{
        .pre_barrier = &pre_barrier,
        .post_barrier = &post_barrier,
    });

    // Wait for all workers to reach the barrier.
    pool.waitForStart();

    // At this point, all workers have incremented pre_barrier.
    try std.testing.expectEqual(@as(u32, n_workers), pre_barrier.load(.monotonic));

    // Workers are now released; wait for them to finish.
    pool.deinit();

    // All workers should have incremented post_barrier.
    try std.testing.expectEqual(@as(u32, n_workers), post_barrier.load(.monotonic));
}

test "cpuCount: at least 1" {
    try std.testing.expect(cpuCount() >= 1);
}
