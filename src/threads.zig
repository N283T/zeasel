// Threading infrastructure for parallel pipelines.
//
// Provides a thread pool (master/worker model) and a work queue
// (producer-consumer mailbox pattern) for HMMER-style parallel search.

const std = @import("std");
const Allocator = std.mem.Allocator;

/// A simple thread pool that manages a fixed number of worker threads.
pub fn ThreadPool(comptime Context: type) type {
    return struct {
        const Self = @This();

        threads: []std.Thread,
        allocator: Allocator,

        /// Create a thread pool with `n_workers` threads, each running `worker_fn`.
        /// The context is passed to each worker. Workers run until they return.
        pub fn init(
            allocator: Allocator,
            n_workers: usize,
            worker_fn: *const fn (Context) void,
            context: Context,
        ) !Self {
            const threads = try allocator.alloc(std.Thread, n_workers);
            var spawned: usize = 0;
            errdefer {
                for (0..spawned) |i| threads[i].join();
                allocator.free(threads);
            }

            for (0..n_workers) |i| {
                threads[i] = try std.Thread.spawn(.{}, worker_fn, .{context});
                spawned += 1;
            }

            return Self{
                .threads = threads,
                .allocator = allocator,
            };
        }

        /// Wait for all worker threads to complete and free resources.
        pub fn deinit(self: *Self) void {
            for (self.threads) |t| t.join();
            self.allocator.free(self.threads);
        }
    };
}

/// A thread-safe work queue for producer-consumer patterns.
/// Items are sent by the producer and received by consumers.
pub fn WorkQueue(comptime T: type) type {
    return struct {
        const Self = @This();

        items: std.ArrayList(T),
        mutex: std.Thread.Mutex,
        not_empty: std.Thread.Condition,
        closed: bool,
        allocator: Allocator,

        pub fn init(allocator: Allocator) Self {
            return Self{
                .items = .empty,
                .mutex = .{},
                .not_empty = .{},
                .closed = false,
                .allocator = allocator,
            };
        }

        /// Send an item to the queue. Wakes one waiting consumer.
        pub fn send(self: *Self, item: T) !void {
            self.mutex.lock();
            defer self.mutex.unlock();
            try self.items.append(self.allocator, item);
            self.not_empty.signal();
        }

        /// Receive an item from the queue. Blocks until an item is available
        /// or the queue is closed. Returns null if the queue is closed and empty.
        pub fn receive(self: *Self) ?T {
            self.mutex.lock();
            defer self.mutex.unlock();
            while (self.items.items.len == 0) {
                if (self.closed) return null;
                self.not_empty.wait(&self.mutex);
            }
            // Pop from front (FIFO)
            const item = self.items.orderedRemove(0);
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
            return self.items.items.len;
        }

        pub fn deinit(self: *Self) void {
            self.items.deinit(self.allocator);
        }
    };
}

/// Return the number of available hardware threads.
pub fn cpuCount() usize {
    return std.Thread.getCpuCount() catch 1;
}

// --- Tests ---

test "WorkQueue: send and receive" {
    var q = WorkQueue(u32).init(std.testing.allocator);
    defer q.deinit();

    try q.send(42);
    try q.send(99);

    try std.testing.expectEqual(@as(?u32, 42), q.receive());
    try std.testing.expectEqual(@as(?u32, 99), q.receive());
}

test "WorkQueue: close returns null" {
    var q = WorkQueue(u32).init(std.testing.allocator);
    defer q.deinit();

    try q.send(1);
    q.close();

    try std.testing.expectEqual(@as(?u32, 1), q.receive());
    try std.testing.expectEqual(@as(?u32, null), q.receive());
}

test "WorkQueue: threaded producer-consumer" {
    var q = WorkQueue(u32).init(std.testing.allocator);
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

test "cpuCount: at least 1" {
    try std.testing.expect(cpuCount() >= 1);
}
