#include "ThreadPool.hpp"

// Constructor
// The stop_flag is initialized with false and unfinished_tasks is initialized with 0
ThreadPool::ThreadPool(size_t num_threads) : stop_flag(false), unfinished_tasks(0)
{
    // Create worker threads
    for (size_t i = 0; i < num_threads; ++i)
    {
        // Each worker thread runs a lambda function that captures "this" pointer, allowing the lambda
        // to access the ThreadPool's members
        threads.emplace_back([this]
            {
                // Keep processing tasks until the stop flag is set and the task queue is empty
                while(true) {

                    // acquire lock on mutex, ensure thread safety when accessing task_queue and stop_flag
                    std::unique_lock<std::mutex> lock(mutex);

                    // Wait for a task to be added to the queue or the stop flag to be set
                    // the worker thread will be blocked until either of these conditions is true
                    condition.wait(lock, [this] {return !task_queue.empty() || stop_flag;});
                    
                    // after the condition variable is signaled and the thread wakes up, the thread checks if the stop_flag is set and task_queue is empty
                    // If so, break the loop and exit the thread
                    if (stop_flag && task_queue.empty()) {return;}

                    // If the worker thread continues, it retrieves the next task from the front of the task_queue and removes it from the queue
                    auto task = std::move(task_queue.front());
                    task_queue.pop();

                    // Unlock the mutex and run the task
                    lock.unlock();
                    task();
                } 
            });
    }
}

// Destructor
// Destructor is responsible for properly shutting down the ThreadPool when it is being destroyed
// It ensures that all worker threads are signaled to stop and that the destructor waits for them to finish 
// before the ThreadPool Object is destroyed 
ThreadPool::~ThreadPool()
{
    {
        // acquire a lock on mtex to ensure thread safety when modifying stop_flag
        std::unique_lock<std::mutex> lock(mutex);

        // Set the stop flag to true to signal to the worker threads to exit
        stop_flag = true;
    }
    // the lock on mutex is released

    // Wake up all worker threads.  So that the worker threads can check the updated stop_flag and exit if necessary
    condition.notify_all();

    // Join all worker threads to wait for them to finish executing
    for (auto &thread : threads)
    {
        thread.join();
    }
}

// Enqueue function template
// The function is a template functions that takes a callable F and a variadic template paramter pack Args
// This allows any function or callable object to be enqueued with its respective arguments
template <typename F, typename... Args>
void ThreadPool::enqueue(F &&f, Args &&...args)
{
    {
        // acquire a lock on mutex to ensure thread safety when modifying the task queue and unifnished_tasks
        std::unique_lock<std::mutex> lock(mutex);

        // A lambda function is created which wraps the provided callable f and its arguments in std::bind expression
        // This lambda function is when added to the task_queue
        // The lambda function executes the bound callable task(), acquires a lock on the mutex to ensure thread safety
        // Decrements the unifnished_tasks and notifies the finished_condition that a task has been completed
        task_queue.emplace([this, task = std::bind(std::forward<F>(f), std::forward<Args>(args)...)]()
            {
                task();
                std::unique_lock<std::mutex> lock(mutex);
                --unfinished_tasks;
                finished_condition.notify_one(); 
            });

        // The unifhsed_tasks is incremented since a new task has been added to the queue
        ++unfinished_tasks;
    }
    // The lock on mutex is released

    // Notify one worker thread that a new task has been added to the queue, which will cause one of the worker threads to wake up and process the new task
    condition.notify_one();
}

void ThreadPool::wait()
{
    // This code waits until all tasks in the task queue has been completed.
    // Unique lock is used to lock mutex that is associated with the task queue
    // The wait() on finished_condition takes 2 arguments, an unique lock and a predicate function that returns a boolean value
    // The predicate function is called repeatedly while the condition variable is waiting.
    // The conditional variable is notified when the predicate returns true.  
    // The wait() function then returns, and the unique_lock object is automatically unlocked
    // THe predicate function captures the "this" pointer to access the task_queue and unfinished_tasks member variables of
    // the 'Threadpool' class
    std::unique_lock<std::mutex> lock(mutex);
    finished_condition.wait(lock, [this] { return task_queue.empty() && unfinished_tasks == 0; });
}

// Explicitly instantiate enqueue for the supported types
template void ThreadPool::enqueue<void (&)(int, int, std::__1::atomic<int> &), unsigned int &, unsigned int &, std::__1::reference_wrapper<std::__1::atomic<int>>>(void (&)(int, int, std::__1::atomic<int> &), unsigned int &, unsigned int &, std::__1::reference_wrapper<std::__1::atomic<int>> &&);

void print_progress(int finished, int total, std::mutex &mtx)
{
    const int bar_width = 50;
    float progress = static_cast<float>(finished) / total;
    int pos = static_cast<int>(bar_width * progress);

    // acquire the lock on mtx
    std::unique_lock<std::mutex> lock(mtx);
    std::cout << "Rendering Progress: [";

    // print the progress bar based on the  finished / toal task ratios
    for (int i = 0; i < bar_width; ++i)
    {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }

    std::cout << "] " << finished << "/" << total << " (" << static_cast<int>(progress * 100.0) << "%)\r";
}