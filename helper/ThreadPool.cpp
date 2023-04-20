#include "ThreadPool.hpp"

// Constructor
ThreadPool::ThreadPool(size_t num_threads) : stop_flag(false), unfinished_tasks(0)
{
    // Create worker threads
    for (size_t i = 0; i < num_threads; ++i)
    {
        threads.emplace_back([this]
                             {
            // Keep processing tasks until the stop flag is set and the task queue is empty
            while(true) {
                std::unique_lock<std::mutex> lock(mutex);

                // Wait for a task to be added to the queue or the stop flag to be set
                condition.wait(lock, [this] {return !task_queue.empty() || stop_flag;});
                
                // If the stop flag is set and the task queue is empty, break the loop and exit the thread
                if (stop_flag && task_queue.empty()) { return;}

                // Retrieve the next task from the queue
                auto task = std::move(task_queue.front());
                task_queue.pop();

                // Unlock the mutex and run the task
                lock.unlock();
                task();
            } });
    }
}

// Destructor
ThreadPool::~ThreadPool()
{
    {
        std::unique_lock<std::mutex> lock(mutex);

        // Set the stop flag to true to signal to the worker threads to exit 
        stop_flag = true;
    }

    // Notify all worker threads that the stop flag has been set
    condition.notify_all();

    // Join all worker threads to wait for them to finish executing
    for (auto &thread : threads)
    {
        thread.join();
    }
}


// Enqueue function template
template <typename F, typename... Args>
void ThreadPool::enqueue(F &&f, Args &&...args)
{
    {
        std::unique_lock<std::mutex> lock(mutex);

        // Add a new task to the task queue with a lambda function that runs the task
        // and decrements the unfinished_tasks counter when the task is completed
        task_queue.emplace([this, task = std::bind(std::forward<F>(f), std::forward<Args>(args)...)]()
                           {
            task();
            std::unique_lock<std::mutex> lock(mutex);
            --unfinished_tasks;
            finished_condition.notify_one(); });
        ++unfinished_tasks;
    }

    // Notify one worker thread that a new task has been added to the queue
    condition.notify_one();
}

void ThreadPool::wait()
{
    std::unique_lock<std::mutex> lock(mutex);
    finished_condition.wait(lock, [this]
                            { return task_queue.empty() && unfinished_tasks == 0; });
}

// Explicitly instantiate enqueue for the supported types
template void ThreadPool::enqueue<void (&)(int, int, std::__1::atomic<int> &), unsigned int &, unsigned int &, std::__1::reference_wrapper<std::__1::atomic<int>>>(void (&)(int, int, std::__1::atomic<int> &), unsigned int &, unsigned int &, std::__1::reference_wrapper<std::__1::atomic<int>> &&);


void print_progress(int finished, int total, std::mutex &mtx) {
    const int bar_width = 50;
    float progress = static_cast<float>(finished) / total;
    int pos = static_cast<int>(bar_width * progress);

    std::unique_lock<std::mutex> lock(mtx);
    std::cout << "Rendering Progress: [";

    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }

    std::cout << "] " << finished << "/" << total << " (" << static_cast<int>(progress * 100.0) << "%)\r";
    
}