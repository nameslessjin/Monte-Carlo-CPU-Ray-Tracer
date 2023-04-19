#include "ThreadPool.hpp"

ThreadPool::ThreadPool(size_t num_threads) : stop_flag(false), unfinished_tasks(0)
{
    for (size_t i = 0; i < num_threads; ++i)
    {
        threads.emplace_back([this]
                                {
            while(true) {
                std::unique_lock<std::mutex> lock(mutex);
                condition.wait(lock, [this] {return !task_queue.empty() || stop_flag;});
                
                if (stop_flag && task_queue.empty()) { return;}

                auto task = std::move(task_queue.front());
                task_queue.pop();
                lock.unlock();
                task();
            } });
    }
}

ThreadPool::~ThreadPool() {
    {
        std::unique_lock<std::mutex> lock(mutex);
        stop_flag = true;
    }
    condition.notify_all();

    for (auto &thread: threads) {
        thread.join();
    }
}

template <typename F, typename... Args>
void ThreadPool::enqueue(F&& f, Args&&... args) {
    {
        std::unique_lock<std::mutex> lock(mutex);
        task_queue.emplace([this, task = std::bind(std::forward<F>(f), std::forward<Args>(args)...)]() {
            task();
            std::unique_lock<std::mutex> lock(mutex);
            --unfinished_tasks;
            finished_condition.notify_one();
        });
        ++unfinished_tasks;
    }
    condition.notify_one();
}

void ThreadPool::wait() {
    std::unique_lock<std::mutex> lock(mutex);
    finished_condition.wait(lock, [this] {return task_queue.empty() && unfinished_tasks == 0;});
}

template void ThreadPool::enqueue<void (&)(int, int, std::__1::atomic<int>&), unsigned int&, unsigned int&, std::__1::reference_wrapper<std::__1::atomic<int>>>(void (&)(int, int, std::__1::atomic<int>&), unsigned int&, unsigned int&, std::__1::reference_wrapper<std::__1::atomic<int>>&&);