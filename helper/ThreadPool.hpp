#ifndef THREAD_POOL_HPP
#define THREAD_POOL_HPP

#include <condition_variable>
#include <functional>
#include <iostream>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

class ThreadPool {
public:
    ThreadPool(size_t num_threads);
    ~ThreadPool();

    template<typename F, typename... Args>
    void enqueue(F&& f, Args&&... args);

    void wait();

private:
    std:: vector<std::thread> threads;
    std::queue<std::function<void()>> task_queue;
    std::mutex mutex;
    std::condition_variable condition;
    std::condition_variable finished_condition;
    bool stop_flag;
    size_t unfinished_tasks = 0;

};

#endif 