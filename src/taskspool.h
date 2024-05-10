#pragma once

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>

class TasksPool {
public:
    void init(size_t nt);
    template<class Func, class... Args>
    auto add(Func&& f, Args&&... args) -> std::future<typename std::result_of<Func(Args...)>::type>;
    void exit();
    ~TasksPool();
private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool done = false;
};
 
inline void TasksPool::init(size_t nt) {
    for(size_t i=0; i<nt; i++) {
        workers.emplace_back([this]{
            while(true) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(this->queue_mutex);
                    this->condition.wait(lock, [this]{ return this->done || !this->tasks.empty(); });
                    if(this->done && this->tasks.empty()) return;
                    task = std::move(this->tasks.front());
                    this->tasks.pop();
                }
                task();
            }
        });
    }
}

template<class Func, class... Args>
auto TasksPool::add(Func&& f, Args&&... args) -> std::future<typename std::result_of<Func(Args...)>::type> {
    using return_type = typename std::result_of<Func(Args...)>::type;
    auto task = std::make_shared<std::packaged_task<return_type()>>(std::bind(std::forward<Func>(f), std::forward<Args>(args)...));
    std::future<return_type> res = task->get_future();
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        if(done) throw std::runtime_error("tasks exit, don't not add any more.");
        tasks.emplace([task](){ (*task)(); });
    }
    condition.notify_one();
    return res;
}

inline void TasksPool::exit() {
    if(done) return;
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        done = true;
    }
    condition.notify_all();
    for(auto &worker : workers) worker.join();
}

inline TasksPool::~TasksPool() {
    exit();
}
