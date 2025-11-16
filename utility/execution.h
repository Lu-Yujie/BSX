#include <iostream>
#include <future>
#include <chrono>
#include <thread>
#include <stdexcept>

// execution with time limited by detaching it, but it will not stop the thread
// for execution a function safely, write a task status detection in func is suggestioned
// stopping the func thread forcely may cause data safety and other problems
template<typename ReturnType, typename Function, typename... Args>
ReturnType execute_with_timeout(std::chrono::milliseconds timeout_duration,
                                ReturnType default_value,
                                Function&& func, Args&&... args) {
    // use std::promise to receive execution results
    std::promise<ReturnType> resultPromise;

    // start a worker
    std::thread worker_thread([&]() {
        try {
            ReturnType result = func(std::forward<Args>(args)...);
            resultPromise.set_value(result); // write result promise
        } catch (...) {
            resultPromise.set_exception(std::current_exception()); // set exception promise
        }
    });

    // wait for limited_time and detect status
    auto future = resultPromise.get_future();
    if (future.wait_for(timeout_duration) == std::future_status::timeout) {
        // if overtime, cancel the task
        worker_thread.detach();  // just detach the worker_thread
        throw std::runtime_error("Task timed out and was forcefully terminated.");
    }

    ReturnType result = default_value;
    try {
        result = future.get();  // get result or throw exception
    } catch (const std::exception& e) {
        throw;
    }

    worker_thread.join(); // wait for worker_thread ends

    return result;
}

int sleep10(int a, int b) {
    std::this_thread::sleep_for(std::chrono::seconds(10));
    return a + b;
}

int exe_test() {
    std::chrono::milliseconds timeout_duration(1000);
    try {
        int result = execute_with_timeout(timeout_duration, -1, sleep10, 10, 20);
        std::cout << "Result: " << result << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
