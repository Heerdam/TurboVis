#ifndef ASYNCRENDERER_HPP
#define ASYNCRENDERER_HPP

#include <vector>
#include <thread>
#include <condition_variable>
#include <atomic>
#include <string>
#include <queue>
#include <functional>
#include <assert.h>
#include <exception>
#include <chrono>

using Func =  std::function<void(size_t /*_index*/, size_t /*_x*/, size_t /*_y*/)>;
using Now = std::chrono::high_resolution_clock::time_point;

class AsyncRenderer {

    Func cb = 0;

    size_t index = 0;
    size_t width, height;
    size_t threads;

    bool b_isShutdown = true;
    std::vector<bool> brestart;
    std::vector<int64_t> done;

    Now started;
    Now finished;

    std::atomic<int64_t> finishedCounter;

    std::vector<std::thread> worker;

    std::mutex mutex;
    std::condition_variable cv;

    [[nodiscard]] std::vector<size_t> sections(size_t /*_threads*/, size_t /*_width*/, size_t /*_height*/) const noexcept;

public:
    AsyncRenderer() noexcept;
    void start(size_t /*_numThreads*/, size_t /*_width*/, size_t /*_height*/, Func /*_cb*/) noexcept;
    void stop() noexcept;
    void restart() noexcept;
    [[nodiscard]] bool isShutdown() noexcept;
    [[nodiscard]] bool isRestart(size_t /*index*/) noexcept;
    [[nodiscard]] std::mutex& getMutex() noexcept;
    [[nodiscard]] double getProgress() noexcept;

private:
    void work(size_t /*_index*/, size_t /*_min*/, size_t /*_max*/) noexcept;

};

// ----------------- PUBLIC -----------------

inline double AsyncRenderer::getProgress() noexcept {
    const int64_t sum = std::accumulate(done.begin(), done.end(), 0ll);
    return double(sum) / (double)(width * height);
}

inline AsyncRenderer::AsyncRenderer() noexcept {}

inline std::mutex& AsyncRenderer::getMutex() noexcept {
    return mutex;
}

inline void AsyncRenderer::start(size_t _numThreads, size_t _width, size_t _height, Func _cb) noexcept { 

    stop();

    cb = _cb;
    width = _width;
    height = _height;
    index = 0;
    b_isShutdown = false;
    threads = _numThreads;
    finishedCounter = _numThreads;

    done.resize(_numThreads);

    const auto s = sections(_numThreads, _width, _height);

    for(size_t i = 0; i < _numThreads; ++i){
        brestart.push_back(false);
        done[i] = 0;
    }

    for(size_t i = 0; i < _numThreads; ++i){
        worker.emplace_back(&AsyncRenderer::work, this, i, s[i], s[i+1]);
    }

    started = std::chrono::high_resolution_clock::now();

    std::cout << "started: " << std::endl;

};

inline void AsyncRenderer::stop() noexcept {
    {
        std::lock_guard<std::mutex> lock(mutex);
        b_isShutdown = true;
        std::cout << "shutdown recieved" << std::endl;
    }

    cv.notify_all();
    for(auto& t : worker)
        t.join();

    worker.clear();
    brestart.clear();
};

inline void AsyncRenderer::restart() noexcept {
    index = 0;
    for(size_t i = 0; i < brestart.size(); ++i)
        brestart[i] = true;
    cv.notify_all();
};

inline bool AsyncRenderer::isShutdown() noexcept {
    std::lock_guard<std::mutex> lock(mutex);
    return b_isShutdown;
};

inline bool AsyncRenderer::isRestart(size_t _index) noexcept {
    std::lock_guard<std::mutex> lock(mutex);
    assert(_index < brestart.size());
    return brestart[_index];
};

// ----------------- PRIVATE -----------------
inline std::vector<size_t> AsyncRenderer::sections(size_t _threads, size_t _width, size_t _height) const noexcept {

    const size_t s = (_width * _height) / _threads;

    // [start, end)
    std::vector<size_t> out;
    out.push_back(0);

    for(size_t i = 1; i < _threads; ++i)
        out.push_back(i*s);
    
    out.push_back(_width * _height);

    return out;
}

inline void AsyncRenderer::work(size_t _index, size_t _min, size_t _max) noexcept {

    std::unique_lock<std::mutex> lock (mutex, std::defer_lock);
    bool endOfWork = false;  

    while(!b_isShutdown){

        //while(true){
        for(size_t ix = _min; ix < _max; ++ix){

            

            //aquire new index
            size_t curIdx = ix;

            // break if we have reached the end
            if(endOfWork) break;

            //do the work
            const size_t x = curIdx % width;
            const size_t y = curIdx / width;  
            
            //if(_index == 0 && curIdx % 1000 == 0){

                //std::cout << curIdx << " " << x << " " << y << std::endl;
                //const float frac = (float(curIdx) / float(width * height)) * 100.f;
             //   std::cout << curIdx << "|" << width * height << std::endl;
            //}
            Now st = std::chrono::high_resolution_clock::now();
            cb(_index, x, y);
            const double ray = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - st).count();
            //if(_index == 0)
                //std::cout << ray << std::endl;

            //break if we have a restart
            //lock.lock();
            const bool br = brestart[_index];
            //lock.unlock();
            if(br) break;

            done[_index] += 1;

        }

        //wait for restart or termination
        if(!brestart[_index]){ 
            finishedCounter -= 1; 
            if(finishedCounter == 0){
                finished = std::chrono::high_resolution_clock::now();
                std::cout << "finished after " << std::chrono::duration<double>(finished - started).count() << "s" << std::endl;
            }
            //if(_index == 0)
                //std::cout << "waiting for restart [" << _index << "]" << std::endl;   
            lock.lock();
            cv.wait(lock, [&](){ return brestart[_index] || b_isShutdown; });       
            brestart[_index] = false;
            lock.unlock();
            if(_index == 0)
                finishedCounter = threads;   
        }else if(brestart[_index])
            brestart[_index] = false;

        //if(_index == 0)
           //std::cout << std::endl;

    }

};

#endif /* ASYNCRENDERER_HPP */