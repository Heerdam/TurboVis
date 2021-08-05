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

using Func =  std::function<void(size_t /*_index*/, size_t /*_x*/, size_t /*_y*/)>;

class AsyncRenderer {

    Func cb = 0;

    size_t index = 0;
    size_t width, height;

    bool b_isShutdown = true;
    std::vector<bool> brestart;

    std::vector<std::thread> worker;

    std::mutex mutex;
    std::condition_variable cv;

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
    void work(size_t /*_index*/) noexcept;

};

// ----------------- PUBLIC -----------------

inline double AsyncRenderer::getProgress() noexcept {
    std::lock_guard<std::mutex> lock(mutex);
    return double(index) / (double)(width * height);
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

    for(size_t i = 0; i < _numThreads; ++i){
        brestart.push_back(false);
    }

    for(size_t i = 0; i < _numThreads; ++i){
        worker.emplace_back(&AsyncRenderer::work, this, i);
    }

};

inline void AsyncRenderer::stop() noexcept {
    {
        std::lock_guard<std::mutex> lock(mutex);
        b_isShutdown = true;
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
inline void AsyncRenderer::work(size_t _index) noexcept {

    std::unique_lock<std::mutex> lock (mutex, std::defer_lock);
    size_t curIdx = 0;
    bool endOfWork = false;

    if(_index == 0)
        std::cout << "started" << std::endl;;

    while(!b_isShutdown){

        while(true){

            //aquire new index
            lock.lock();
            curIdx = index;
            index += 1;
            endOfWork = curIdx >= width * height;
            lock.unlock();
            
 
            

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
            cb(_index, x, y);

            //break if we have a restart
            lock.lock();
            const bool br = brestart[_index];
            lock.unlock();
            if(br) break;

        }

        //wait for restart or termination
        if(!brestart[_index]){  
            if(_index == 0)
                std::cout << "waiting for restart [" << _index << "]" << std::endl;   
            lock.lock();
            cv.wait(lock, [&](){ return brestart[_index] || b_isShutdown; });       
            brestart[_index] = false;
            lock.unlock();
        }else if(brestart[_index])
            brestart[_index] = false;

        //if(_index == 0)
           //std::cout << std::endl;

    }

};

#endif /* ASYNCRENDERER_HPP */