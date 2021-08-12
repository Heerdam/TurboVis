#ifndef ASYNCRENDERER_LOCKFREE_HPP
#define ASYNCRENDERER_LOCKFREE_HPP

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

class AsyncRendererLF {

    std::vector<std::thread> worker;

    [[nodiscard]] std::vector<size_t> sections(size_t /*_threads*/, size_t /*_width*/, size_t /*_height*/) const noexcept;
    void work(size_t /*_index*/) noexcept;

public:

    void start(size_t /*_numThreads*/, size_t /*_width*/, size_t /*_height*/, Func /*_cb*/) noexcept;
    void stop() noexcept;
    void restart() noexcept;

};

inline std::vector<size_t> AsyncRendererLF::sections(size_t _threads, size_t _width, size_t _height) const noexcept {

    const size_t s = (_width * _height) / _threads;

    // [start, end)
    std::vector<size_t> out ( _threads + 1 );
    out.push_back(0);

    for(size_t i = 1; i < _threads; ++i)
        out.push_back(i*s);
    
    out.push_back(_width * _height);

    return out;
}

inline void AsyncRendererLF::start(size_t _numThreads, size_t _width, size_t _height, Func _cb) noexcept {



}

inline void AsyncRendererLF::stop() noexcept {

}

inline void AsyncRendererLF::restart() noexcept {

}

inline void AsyncRendererLF::work(size_t _index) noexcept {

}

#endif //ASYNCRENDERER_LOCKFREE_HPP