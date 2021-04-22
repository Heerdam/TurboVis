#ifndef BENCH_HPP
#define BENCH_HPP
#include <chrono>
#include <numeric>
#include <utility>
namespace std{
        unsigned long long nanoTime(){
                using namespace std;
                using namespace chrono;
                return duration_cast<nanoseconds>(high_resolution_clock::now().time_since_epoch()).count();
        }
}
template<size_t run_length = 300, typename Function, typename... Ts>
double bench_function(Function f, Ts&&... args){
        static_assert(run_length > 0, "run_length must be positive");
        unsigned int since_best = 0;
        unsigned long long best_time = std::numeric_limits<unsigned long long>::max();
        unsigned long long accum = 0;
        while(true){
                const auto t1 = std::nanoTime();
                f(std::forward<Ts>(args)...);
                const auto t2 = std::nanoTime();
                if(t2 - t1 < best_time){
                        best_time = t2 - t1;
                        since_best = 0;
                        accum = 0;
                }
                else{
                        since_best++;
                }
                accum += (t2 - t1);
                if(since_best >= run_length){
                        break;
                }
        }
        return accum / double(run_length);
}

#endif
