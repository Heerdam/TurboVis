
#include "../../include/defines.hpp"
#include "../../include/data.hpp"

int main() {

    using Vector = Eigen::Matrix<float, 2, 1>; //rows, cols
    using namespace std::chrono_literals;

    spdlog::set_level(spdlog::level::trace);
    

    std::mt19937_64 gen;
    gen.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    //uniform distributed signal 
    {
        const std::string testName = "small normal distributed cluster";
        spdlog::info("Test: {}", testName);

        constexpr size_t samples = 1000;

        //create the signal
        std::vector <Vector> signal(samples);
        std::uniform_real_distribution<float> distf(-100.f, 100.f);
        for(size_t i = 0; i < samples; ++i)
            signal[i] = { distf(gen), distf(gen) };

        constexpr size_t tests = 1;

        std::vector<std::vector<float>> time (2);
        time[0].reserve(tests);
        time[1].resize(tests);
        std::vector<std::vector<size_t>> iterations(2);
        iterations[0].reserve(tests);
        iterations[1].reserve(tests);
        std::vector<std::vector<std::string>> testNames (2);
        testNames[0].reserve(tests);
        testNames[1].reserve(tests);
        std::vector<std::vector<std::vector<std::vector<Vector>>>> clusters (2); // meta - tests - #clusters - clusters
        clusters[0].resize(tests);
        clusters[1].resize(tests);

        for(size_t c = 0; c < tests; ++c){
            const size_t cs = 5+c*5; //clusters
            auto start = std::chrono::high_resolution_clock::now();
            const auto pca = Data::PCA::PCA<Vector, decltype(signal.begin()), float>(signal.begin(), signal.end(), cs);
            auto end = std::chrono::high_resolution_clock::now() - start;
            
        }

        

    }

    return 0;

}