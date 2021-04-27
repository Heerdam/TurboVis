
#include "../../include/defines.hpp"
#include "../../include/data.hpp"

int main() {

    constexpr size_t n = 3;

    using Vector = Eigen::Matrix<float, n, 1>; //rows, cols
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
        std::uniform_real_distribution<float> distf(-1000.f, 1000.f);
        for(size_t i = 0; i < samples; ++i){
            Vector v;
            for(size_t d = 0; d < n; ++d)
                v[d] = distf(gen);
            signal[i] = v;
        }             

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
            const size_t cs = 10+c*5; //clusters
            auto start = std::chrono::high_resolution_clock::now();
            const auto pca = Data::PCA::PCA<Vector, decltype(signal.begin()), float>(signal.begin(), signal.end(), cs);
            auto end = std::chrono::high_resolution_clock::now() - start;

            //recover solution
            double result = 0.;
            for(size_t cl = 0; cl < pca->vq->clusters.size(); ++cl){
                const size_t m = pca->vq->clusters[cl].size();
                for(size_t p = 0; p < m; ++p){

                    Vector sol = pca->vq->means[cl];
                    for(size_t i = 0; i < n; ++i){
                        sol += pca->weights[cl](i, p) * pca->pca[cl].col(i);
                    }

                    //error
                    const double error = (signal[pca->vq->clusters[cl][p]] - sol).squaredNorm();
                    //std::cout << error << std::endl;
                    result += error;

                }
            }

            spdlog::info("Error: {}", std::sqrt(result));
            
        }

        

    }

    return 0;

}