
#include "../../include/data.hpp"

#include <sciplot/sciplot.hpp>

#include <iostream>

template<class Vector>
void plot(
        std::vector<std::vector<float>> _time, 
        std::vector<std::vector<std::string>> _testName,
        std::vector<std::vector<std::vector<Vector>>> _clusters)
    {

    using namespace sciplot;

    //assert(_time.size() == _testName.size() == _clusters.size());
    const size_t n = _time.size();


    for(size_t i = 0; i < 1; ++i){

        //plot time to convergence for method and clusters
        const auto& time = _time[i];
        const auto& names = _testName[i];
        std::vector<float> xaxis(_testName[0].size());

        for(size_t i = 0; i < xaxis.size(); ++i)
            xaxis[i] = 0.5f;

        Plot plot;
        plot.size(1000, 1000);
        plot.fontSize(16);
        plot.legend().hide();
        plot.ylabel("ms");
        plot.xtics().rotate(true);
        plot.yrange(0, 150);
        plot.drawBoxes(names, time, xaxis)
            .fillSolid()
            .fillColor("blue")
            .fillIntensity(0.9)
            .borderShow()            
            .labelNone();

        plot.boxWidthRelative(0.5f);

        plot.autoclean(true);

        plot.save("example-boxes-ticklabels.png");

        //plot points with colors for clusters

    }
}

/*

*/
int main() {

    using Vector = Eigen::Matrix<float, 3, 1>; //rows, cols
    using namespace std::chrono_literals;

    spdlog::set_level(spdlog::level::trace);
    

    std::mt19937_64 gen;
    gen.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    //uniform distributed signal 
    {
        const std::string testName = "small normal distributed cluster";
        spdlog::info("Test: {}", testName);
        std::vector <Vector> signal(1000);
        std::uniform_real_distribution<float> distf(-1000.f, 1000.f);
        for(size_t i = 0; i < signal.size(); ++i){
            signal[i] = { distf(gen), distf(gen), distf(gen) };
        }

        std::vector<std::vector<float>> time (2);
        time[0].reserve(10);
        std::vector<std::vector<std::string>> testNames (2);
        testNames[0].reserve(10);
        std::vector<std::vector<std::vector<Vector>>> clusters (2);
        clusters[0].reserve(10);

        for(size_t c = 0; c < 10; ++c){
            auto start = std::chrono::high_resolution_clock::now();
            const auto res = Data::PCA::VQ(signal, 5+c*5);
            auto end = std::chrono::high_resolution_clock::now() - start;
            time[0][c] = static_cast<float>(std::chrono::duration_cast<std::chrono::milliseconds>(end).count()) / 1000.f;
            clusters[0][c].reserve(res.second.size());
            std::memcpy(clusters[0][c].data(), res.second.data(), res.second.size() * sizeof(Vector));
            testNames[0][c] = std::to_string(c);
        }

        plot(time, testNames, clusters);

    }

    return 0;

    //normal distributed signal
    {

    }

   




        /* 
        //# of clusters
        std::uniform_int_distribution<size_t> dist(50, 100);
        const size_t clusterCount = dist(gen);

        //create clusters
        size_t count = 0;
        std::vector<std::vector<Vector>> clusters (clusterCount);
        std::vector<bool> cmap(clusterCount*clusterCount);
        std::fill(cmap.begin(), cmap.end(), false);
        dist = std::uniform_int_distribution<size_t>(0, clusterCount-1);
        std::uniform_real_distribution<float> distf(-10.f, 10.f);
        for(size_t j = 0; j < clusterCount; ++j){ //assume grid with size clusterCount x clusterCount

            for(size_t k = 0; k < 25; ++k){
                const size_t x = dist(gen);
                const size_t y = dist(gen);
                if(cmap[y * clusterCount + x]){
                    --j;
                    continue;
                }
                Vector v = Vector(distf(gen) + x * 100.f, distf(gen) + y * 100.f, distf(gen));
                clusters[j].push_back(v);
                ++count;
            }
        }

        //create random signal from clusters.
        Eigen::Matrix<float, 3, -1> signal; //rows, cols
        {
            std::vector<Vector> lin (count);
           
            for(size_t j = 0, ptr = 0; j < clusterCount; ++j){
                std::memcpy(lin.data() + ptr, clusters[j].data(), clusters[j].size() * sizeof(Vector));
                ptr += clusters[j].size();
            }

            std::shuffle(lin.begin(), lin.end(), gen);

            for(size_t j = 0; j < count; ++j)
                signal << lin[j];
        }

        //call VQ
        const auto vq = Data::PCA::VQ(signal, clusterCount, 25);

        auto findMatch = [](const std::vector<Vector>& _data, const Vector& _toFind)-> size_t {
            for(const auto& v : _data)
                if(v.isApprox(_toFind)) return 1;
            return 0;
        };

        //compare result with initial clusters
        for(size_t j = 0; j < clusterCount; ++j){

            const auto& vals_out = vq.second[j];
            const auto& vals_expt = clusters[j];
            size_t found = 0;
            for(const auto& v : vals_expt)
                found += findMatch(vals_out, v);

            spdlog::debug("Cluster {}: {} / {} found", j, found, vals_expt.size());

        }

    }
    */
    

    return 0;
}