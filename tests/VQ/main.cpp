
#include "../../include/data.hpp"

#include <sciplot/sciplot.hpp>

#include <iostream>

using namespace sciplot;

template <class Vector>
void toPlot(Plot& _plot, const std::vector<std::vector<Vector>>& data) {

    std::string cols[] = {
        "black", "dark-grey", "red", "web-green", "web-blue", "dark-magenta", "dark-cyan", "dark-orange", "dark-yellow", "royalblue", 
        "goldenrod", "dark-spring-green", "purple", "steelblue", "dark-red", "dark-chartreuse", "orchid", "aquamarine", "brown", "yellow", 
        "turquoise", "grey", "light-red", "light-green", "light-blue", "light-magenta", "light-cyan", "light-goldenrod", "light-pink", "light-turquoise", 
        "gold", "green", "dark-green", "spring-green", "forest-green", "sea-green", "blue", "dark-blue", "midnight-blue", "navy", 
        "medium-blue", "skyblue", "cyan", "magenta", "dark-turquoise", "dark-pink", "coral", "light-coral",
        "dark-khaki", "dark-goldenrod", "beige", "olive", "orange", "violet" 
    };


    _plot.legend().hide();

    for (size_t k = 0; k < data.size(); ++k) {
        const auto& cl = data[k];
        std::vector<float> x(cl.size());
        std::vector<float> y(cl.size());
        for (size_t l = 0; l < cl.size(); ++l) {
            x[l] = cl[l](0, 0);
            y[l] = cl[l](1, 0);
        }
        _plot.drawPoints(x, y).lineColor(cols[k]);
    }
    _plot.xlabel("");
    _plot.ylabel("");

    std::stringstream title;
    title << "set title \"" << "Clusters: " << data.size() << ", Size: 1000" << "\" offset 0,-2 font \",25\"";
    _plot.gnuplot(title.str());
    _plot.gnuplot("show title");

    _plot.autoclean(true);
};

template <class Vector>
void plot(
        std::vector<std::vector<float>> _time, 
        std::vector<std::vector<size_t>> _iterations,
        std::vector<std::vector<std::string>> _testName,
        std::vector<std::vector<std::vector<std::vector<Vector>>>> _clusters)
    {

    const size_t n = _time.size();



    for(size_t i = 0; i < 1; ++i){

        //-------------------- TIME --------------------
        Plot p_time;
        {
            const auto& time = _time[i];
            const auto& names = _testName[i];
            std::vector<float> xaxis(_testName[0].size());

            for(size_t i = 0; i < xaxis.size(); ++i)
                xaxis[i] = 0.5f;

            p_time.legend().hide();
            p_time.ylabel("ms");
            p_time.xlabel("clusters");
            p_time.xtics().rotate(true);
            p_time.drawBoxes(names, time, xaxis)
                .fillSolid()
                .fillColor("coral")
                .fillIntensity(0.9)
                .borderShow()            
                .labelNone();

            std::stringstream title;
            title << "set title \"runtime in ms\" offset 0,-3 font \",25\"";
            p_time.gnuplot(title.str());
            p_time.gnuplot("show title");

            p_time.boxWidthRelative(0.5f);
            p_time.autoclean(true);
        }

        //-------------------- ITERATIONS --------------------
        Plot p_it;
        {
            const auto& iterations = _iterations[i];
            const auto& names = _testName[i];
            std::vector<float> xaxis(_testName[0].size());

            for(size_t i = 0; i < xaxis.size(); ++i)
                xaxis[i] = 0.5f;

            p_it.legend().hide();
            p_it.ylabel("iterations");
            p_time.xlabel("clusters");
            p_it.xtics().rotate(true);
            p_it.drawBoxes(names, iterations, xaxis)
                .fillSolid()
                .fillColor("coral")
                .fillIntensity(0.9)
                .borderShow()            
                .labelNone();

            std::stringstream title;
            title << "set title \"iterations for clusters to convergence\" offset 0,-3 font \",25\"";
            p_it.gnuplot(title.str());
            p_it.gnuplot("show title");

            p_it.boxWidthRelative(0.5f);
            p_it.autoclean(true);
        }

        //-------------------- CLUSTERS --------------------
        Plot p_d_1;
        toPlot(p_d_1, _clusters[i][0]);

        Plot p_d_2;
        toPlot(p_d_2, _clusters[i][1]);

        Plot p_d_3;
        toPlot(p_d_3, _clusters[i][2]);

        Plot p_d_4;
        toPlot(p_d_4, _clusters[i][3]);

        Plot p_d_5;
        toPlot(p_d_5, _clusters[i][4]);

        Plot p_d_6;
        toPlot(p_d_6, _clusters[i][5]);

        Plot p_d_7;
        toPlot(p_d_7, _clusters[i][6]);

        Plot p_d_8;
        toPlot(p_d_8, _clusters[i][7]);

        Plot p_d_9;
        toPlot(p_d_9, _clusters[i][8]);

        Plot p_d_10;
        toPlot(p_d_10, _clusters[i][9]);

        //-------------------- FIGURE --------------------
        Figure fig = 
        {
            { p_time, p_it },
            { p_d_1, p_d_2 },
            { p_d_3, p_d_4 },
            { p_d_5, p_d_6 },
            { p_d_7, p_d_8 },
            { p_d_9, p_d_10 }
        };

        fig.size(2*1000, 6*1000);
        fig.fontSize(16);        
       
        std::stringstream ss;
        ss << "vq_" << i + 1 << ".png";
        fig.save(ss.str());

    }
}

/*

*/
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

        constexpr size_t tests = 10;

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
            const size_t cs = 5+c*5;
            auto start = std::chrono::high_resolution_clock::now();
            size_t it;
            const auto res = Data::PCA::VQ(signal, cs, it);
            auto end = std::chrono::high_resolution_clock::now() - start;
            
            time[0].push_back(static_cast<float>(std::chrono::duration_cast<std::chrono::milliseconds>(end).count()) / 1000.f);
            testNames[0].push_back(std::to_string(cs));         
            iterations[0].push_back(it);

            clusters[0][c].resize(res.second.size()); //# clusters
            for(size_t i = 0; i < res.second.size(); ++i){      
                clusters[0][c][i].resize(res.second[i].size());        
                const size_t s = res.second[i].size() * sizeof(Vector);
                std::memcpy(clusters[0][c][i].data(), res.second[i].data(), s);             
            }
        }

        plot(time, iterations, testNames, clusters);

    }

    return 0;
}