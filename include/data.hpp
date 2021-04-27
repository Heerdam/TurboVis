#ifndef DATA_HPP
#define DATA_HPP

#include "defines.hpp"

namespace Data {

    //deprecated
    using Triplet = std::tuple<float, uint16_t, uint16_t, uint16_t>;

    namespace PCA {

        namespace detail {

            template<class Vector, class Iterator, class T>
            [[nodiscard]] size_t nn(Iterator _begin, const Iterator _end, const Iterator _signal, const Vector& _e) noexcept {
                size_t idx = 0;
                T dist2_best = std::numeric_limits<T>::max();
                for (size_t i = 0; _begin != _end; ++_begin, ++i) {
                    const Vector& w = *_begin;
                    const T dist2 = (w - _e).squaredNorm();
                    if (dist2 < dist2_best) {
                        dist2_best = dist2;
                        idx = i;
                    }
                }
                return idx;
            };

            template <class Vector>
            struct VQResult {
                const std::vector<Vector> means;
                const std::vector<std::vector<size_t>> clusters;
                const std::vector<size_t> words;
                const size_t iterations;
                const std::vector<double> time;

                VQResult(
                    std::vector<Vector>&& _means, 
                    std::vector<std::vector<size_t>>&& _clusters,
                    std::vector<size_t>&& _words,
                    size_t _iterations,
                    std::vector<double>&& _time ) : 
                    iterations(_iterations), 
                    means(std::forward<std::vector<Vector>>(_means)), 
                    clusters(std::forward<std::vector<std::vector<size_t>>>(_clusters)), 
                    words(std::forward<std::vector<size_t>>(_words)), 
                    time(std::forward<std::vector<double>>(_time)) {}
            };

            template <class Vector, class Matrix_PCA, class Matrix_Weights>
            struct PCAResult {
                const std::unique_ptr<VQResult<Vector>> vq;
                const std::vector<Matrix_PCA> pca;
                const std::vector<Matrix_Weights> weights;
                const size_t time;

                PCAResult(std::unique_ptr<VQResult<Vector>>&& _vq, std::vector<Matrix_PCA>&& _pca, std::vector<Matrix_Weights>&& _weights, size_t _time) : 
                    vq(std::forward<std::unique_ptr<VQResult<Vector>>>(_vq)), 
                    pca(std::forward<std::vector<Matrix_PCA>>(_pca)), 
                    weights(std::forward<std::vector<Matrix_Weights>>(_weights)), 
                    time(_time) {}
            };

        }

        /*
            Vector Quantizer to create clusters from an input signal
            iterativly use means for VQ until no points are switched or iteration count is reached
            compute means of clusters
        */
        template<class Vector, class Iterator, class T>
        [[nodiscard]] std::unique_ptr<detail::VQResult<Vector>> VQ(Iterator _begin, Iterator _end, size_t _ccount, size_t _iterations = 100, T _epsilon2 = T(0.5), size_t _threads = 4) {

            const size_t ssize = _end - _begin;
            spdlog::debug("[VQ] Started. Signal size: {}, clusters: {}, max iterations: {}, epsilon2: {}, max threads: {}", ssize, _ccount, _iterations, _epsilon2, _threads);

            // ------------------- OUT -------------------
            std::vector<Vector> means (_ccount);
            std::vector<std::vector<size_t>> clusters (_ccount);
            std::vector<size_t> words (_ccount);
            size_t iterations = 0;
            std::vector<double> time;

            // ------------------- TEMP -------------------
            std::vector<Vector> oldMeans (_ccount);
            std::vector<size_t> count (_ccount);

            for(size_t it = 0; it < _iterations; ++it){
                auto now = std::chrono::high_resolution_clock::now();
                spdlog::debug("[VQ] Iteration {} start.", it);

                for(auto& v : clusters)
                    v.clear();

                //create words
                if(it == 0){
                    std::mt19937_64 gen;
                    gen.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
                    std::uniform_int_distribution<size_t> dist(0, ssize - 1);

                    std::vector<bool> map(ssize);
                    std::fill(map.begin(), map.end(), false);

                    for (size_t i = 0; i < _ccount; ++i) {
                        const size_t c = dist(gen);
                        //avoid inserting twice the same
                        if (map[c]) {
                            --i;
                            continue;
                        }
                        map[c] = true;
                        words[i] = c;
                    }

                    for(size_t i = 0; i < _ccount; ++i)
                        oldMeans[i] = *(_begin + words[i]);

                } else {

                    
                    for(size_t i = 0; i < _ccount; ++i){;
                        oldMeans[i] = means[i];
                        means[i].setZero();
                    }
                    

                }

                //insert signal into the clusters (including words)
                std::memset(count.data(), 0, _ccount * sizeof(size_t));

                for(Eigen::Index i = 0; i < ssize; ++i){
                    //bool isSet = false;                
                    const size_t idx = detail::nn<Vector, Iterator, T>(oldMeans.begin(), oldMeans.end(), _begin, *(_begin + i)); //cluster idx
                    clusters[idx].push_back(i);
                    means[idx] = means[idx] + (*(_begin + i) - means[idx]) / ++count[idx];
                }

                bool done = true;
                for(size_t i = 0; i < _ccount; ++i){
                    //std::cout << means[i] << std::endl << std::endl;
                    //std::cout << words[i] << std::endl;
                    const T dist2 = (means[i] - oldMeans[i]).squaredNorm();
                    //std::cout << dist2 << std::endl << std::endl << std::endl;
                     if(dist2 > _epsilon2)
                        done = false;
                }
                const std::chrono::duration<double, std::milli> dtime = (std::chrono::high_resolution_clock::now() - now);
                spdlog::debug("[VQ] Iteration {} ended. ({} ms)", it, dtime.count());
                time.push_back(dtime.count());

                if(done) break;

            }
            return std::make_unique<detail::VQResult<Vector>>(std::move(means), std::move(clusters), std::move(words), iterations, std::move(time));
        };

        /*
            for every cluster
            calc mean x0
            compute residual matrix C = X - x0
            SVD(C)
            PCA vecs: first n rows of V
            projection weight of point j: first n columns of row j of UD
        */
        template <class Vector, class Iterator, class T>
        [[nodiscard]] auto PCA(Iterator _begin, Iterator _end, size_t _ccount){
            auto now = std::chrono::high_resolution_clock::now();

            constexpr size_t n = Vector::RowsAtCompileTime;

            using Matrix_MxN = Eigen::Matrix<T, -1, n>;
            using Matrix_NxM = Eigen::Matrix<T, n, -1>;

            auto vq = VQ<Vector, Iterator, T>(_begin, _end, _ccount);

            std::vector<Matrix_MxN> pca (_ccount);
            std::vector<Matrix_NxM> weights (_ccount);

            for(size_t k = 0; k < vq->clusters.size(); ++k){ 
                const auto& cluster = vq->clusters[k];
                const size_t m = cluster.size();

                Matrix_NxM C (n, m); //residiual matrix m x n

                for(size_t i = 0; i < m; ++i)   
                    C.col(i) = *(_begin + cluster[i]) - vq->means[k];

                Matrix_MxN CC (m, n);
                CC = C.transpose();

                //std::cout << CC << std::endl;

                Eigen::BDCSVD<Matrix_MxN> solver (CC, Eigen::ComputeFullU | Eigen::ComputeFullV);
                pca[k] = solver.matrixV(); //n x n
                
                Matrix_MxN D (m, n);
                D.setZero();
                const auto SV = solver.singularValues();
                //std::cout << SV.cols() << " " << SV.rows() << std::endl;
                for(Eigen::Index i = 0; i < SV.rows(); ++i)
                    D(i, i) = SV(i);
            
                const auto UD = (solver.matrixU() * D).transpose();
                //std::cout << UD.cols() << " " << UD.rows() << std::endl;
                weights[k] = UD;

            }          
            return std::make_unique<detail::PCAResult<Vector, Matrix_MxN, Matrix_NxM>>(std::move(vq), std::move(pca), std::move(weights), size_t((std::chrono::high_resolution_clock::now() - now).count()));
        };

        
        //template<class Vector>
        //[[nodiscard]] std::pair<std::vector<Vector> /*means*/, std::vector<std::vector<Vector>> /*clusters*/> VQ(const std::vector<Vector>& _signal, size_t _clusters, size_t& _iterations) {
/*
            std::unordered_map<std::uintptr_t , size_t> clusterMap(_signal.size());
            std::vector<std::vector<Vector>> out (_clusters);
            std::vector<Vector> words (_clusters);
            std::vector<Vector> means (_clusters);

            for(size_t it = 0; it < 100; ++it){

                for(auto& v : out)
                    v.clear();
                bool hasChanged = false;

                //create words
                if(it == 0){
                    std::mt19937_64 gen;
                    gen.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
                    std::uniform_int_distribution<size_t> dist(0, _signal.size() - 1);

                    std::vector<bool> map(_signal.size());
                    std::fill(map.begin(), map.end(), false);

                    for (size_t i = 0; i < _clusters; ++i) {
                        const size_t c = dist(gen);
                        //avoid inserting twice the same
                        if (map[c]) {
                            --i;
                            continue;
                        }
                        map[c] = true;
                        Vector v = _signal[c];
                        words[i] = v;
                    }

                } else {

                    for(size_t i = 0; i < _clusters; ++i){
                        //std::cout << means[i] << std::endl;
                        words[i] = means[i];
                        means[i].setZero();
                    }

                    //std::cout << std::endl;

                }

                //insert signal into the clusters (including words)
                std::vector<size_t> count (_clusters);
                std::memset(count.data(), 0, _clusters * sizeof(size_t));
                for(Eigen::Index i = 0; i < _signal.size(); ++i){
                    bool isSet = false;                
                    const size_t idx = nn(words, _signal[i]); //cluster idx
                    out[idx].push_back(_signal[i]);
                    means[idx] = means[idx] + (_signal[i] - means[idx]) / ++count[idx];
                }

                bool done = true;
                    for(size_t i = 0; i < _clusters; ++i){
                        const float dist2 = (means[i] - words[i]).squaredNorm();
                        if(dist2 > 0.5f)
                            done = false;
                    }
                    if(done) {
                        _iterations = it;
                        break;
                    }

            }
 
            return { means, out };
            
        };
*/

       
        /*
        performs a principal component analysis utlising the kernel Hebbian algorithm with stochastic meta-descent
        From: Fast Iterative Kernel Principal Component Analysis, 2007, Günter el al.
        REM: input matrix
        iterations: #of steps it should do
        eta_0, mu and xsi: tuning parameter
        */
        template <class M>
        [[nodiscard]] M KHA_SMD(const M& _REM, size_t _iterations, float _etc_0 = 0.1f, float _mu = 1.f, float _xsi = 0.99f) {
            constexpr size_t n = _REM::RowsAtCompileTime;
            constexpr size_t m = _REM::ColsAtCompileTime; //TODO: kritisch

            using Vector = Eigen::Matrix<float, 1, m>;
            using Matrix = Eigen::Matrix<float, n, m>;

            //initialization
            const float eta_0 = 1.f;  //tuning parameter
            const size_t l = n;       //observations

            const Eigen::Matrix<float, n, n> M(1.f / n);
            const Matrix MK = M * _REM;
            const Matrix K = _REM - MK - MK.transpose() - MK * M;  //K'

            //a
            Matrix MK = M * K;
            Matrix MKM = MK * M;

            //b
            std::random_device rd;
            std::mt19937_64 gen(rd());
            std::normal_distribution<float> dist(1.f, 0.1f);
            Matrix A;  //A_1
            for (size_t i = 0; i < n; ++i)
                for (size_t j = 0; j < m; ++j)
                    A(i, j) = dist(gen);

            //c
            Matrix AK = A * K;  //A_1_K

            //d
            Vector rho;  //ro_0
            rho.setOnes();

            Matrix B;  // B_1
            B.setZero();

            //loop
            for (size_t t = 0; t < _iterations; ++t) {  //t: iteration number, g(t) = floor(n/_steps*t)

                const size_t i = static_cast<size_t>(std::floor(n / _steps * t));
                //a
                const auto d1 = (A * K * (A * K).transpose());
                const auto d2 = (A * A.transpose());
                const auto S1 = d1.diagonal();
                const auto S2 = d2.diagonal();

                Vector delta;  //update every l iterations
                for (size_t i = 0; i < n; ++i)
                    delta[i] = std::sqrt(S1[i] / S2[i]);

                //b
                const float delta_norm = delta.norm();
                Vector eta;  //gain vector - update every l iterations
                for (size_t i = 0; i < n; ++i)
                    eta[i] = (delta_norm / delta[i]) * (l / t + l) * eta_0;

                //c
                const Vector k_i = K.row(i);  //k'i = ith row of K
                const Vector k_i_t = k_i.transpose();

                //d
                const Vector y_t = A * k_i;
                const Vector y_t_t = y_t.transpose();

                //e
                const Matrix yt_yt = y_t * y_t_t.triangularView<Eigen::Lower>();
                const Matrix yt_yt_A = yt_yt * A;
                const Matrix gamma = y_t * Matrix::Identity().transpose() - yt_yt_A;

                //f
                Matrix gamma_K = y_t * k_i_t - yt_yt_A * K;

                //g
                const Matrix gamma_k_B = gamma * K * B;
                rho = rho + mu * gamma_k_B.diagonal();
                Vector e_rho;
                for (std::ptrdiff_t i = 0; i < e_rho.cols(); ++i)
                    e_rho[i] = std::exp(rho[i]);
                const Matrix e_rho_m = Matrix::Identity() * e_rho;

                //h
                const Matrix xsi_B = xsi * B;
                const auto inner = (B * k_i * y_t_t + y_t * k_i_t * B.transpose()).triangularView<Eigen::Lower>().toDenseMatrix();
                B = xsi_B + e_rho_m * gamma.diagonal() * ((A + xsi_B) * k_i * Matrix::Identity().row(i).transpose() - yt_yt * (A + xsi_B) - xsi * inner) * A;

                //i
                AK = A * K + e_rho * gamma.diagonal() * gamma * K; 

                //j
                A = A + e_rho * eta.diagonal() * gamma;
            }
            return A;
        }
    }  // namespace PCA

}


#endif /* DATA_HPP */
