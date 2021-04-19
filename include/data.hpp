#ifndef DATA_HPP
#define DATA_HPP

#include "defines.hpp"

namespace Data {
    using Triplet = std::tuple<float, uint16_t, uint16_t, uint16_t>;
    using Iterator = std::vector<Triplet>::iterator;

    enum class Gridsize : uint16_t {
        x2      =   2, 
        x4      =   4, 
        x16     =   16, 
        x32     =   32, 
        x64     =   64, 
        x128    =   128,
        x256    =   256, 
        x512    =   512, 
        x1024   =   1024, 
        x2048   =   2048, 
        x4096   =   4096
    };

    template<class T>
    inline constexpr size_t operator*(T _v1, Gridsize _v2){
        return _v1 * static_cast<size_t>(_v2);
    }

    template<class T>
    inline constexpr size_t operator*(Gridsize _v2, T _v1){
        return _v1 * static_cast<size_t>(_v2);
    }

    inline constexpr size_t operator*(Gridsize _v1, Gridsize _v2){
        return static_cast<size_t>(_v1) * static_cast<size_t>(_v2);
    }

    class SparseGrid {

        std::mutex mutex;
        std::unordered_map<size_t, float> map;
        Gridsize size;

    public:
        SparseGrid();
        void init(Gridsize);
        void downsample(Gridsize);
        Gridsize getSize() const;
        size_t getSizeOfMapInBytes() const;

        void insert(float, uint16_t, uint16_t, uint16_t);
        void insertBulk(Iterator, Iterator);

        void insert_async(float, uint16_t, uint16_t, uint16_t);
        void insertBulk_async(Iterator, Iterator);

        std::vector<float> toBuffer() const;
        std::vector<float> createBuffers(float) const;
    };

    class FunctionData {

        std::mutex mutex;
        SparseGrid grid;

        static void eval(std::function<float(int16_t, int16_t, int16_t)>, SparseGrid&, bool, int16_t, int16_t, int16_t, int16_t, int16_t, int16_t);
        
    public:

        void createDataFromFunction(std::function<float(int16_t, int16_t, int16_t)>, Gridsize);
        void downsample(Gridsize);

        std::vector<float> toBuffer() const;
        std::vector<float> createBuffers() const;

        size_t getDataSizeinBytes() const;
        Gridsize currentGridSize() const;

    };

    namespace PCA {

        template<class Vector>
        [[nodiscard]] size_t nn(const std::vector<Vector>& _vals, const Vector& _e){
            size_t idx = 0;
            float dist2_best = std::numeric_limits<float>::max();
            for(size_t i = 0; i < _vals.size(); ++i){
                const float dist2 = (_vals[i] - _e).squaredNorm();
                if(dist2 < dist2_best){
                    dist2_best = dist2;
                    idx = i;
                }
            }
            return idx;
        };

        /*
            Vector Quantizer to create clusters from an input signal
            iterativly use means for VQ until no points are switched or iteration count is reached
            compute means of clusters
        */
        template<class Vector>
        [[nodiscard]] std::pair<std::vector<Vector>, std::vector<std::vector<Vector>>> VQ(const std::vector<Vector>& _signal, size_t _clusters, size_t& _iterations) {

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

        /*
            for every cluster
            calc mean x0
            compute residual matrix C = X - x0
            SVD(C)
            PCA vecs: first n rows of V
            projection weight of point j: first n columns of row j of UD
        */
        template <class Signal>
        std::pair<std::vector<Eigen::Matrix<float, -1, Signal::RowsAtCompileTime>>, std::vector<Eigen::Matrix<float, -1, Signal::RowsAtCompileTime>>> PCA(const Signal& _signal, size_t _clusters){

            constexpr size_t n = Signal::RowsAtCompileTime;

            using Matrix = Eigen::Matrix<float, -1, Signal::RowsAtCompileTime>;
            using Vector = Eigen::Matrix<float, 1, n>;

            const auto clusters = VQ(_signal, _clusters);

            std::vector<Matrix> pca (_clusters);
            std::vector<Matrix> weights (_clusters);

            for(size_t k = 0; k < clusters.size(); ++k){ 
                const auto& c = clusters[k];
                constexpr size_t m = c.size();
                Eigen::Matrix<float, m, n> C; //residiual matrix m x n
                for(size_t i = 0; i < m; ++i)   
                    C.col(i) = c[i] - x_0;

                Eigen::JacobiSVD<decltype(C)> solver (C, Eigen::ComputeThinU | Eigen::ComputeThinV);
                pca[k] = solver.matrixV().transpose().block(0, m, 0, n);
                weights[k] = solver.matrixU() * solver.singularValues().block(0, m, 0, n);
            }

            return { pca, weights };

        };

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
