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

        template<class Vector, size_t k = 3>
        class KD_TREE {

            struct Node {
                Vector* location;
                std::unique_ptr<Node> leftChild;
                std::unique_ptr<Node> rightChild;
            };

            std::unique_ptr<Node> root;

        public:

            template<class Iterator>
            void build(Iterator _begin, Iterator _end) {
                std::sort(_begin, _end, [&](const Vector& _v1, const Vector& _v2){
                    return _v1[k] > _v2[k];
                });
                root = insert(_begin, _end, 0);

            }         

            [[nodiscard]] Vector* nn(const Vector& _in) const noexcept {
                const float dist2 = (root->location - _in).squaredNorm();
                return nn(root.get()->location[0] >= _in[0] = root->leftChild : root->rightChild, root.get(), dist2, _in, 0);
            }

        private:
            template<class Iterator>
            [[nodiscard]] std::unique_ptr<Node> insert(Iterator _begin, Iterator _end, size_t _depth) noexcept {
                constexpr size_t axis = _depth % k;
                Iterator median = (_begin + (_end - _begin) / 2);
                std::unique_ptr<Node> node = std::make_unique<Node>();
                node->location = &(*median);
                node->leftChild = insert(_begin, median, _depth + 1);
                node->rightChild = insert(median, _end, _depth + 1);
                return node;
            }

            [[nodiscard]] Vector* nn(const Node* _node, const Node* best, float _dist2, const Vector& _in, size_t _depth) const noexcept {
                if(_node == nullptr) return best;
                constexpr size_t axis = _depth % k;
                const float dist2 = (_node->location - _in).squaredNorm();
                return nn(_node.get()->location[k] >= _in[k] = _node->leftChild : root->rightChild, dist2 < _dist2 ? _node : best, dist2, _in, _depth + 1);
            }
            
        };

        /*
            Vector Quantizer to create clusters from an input signal

        */
        template<class Signal>
        [[nodiscard]] inline std::vector<std::vector<Eigen::Matrix<float, Signal::RowsAtCompileTime, 1>>> VQ(const Signal& _signal, size_t _clusters) {

            using Vector = Eigen::Matrix<float, _signal::RowsAtCompileTime, 1>;

            std::mt19937_64 gen;
            std::uniform_int_distribution<size_t> dist (0, _signal.cols() - 1);

            std::vector<bool> map(_signal.cols());
            std::fill(map.begin(), map.end(), false);

            std::vector<Vector> words (_signal.cols());

            std::unordered_map<Vector*, size_t> lookUpTable;

            for(size_t i = 0; i <_clusters; ++i){
                const size_t c = dist(gen);
                if(map[c]){
                    --i;
                    continue;
                }
                map[c] = true;
                Vector v = _signal.col(c);
                words.push_back(v);
                lookUpTable[&v] = i; //falsch
            }

            KD_TREE<Vector, _signal::RowsAtCompileTime> tree;
            tree.insert(words);

            std::vector<std::vector<Vector>> out (_clusters);

            for(size_t i = 0; i < _signal.cols(); ++i){
                if(!map[i]){

                    auto word = tree.nn(_signal.row(i));
                    assert(lookUpTable.count(word) != 0); //falsch
                    size_t index = lookUpTable[word]; //falsch
                    out[index].push_back(words[index]); //semnatisch falsch will words sortiert wird
                }
            }
            

            return out;
            
        };

        inline void staticPCA(const Eigen::Matrix<float, 3, 3>& _signal, size_t _clusters){

            //first pick of vectors for VQ
            

            Eigen::Matrix<float, 3, 3> words = _signal(Eigen::all, clusterIndicies);
            Eigen::Vector2f f;
            

            //calc first vq
            auto clusters = VQ(_signal, clusterIndicies);

            //calc mean of cluster

            //iterativly use means for VQ until no points are switched or iteration count is reached

            /*
                for every cluster
                    calc mean x0
                    compute residual matrix C = X - x0
                    SVD(C)
                    PCA vecs: first n rows of V
                    projection weight of point j: first n columns of row j of UD
            */
        };

        /*
        performs a principal component analysis utlising the kernel Hebbian algorithm with stochastic meta-descent
        From: Fast Iterative Kernel Principal Component Analysis, 2007, Günter el al.
        REM: input matrix
        iterations: #of steps it should do
        eta_0, mu and xsi: tuning parameter
        */
        template <class M>
        [[nodiscard]] inline M KHA_SMD(const M& _REM, size_t _iterations, float _etc_0 = 0.1f, float _mu = 1.f, float _xsi = 0.99f) {
            constexpr size_t n = _REM::RowsAtCompileTime;
            constexpr size_t m = _REM::ColsAtCompileTime;

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