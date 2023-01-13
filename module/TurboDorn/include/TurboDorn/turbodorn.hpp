#ifndef TURBODORN_HPP
#define TURBODORN_HPP

#include <complex>
#include <numeric>
#include <vector>
#include <array>
#include <filesystem>
#include <fstream>
#include <cassert>
#include <cstring>
#include <mutex>
#include <cmath>
#include <chrono>

#include <omp.h>

#include <Eigen/Eigen>
#include <tsl/robin_map.h>

#include <highfive/H5File.hpp>
#include <highfive/h5easy_bits/H5Easy_Eigen.hpp>

#include <nlohmann/json.hpp>

#include <libmorton/morton.h>

#include <lodepng.h>

using namespace nlohmann;
using namespace std::chrono_literals;

namespace TurboDorn {

    namespace Hagedorn::Detail {
        template<class T>
        EIGEN_STRONG_INLINE uint_fast64_t morton_index(const Eigen::Vector3<T>& _x);
    }

    namespace Geometry {

        namespace Detail {

            //---------------------------------------------------------------------------------------//
            //                                        Chunk
            //---------------------------------------------------------------------------------------//

            /*
                A chunk holds a certain number of values of the grid. Every chunk has a file on the disk
                associated and can be loaded and unloaded on demand.
            */

            template<class T, bool IS_SAMPLER>
            class Chunk {

                const std::filesystem::path path;

                tsl::robin_map<uint_fast64_t, std::complex<T>> data;

                const Eigen::Vector3i cardinal;
                const Eigen::Vector3<T> cell_size;
                const Eigen::Vector3<T> cell_half;

                bool is_loaded = false;

            public:
                /*
                  _file - the path to the chunk data of this chunk
                  _cardinal - the 3 main dimensions
                  _cell_size - the size of a single cell of the lattice in world unit
                */
                Chunk(std::filesystem::path _file, const Eigen::Vector3i& _cardinal, const Eigen::Vector3<T>& _cell_size) noexcept 
                    : path(std::move(_file)), cardinal(_cardinal), cell_size(_cell_size), cell_half(_cell_size * T(0.5)) { 

#ifdef DEBUG_TRACE
                    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
                    std::cout << "Chunk:\t\t\t\tNew \nPath:\t\t\t\t" << path << "\nCardinal:\t\t\t[" << 
                    _cardinal(0) << ", " << _cardinal(1) << ", " << _cardinal(2) << 
                    "]\nCell Size:\t\t\t[" << cell_size(0) << ", " << cell_size(1) << ", " << cell_size(2) << "]\n";
                    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";  
#endif
                }

                bool load() {
                    is_loaded = true;
                    data = tsl::robin_map<uint_fast64_t, std::complex<T>>();
                    if(!std::filesystem::exists(path)){
#ifdef DEBUG_TRACE
                        std::cout << "Chunk:\t\t\t\tLoad\nStatus:\t\t\t\tNot found" << "\nPath:\t\t\t\t" << path <<
                            "\nCount:\t\t\t\t" << data.size() << std::endl << std::endl;
#endif
                        return true;
                    } 
                    std::fstream in (path, in.binary | in.in);
                    if(in.good()){
                        size_t s = 0;
                        in.read(reinterpret_cast<char*>(&s), sizeof(s));
                        data.reserve(s);
                        for(size_t i = 0; i < s; ++i){
                            uint_fast64_t idx;
                            std::complex<T> val; 
                            in.read(reinterpret_cast<char*>(&idx), sizeof(idx));
                            in.read(reinterpret_cast<char*>(&val), sizeof(val));
                            data[idx] = val;
                        }
#ifdef DEBUG_TRACE
                        std::cout << "Chunk:\t\t\t\tLoad\nStatus:\t\t\t\tSuccess" << "\nPath:\t\t\t\t" << path <<
                            "\nCount:\t\t\t\t" << data.size() << std::endl << std::endl;
#endif
                        return true;
                    }
#ifdef DEBUG_TRACE
                    std::cout << "Chunk:\t\t\t\tLoad\nStatus:\t\t\t\tFailed" << "\nPath:\t\t\t\t" << path << std::endl << std::endl;
#endif
                    return false;
                }

                bool unload() {
                    if(!is_loaded) return true;
                    is_loaded = false;
                    if constexpr(IS_SAMPLER){
                        std::fstream in (path, in.binary | in.out);
                        if(in.good()){
                            size_t s = data.size();
                            in.write(reinterpret_cast<char*>(&s), sizeof(s));
                            for(auto it = data.begin(); it != data.end(); ++it){
                                uint_fast64_t k = it.key();
                                std::complex<T> v = it.value();
                                in.write(reinterpret_cast<char*>(&k), sizeof(k));
                                in.write(reinterpret_cast<char*>(&v), sizeof(v));
                            }
#ifdef DEBUG_TRACE
                            std::cout << "Chunk:\t\t\t\tUnload\nStatus:\t\t\t\tSuccess" << "\nPath:\t\t\t\t" << path <<
                                "\nCount:\t\t\t\t" << data.size() << std::endl << std::endl;
#endif                        
                            data = tsl::robin_map<uint_fast64_t, std::complex<T>>();
                            return true;
                        }
                    } else {
                        data = tsl::robin_map<uint_fast64_t, std::complex<T>>();
                        return true;
                    }
                    
#ifdef DEBUG_TRACE
                    std::cout << "Chunk:\t\t\t\tUnload\nStatus:\t\t\t\tFailed" << "\nPath:\t\t\t\t\t" << path << std::endl << std::endl;
#endif
                    return false;
                }

                uint_fast64_t idx(const Eigen::VectorX<T>& _pos) const {
                    const Eigen::Vector3<T> pos = Eigen::Vector3<T>(_pos(cardinal(0)), _pos(cardinal(1)), _pos(cardinal(2))) - cell_half;
                    const T x = std::floor(pos(0) / cell_size(0));
                    const T y = std::floor(pos(1) / cell_size(1));
                    const T z = std::floor(pos(2) / cell_size(2));
                    return ::TurboDorn::Hagedorn::Detail::morton_index(Eigen::Vector3<T>(x, y, z));
                }

                bool insert(const Eigen::VectorX<T>& _pos, const std::complex<T>& _val) {
                    if(!is_loaded) return false;
                    const auto id = idx(_pos);
                    data[id] =_val;
#ifdef DEBUG_TRACE
                    std::cout << "Insert Chunk\t\tPosition: [" << _pos(0) << "," << _pos(1) << "," << _pos(2) << "]\t\t" <<
                     "Val: " << _val << "\t\tIdx: " << id << std::endl;
#endif
                    return true;
                }

                std::complex<T> sample(const Eigen::VectorX<T>& _pos) const {
                    const uint_fast64_t id = idx(_pos);
                    const auto it = data.find(id);
                    if(it == data.end()) return {0., 0.};
                    else return it.value();
                }

                bool isLoaded() const noexcept {
                    return is_loaded;
                }

                size_t get_entry_count() const noexcept {
                    return data.size();
                }

                size_t get_byte_size() const noexcept {
                    return data.size() * (sizeof(uint_fast64_t) + sizeof(std::complex<T>));
                }

            };//Chunk

            //---------------------------------------------------------------------------------------//
            //                                        ChunkGrid
            //---------------------------------------------------------------------------------------//

            /*
                A 3d hashed infinite grid that automatically loads and unloads
                chunks to keep the memory footprint under a certain value.
                Uses a simple LRU policy for chunk loading and morton keys for indices.
            */

            template<class T, bool IS_SAMPLER>
            class ChunkGrid {

                Eigen::Vector3i cardinal;
                Eigen::Vector3<T> cell_size;
                Eigen::Vector3<T> chunk_size;
                Eigen::Vector3<T> chunk_half;

                size_t chunk_extend;
                size_t chunk_byte_size;
                size_t max_byte_size_dense;

                tsl::robin_map<uint_fast64_t, std::unique_ptr<Detail::Chunk<T, IS_SAMPLER>>> chunks;

                size_t max_chunks_loaded;
                size_t max_chunk_size;

                std::vector<uint_fast64_t> lru;

                std::filesystem::path path;

                std::mutex mutex;

                bool ensure_loaded(uint_fast64_t _idx) {

                    //check if chunk exists
                   
                    if(!chunks.contains(_idx)) {
                        //if constexpr(IS_SAMPLER) {
                            if(!std::filesystem::exists(path)){
                                std::filesystem::create_directories(path / "chunks");
                            }
                            std::stringstream ss;
                            ss << "chunks/chunk_" << _idx;
                            const auto np = path / ss.str();

                            std::unique_ptr<Detail::Chunk<T, IS_SAMPLER>> nc = std::make_unique<Detail::Chunk<T, IS_SAMPLER>>(std::move(np), cardinal, cell_size);
                            //nc->unload();
                            chunks[_idx] = std::move(nc);
                        //} else return false;
                    }

                    //check if already loaded
                    auto it = lru.begin();
                    for(; it != lru.end(); ++it){
                        if(*it == _idx) break;
                    }

                    //if already loaded just push to the front
                    if(it != lru.end()){
                        lru.erase(it);
                        lru.insert(lru.begin(), _idx);
                        return true;
                    }

                    //if new but lru has still space
                    if(lru.size() < max_chunks_loaded){
                        lru.insert(lru.begin(), _idx);
                        return chunks[_idx]->load();
                    }

                    //if new but lru is full
                    const uint_fast64_t toUnload = *(lru.end() - 1);
                    lru.erase(lru.end() - 1);
                    lru.insert(lru.begin(), _idx);
                    return chunks[toUnload]->unload() && chunks[_idx]->load();

                }

                uint_fast64_t idx(const Eigen::VectorX<T>& _pos) const {
                    const Eigen::Vector3<T> pos = Eigen::Vector3<T>(_pos(cardinal(0)), _pos(cardinal(1)), _pos(cardinal(2))) - chunk_half;
                    const T x = std::floor(pos(0) / chunk_size(0));
                    const T y = std::floor(pos(1) / chunk_size(1));
                    const T z = std::floor(pos(2) / chunk_size(2));
                    return ::TurboDorn::Hagedorn::Detail::morton_index(Eigen::Vector3<T>(x, y, z));
                }

            public:

                /*
                    This is the constructor used for sampling.
                    _path - the path to the folder to store chunk data
                    _cardinal - the 3 main dimensions
                    _cell_size - the size of a single cell in world units
                    _max_chunk_size - the max size of a single chunk in byte
                    max_gridsize_byte - the maximum size the grid in byte. 
                */
                ChunkGrid(
                    std::filesystem::path _path, 
                    const Eigen::Vector3i& _cardinal, 
                    const Eigen::Vector3<T>& _cell_size, 
                    size_t _max_chunk_size = 50000000,
                    size_t max_gridsize_byte = 50000000 * 10
                ) : cardinal(_cardinal), cell_size(_cell_size), chunk_half(cell_size * T(0.5)), max_chunk_size(_max_chunk_size), path(std::move(_path)) {
                    static_assert(IS_SAMPLER == true, "Don't use sampler constructor for rendering.");
                    if(std::filesystem::exists(path)){
                        std::cerr << "[Samplergrid] Chunk data already exists at this path: " << path << std::endl;
                        std::terminate();
                    }
                    chunk_extend = size_t(std::floor(std::sqrt(max_chunk_size / sizeof(std::complex<T>))));
                    chunk_size = cell_size * chunk_extend;
                    max_chunks_loaded = max_gridsize_byte / max_chunk_size;
#ifdef DEBUG_TRACE
                    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
                    std::cout << "SamplerGrid:\t\t\tNew \nPath:\t\t\t\t" << path << "\nCardinal:\t\t\t[" << 
                    _cardinal(0) << ", " << _cardinal(1) << ", " << _cardinal(2) << 
                    "]\nCell Size:\t\t\t[" << 
                    cell_size(0) << ", " << cell_size(1) << ", " << cell_size(2) << 
                    "]\nMax. Chunk Size:\t\t" << max_chunk_size << " bytes \nChunk Extent:\t\t\t" << 
                    chunk_extend << "\nChunk Size:\t\t\t[" << 
                    chunk_size(0) << ", " << chunk_size(1) << ", " << chunk_size(2) << 
                    "]\nMax Chunks Loaded:\t\t" << max_chunks_loaded << std::endl;
                    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;              
#endif
                }

                /*
                    This is the constructor used for rendering.
                    _path - the path to the folder to store chunk data
                    max_gridsize_byte - the maximum size the grid in byte. 
                */
                ChunkGrid(std::filesystem::path _path, size_t max_gridsize_byte = 50000000 * 10) : path(std::move(_path)){
                    static_assert(IS_SAMPLER == false, "Don't use render constructor for sampling.");
                    if(!std::filesystem::exists(path)){
                        std::cerr << "[Rendergrid] Chunk data does not exists at this path: " << path << std::endl;
                        std::terminate();
                    }
                    //load meta data
                    std::fstream in (path / "grid", in.binary | in.in);
                    if(in.good()){
                        in.read(reinterpret_cast<char*>(cardinal.data()), sizeof(cardinal));
                        in.read(reinterpret_cast<char*>(cell_size.data()), sizeof(cell_size));
                        in.read(reinterpret_cast<char*>(&chunk_extend), sizeof(chunk_extend));
                        in.read(reinterpret_cast<char*>(chunk_size.data()), sizeof(chunk_size));
                        in.read(reinterpret_cast<char*>(&max_chunk_size), sizeof(max_chunk_size));
                        max_chunks_loaded = max_gridsize_byte / max_chunk_size;
#ifdef DEBUG_TRACE
                        std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
                        std::cout << "RenderGrid:\t\t\tNew \nPath:\t\t\t\t" << path << "\nCardinal:\t\t\t[" << 
                        cardinal(0) << ", " << cardinal(1) << ", " << cardinal(2) << 
                        "]\nCell Size:\t\t\t[" << 
                        cell_size(0) << ", " << cell_size(1) << ", " << cell_size(2) << 
                        "]\nMax. Chunk Size:\t\t" << max_chunk_size << " bytes \nChunk Extent:\t\t\t" << 
                        chunk_extend << "\nChunk Size:\t\t\t[" << 
                        chunk_size(0) << ", " << chunk_size(1) << ", " << chunk_size(2) << 
                        "]\nMax Chunks Loaded:\t\t" << max_chunks_loaded << std::endl;
                        std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;              
#endif
                    } else {
                        std::cerr << "[RenderGrid] No data found at provided path: " << _path << std::endl;
                        std::terminate();
                    }

                }

                ~ChunkGrid() {
                    //write meta data
                    if constexpr(IS_SAMPLER){
                        std::fstream in (path / "grid", in.binary | in.out);
                        if(in.good()){
                            in.write(reinterpret_cast<char*>(cardinal.data()), sizeof(cardinal));
                            in.write(reinterpret_cast<char*>(cell_size.data()), sizeof(cell_size));
                            in.write(reinterpret_cast<char*>(&chunk_extend), sizeof(chunk_extend));
                            in.write(reinterpret_cast<char*>(chunk_size.data()), sizeof(chunk_size));
                            in.write(reinterpret_cast<char*>(&max_chunk_size), sizeof(max_chunk_size));
                        }
                   
                        size_t byte_out = 0;
                        size_t chunks_out = chunks.size();
                        size_t entries_out = 0;
                        for(auto it = chunks.begin(); it != chunks.end(); ++it){
                            byte_out += it.value()->get_byte_size();
                            entries_out += it.value()->get_entry_count();
                            it.value()->unload();      
                        }
                    std::cout << "Data written:\n" << "Chunks:\t\t\t" << chunks_out << "\nTotal Byte:\t\t" << byte_out << "\nEntries:\t\t" << entries_out << std::endl;
                    }
                }

                //snychronized
                std::complex<T> operator[](const Eigen::VectorX<T>& _pos) {
                    std::lock_guard<std::mutex> lock (mutex);
                    if(ensure_loaded(idx(_pos))){
                        const std::complex<T> val = chunks[idx(_pos)]->sample(_pos);
#ifdef DEBUG_TRACE
                        std::cout << "Sample: \tPosition: [" << _pos(0) << "," << _pos(1) << "," << _pos(2) << "]\t" <<
                            "Val: " << val << std::endl;
#endif                        
                        return val;
                    }
                    return {0., 0.};                
                }

                //snychronized
                bool insert(const Eigen::VectorX<T>& _pos, const std::complex<T>& _val){
                    if(std::abs(_val.real()) < std::numeric_limits<T>::epsilon() && std::abs(_val.imag()) < std::numeric_limits<T>::epsilon()) return true;
                    std::lock_guard<std::mutex> lock (mutex);
                    const auto id = idx(_pos);
                    if(ensure_loaded(id)){
                        const bool r = chunks[id]->insert(_pos, _val);
#ifdef DEBUG_TRACE
                        if(r)
                            std::cout << "Insert Grid: Success\tPosition: [" << _pos(0) << "," << _pos(1) << "," << _pos(2) << "]\t" <<
                            "Val: " << _val << "\tIdx: " << id << std::endl;
                        else
                            std::cout << "Insert Grid: Fail\tPosition: [" << _pos(0) << "," << _pos(1) << "," << _pos(2) << "]\t" <<
                            "Val: " << _val << "\tIdx: " << id << std::endl;
#endif
                        return r;
                    }
#ifdef DEBUG_TRACE
                    std::cout << "Insert: Ensure load failed\t\tPosition: [" << _pos(0) << "," << _pos(1) << "," << _pos(2) << "]\t" <<
                        "Val: " << _val << std::endl;
#endif
                    return false;
                }

            };

        }//Detail

        template<class T>
        using SamplerGrid = Detail::ChunkGrid<T, true>;

        template<class T>
        using RenderGrid = Detail::ChunkGrid<T, false>;

    }
    
    namespace Hagedorn::Detail {
        std::vector<Eigen::VectorXi> hyperbolic_cut_shape(size_t _dim, size_t _K);
        Eigen::Index index(const Eigen::VectorXi& _i, const Eigen::VectorXi& _e ) noexcept;
    }

    namespace IO {

        namespace Detail {

            //---------------------------------------------------------------------------------------//
            //                                        FilePathResolver
            //---------------------------------------------------------------------------------------//
            /*
                A functional that returns the path to the example files folder
            */
            class FilePathResolver {

                std::filesystem::path path_asset;
            
            public:

                FilePathResolver() {
                    //asset
                    {
                        const auto s_path_1 = std::filesystem::current_path() /= "example_files/";
                        const auto s_path_2 = std::filesystem::current_path() /= "../example_files/";
                        const auto s_path_3 = std::filesystem::current_path() /= "../../example_files/";

                        if(std::filesystem::exists(s_path_1))
                            path_asset = s_path_1;
                        else if(std::filesystem::exists(s_path_2))
                            path_asset = s_path_2;
                        else if(std::filesystem::exists(s_path_3))
                            path_asset = s_path_3;
                        else throw std::runtime_error("folder example_files not found");
                    }
                }

                std::filesystem::path operator()() const {
                    return path_asset;
                }

            };//FilePathResolver

            //---------------------------------------------------------------------------------------//
            //                                        File
            //---------------------------------------------------------------------------------------//

            template<class T>
            struct File {
                size_t dimensions;
                size_t timesteps;
                size_t K;
                T epsilon = 1.;  
                Eigen::VectorXcd S;       
                std::vector<Eigen::VectorX<std::complex<T>>> c_0;
                std::vector<Eigen::VectorX<std::complex<T>>> p, q;
                std::vector<Eigen::MatrixX<std::complex<T>>> P, Q;

                //shapefunction
                std::vector<Eigen::VectorXi> Ks; 
                Eigen::VectorXi k_max; //max k in d direction
                tsl::robin_map<Eigen::Index, bool> b_Ks;
            };//File

            //---------------------------------------------------------------------------------------//
            //                                   load_from_file
            //---------------------------------------------------------------------------------------//

            template<class T, class I, class CUTSHAPE>
            File<T> load_from_file(const std::filesystem::path& _path, I _dims, I _K, CUTSHAPE&& _cs) {

                using Vector = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>;
                using Matrix = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;

                if(!std::filesystem::exists(_path))
                    throw std::runtime_error("file does not exist");

                File<T> out;
                out.dimensions = _dims;
                out.K = _K;

                try {

                    //file
                    HighFive::File file(_path.string(), HighFive::File::ReadOnly);

                    //timesteps
                    {
                        const Eigen::Matrix<T, Eigen::Dynamic, 1> steps = H5Easy::load<Eigen::Matrix<T, Eigen::Dynamic, 1>>(file, "datablock_0/wavepacket/timegrid");
                        out.timesteps = steps.rows();
                    }

                    //S
                    {
                        const H5Easy::DataSet dataset = file.getDataSet("datablock_0/wavepacket/Pi/S");
                        std::vector<std::vector<std::vector<std::complex<T>>>> S;
                        dataset.read(S);
                        out.S.resize(S.size());
                        for(size_t i = 0; i < S.size(); ++i){
                            out.S(i) = S[i][0][0];
                        }
                    }
                
                    //c_0
                    {
                        const Matrix tmp = H5Easy::load<Matrix>(file, "datablock_0/wavepacket/coefficients/c_0");
                        for(size_t t = 0; t < out.timesteps; ++t){
                            Vector c_0(tmp.cols());
                            for(size_t d = 0; d < tmp.cols(); ++d)
                                c_0(d) = tmp(t, d);
                            out.c_0.push_back(std::move(c_0));
                        }
                    }
                
                    //p
                    {
                        const H5Easy::DataSet dataset = file.getDataSet("datablock_0/wavepacket/Pi/p");
                        std::vector<std::vector<std::vector<std::complex<T>>>> p;
                        dataset.read(p);
                        for(size_t t = 0; t < p.size(); ++t){
                            Vector temp(out.dimensions);
                            for(size_t y = 0; y < p[t].size(); ++y){
                                temp(y) = p[t][y][0];
                            }
                            out.p.push_back(std::move(temp));
                        }
                    }
                
                    //q
                    {
                        const H5Easy::DataSet dataset = file.getDataSet("datablock_0/wavepacket/Pi/q");
                        std::vector<std::vector<std::vector<std::complex<T>>>> q;
                        dataset.read(q);
                        for(size_t t = 0; t < q.size(); ++t){
                            Vector temp(out.dimensions);
                            for(size_t y = 0; y < q[t].size(); ++y){
                                temp(y) = q[t][y][0];
                            }
                            out.q.push_back(std::move(temp));
                        }
                    }
            
                    //P
                    {
                        const H5Easy::DataSet dataset = file.getDataSet("datablock_0/wavepacket/Pi/P");
                        std::vector<std::vector<std::vector<std::complex<T>>>> P;
                        dataset.read(P);
                        for(size_t t = 0; t < P.size(); ++t){
                            Matrix temp(out.dimensions, out.dimensions);
                            for(size_t x = 0; x < out.dimensions; ++x){
                                for(size_t y = 0; y < out.dimensions; ++y)
                                    temp(x, y) = P[t][x][y];
                            }
                            out.P.push_back(std::move(temp));
                        }
                    }

                    //Q
                    {
                        const H5Easy::DataSet dataset = file.getDataSet("datablock_0/wavepacket/Pi/Q");
                        std::vector<std::vector<std::vector<std::complex<T>>>> Q;
                        dataset.read(Q);
                        for(size_t t = 0; t < Q.size(); ++t){
                            Matrix temp(out.dimensions, out.dimensions);
                            for(size_t x = 0; x < out.dimensions; ++x){
                                for(size_t y = 0; y < out.dimensions; ++y)
                                    temp(x, y) = Q[t][x][y];
                            }
                            //std::cout << temp << std::endl;
                            out.Q.push_back(std::move(temp));
                        }
                    }
                    

                } catch(const std::exception& _e){
                    throw _e;
                }

                // -------------- cut shape --------------
                {
                    out.Ks = _cs(out.dimensions, out.K);

                    out.k_max = Eigen::VectorXi(out.dimensions);
                    out.k_max.setZero();

                    for(size_t i = 0; i < out.dimensions; ++i){
                        for(const auto& k : out.Ks){
                            out.k_max(i) = std::max(out.k_max(i), k(i));
                        }
                    }
            
                    for(const auto& c : out.Ks){
                        const Eigen::Index ii = ::TurboDorn::Hagedorn::Detail::index(c, out.k_max);
                        //std::cout << c(0) << ", " << c(1) << ", " << c(2) << " | " << out.k_max(0) << ", " << out.k_max(1) << ", " << out.k_max(2) << " | " << ii << std::endl;
                        out.b_Ks.insert( {ii, true} );
                    }
                }

                return out ;
            }; //loadFromFile

        }

        Detail::File<double> simulation_results() {
            try {
                const auto path = Detail::FilePathResolver()() /= "simulation_results_phi.hdf5";
                return Detail::load_from_file<double> (path, 3, 4, ::TurboDorn::Hagedorn::Detail::hyperbolic_cut_shape);
            } catch(const std::exception& _e){
                throw _e;
            }
        }//simulation_results

        Detail::File<double> simulation_results_phi000() {
            try {
                const auto path = Detail::FilePathResolver()() /= "simulation_results_phi000.hdf5";
                return Detail::load_from_file<double>(path, 3, 1, ::TurboDorn::Hagedorn::Detail::hyperbolic_cut_shape);
            } catch(const std::exception& _e){
                throw _e;
            }
        }//simulation_results_phi000

        Detail::File<double> simulation_results_phi100() {
            try {
                const auto path = Detail::FilePathResolver()() /= "simulation_results_phi100.hdf5";
                return Detail::load_from_file<double>(path, 3, 2, ::TurboDorn::Hagedorn::Detail::hyperbolic_cut_shape);
            } catch(const std::exception& _e){
                throw _e;
            }
        }//simulation_results_phi100

        Detail::File<double> simulation_results_phi121() {
            try {
                const auto path = Detail::FilePathResolver()() /= "simulation_results_phi121.hdf5";
                return Detail::load_from_file<double>(path, 3, 12, ::TurboDorn::Hagedorn::Detail::hyperbolic_cut_shape);
            } catch(const std::exception& _e){
                throw _e;
            }
        }//simulation_results_phi121

        Detail::File<double> simulation_results_phi412() {
            try {
                const auto path = Detail::FilePathResolver()() /= "simulation_results_phi412.hdf5";
                return Detail::load_from_file<double>(path, 3, 30, ::TurboDorn::Hagedorn::Detail::hyperbolic_cut_shape);
            } catch(const std::exception& _e){
                throw _e;
            }
        }//simulation_results_phi412

    }//IO

    namespace Hagedorn {

        namespace Detail {

            //---------------------------------------------------------------------------------------//
            //                                        Invariants
            //---------------------------------------------------------------------------------------//

            template<class T>
            struct Invariants {
                size_t dimensions;
                Eigen::VectorXi k; //k extends
                tsl::robin_map<Eigen::Index, bool> k_shape; //lookup for shape
                std::vector<Eigen::VectorX<std::complex<T>>> p;
                std::vector<Eigen::VectorX<std::complex<T>>> q;
                //phi_0
                std::vector<std::complex<T>> pre;
                std::complex<T> i_2_E_2;
                std::vector<Eigen::MatrixX<std::complex<T>>> P_Q_1;
                std::vector<Eigen::RowVectorX<std::complex<T>>> i_E_2_p;
                //std::vector<Eigen::Matrix<std::complex<T>, 1, Eigen::Dynamic>> i_E_2_p;
                //phi
                //sqrt*Q-1
                std::vector<Eigen::MatrixX<std::complex<T>>> Q_1;
                //Q-1*QT
                std::vector<Eigen::MatrixX<std::complex<T>>> Q_1_Q_T;
                //
                std::vector<Eigen::VectorX<std::complex<T>>> c_0;
                std::vector<Eigen::VectorXi> Ks; 
                Eigen::VectorXi k_max; //max k in d direction
                Eigen::VectorXcd S;  
            };//Invariants

            //---------------------------------------------------------------------------------------//
            //                             prepare_invariants_from_file
            //---------------------------------------------------------------------------------------//

            template<class T, template<typename> class FILE>
            EIGEN_STRONG_INLINE std::unique_ptr<Invariants<T>> prepare_invariants_from_file(const FILE<T>& _file) {
                std::unique_ptr<Invariants<T>> out = std::make_unique<Invariants<T>>();
                out->dimensions = _file.dimensions;
                out->k = _file.k_max;
                out->i_2_E_2 = std::complex<T>(0., 1.) / (2. * _file.epsilon * _file.epsilon);
                out->p = _file.p;
                out->q = _file.q;
                out->k_shape = _file.b_Ks;

                for(size_t t = 0; t < _file.timesteps; ++t){
                    //phi0
                    out->pre.push_back( std::pow(T(M_PI) * _file.epsilon * _file.epsilon, - T(_file.dimensions) / 4.) * std::pow(_file.Q[t].determinant(), -0.5) );
                    out->P_Q_1.push_back( _file.P[t] * _file.Q[t].inverse() );

                    out->i_E_2_p.push_back( (std::complex<T>(0., 1.) / _file.epsilon * _file.epsilon) * _file.p[t].transpose() );
                    //phi
                    out->Q_1.push_back( std::sqrt(2. / (_file.epsilon * _file.epsilon)) * _file.Q[t].inverse() );
                    out->Q_1_Q_T.push_back( _file.Q[t].inverse() * _file.Q[t].conjugate() );
                }
                
                // -------------- copy some things --------------
                {
                    out->c_0 = _file.c_0;
                    out->Ks = _file.Ks;
                    out->k_max = _file.k_max;
                    out->S = _file.S;
                }
                return out;
            }

            //---------------------------------------------------------------------------------------//
            //                                            morton_index
            //---------------------------------------------------------------------------------------//

            /*
             * _x: indices to be encoded in signed form
            */
            template<class T>
            EIGEN_STRONG_INLINE uint_fast64_t morton_index(const Eigen::Vector3<T>& _x) {
                const uint_fast64_t out = libmorton::morton3D_64_encode(uint_fast64_t(std::abs(_x(0))), uint_fast64_t(std::abs(_x(1))), uint_fast64_t(std::abs(_x(2))));
                return (out & 0x1FFFFFFFFFFFFFFF) | uint_fast64_t(std::signbit(_x(0))) << 63 | uint_fast64_t(std::signbit(_x(1))) << 62 | uint_fast64_t(std::signbit(_x(2))) << 61;
            }

            //---------------------------------------------------------------------------------------//
            //                                            index
            //---------------------------------------------------------------------------------------//

            /*
            i: index
            e: extents, # units
            */
            EIGEN_STRONG_INLINE Eigen::Index index(const Eigen::VectorXi& _i, const Eigen::VectorXi& _e ) noexcept {
                assert(_i.size() == _e.size());
                Eigen::Index out = _i(0);
                for (Eigen::Index k = 1; k < _i.size(); ++k) {
                    out *= _e(k);
                    out += _i(k);
                }
                return out;
            }; //index

            //---------------------------------------------------------------------------------------//
            //                             hyperbolic_cut_shape
            //---------------------------------------------------------------------------------------//

            EIGEN_STRONG_INLINE std::vector<Eigen::VectorXi> hyperbolic_cut_shape(size_t _dim, size_t _K) {
                std::vector<Eigen::VectorXi> out;

                Eigen::VectorXi index(_dim);
                index.setZero();

                while (true) {
                    for (index(_dim - 1) = 0; index(_dim - 1) <= _K; ++index(_dim - 1)) {
                        Eigen::Index p = 1;
                        for (Eigen::Index d = 0; d < _dim; ++d)
                            p *= (1 + index(d));
                        if (p <= _K)
                            out.push_back(index);
                    }

                    bool done = false;
                    for (Eigen::Index d = _dim - 2; d >= 0; --d) {
                        index(d) += 1;

                        if (index(d) >= _K) {
                            if (d == 0)
                                done = true;
                            else
                                index(d) = 0;
                        } else
                            break;
                    }
                    if (done) break;
                }

                return out;
            }; //hyperbolic_cut_shape

            //---------------------------------------------------------------------------------------//
            //                                    ray_aabb_intersect
            //---------------------------------------------------------------------------------------//

            template<class T>
            EIGEN_STRONG_INLINE bool ray_aabb_intersect(const Eigen::Vector3<T>& _r_o, const Eigen::Vector3<T>& _r_d, const Eigen::Vector3<T>& _low, const Eigen::Vector3<T>& _high, T _tmax, T& _t) noexcept {
                _t = -std::numeric_limits<T>::infinity();
                for (size_t i = 0; i < _r_o.size(); ++i) {
                    if (std::abs(_r_d[i]) < std::numeric_limits<T>::epsilon()){
                        if (_r_o[i] < _low[i] || _r_o[i] > _high[i]) 
                            return false;
                    } else {
                        const T ood = 1. / _r_d[i];
                        T t1 = (_low[i] - _r_o[i]) * ood;
                        T t2 = (_high[i] - _r_o[i]) * ood;
                        if (t1 > t2) std::swap(t1, t2);
                        _t = std::max(_t, t1);
                        _tmax = std::min(_tmax, t2);
                        if (_t > _tmax) return false;
                    }
                }
                return true;
            }//intersect

            //---------------------------------------------------------------------------------------//
            //                                        c_to_HSL
            //---------------------------------------------------------------------------------------//

            template <class T>
            EIGEN_STRONG_INLINE constexpr Eigen::Vector3<T> c_to_HSL(const std::complex<T>& _c, T MAXVALUE = 10.) noexcept {
                const T H = std::clamp(std::abs(std::fmod(std::arg(_c), 2. * M_PI)), 0., 2. * M_PI);
                const T S = 1.;
                const T L = std::clamp(std::abs(MAXVALUE * std::atan(std::abs(_c)) / (0.5 * M_PI)), 0., 1.);
                return { H, S, L };
            }//c_to_HSL

            //---------------------------------------------------------------------------------------//
            //                                        HSL_to_RGB_deg
            //---------------------------------------------------------------------------------------//

            /*.
                h: [0, 360]
                s: [0, 1]
                l: [0, 1]
                rgb: [0, 1]
            */
            template <class T>
            EIGEN_STRONG_INLINE constexpr Eigen::Vector3<T> HSL_to_RGB_deg(const Eigen::Matrix<T, 3, 1>& _hsl) noexcept {

                const T H = _hsl(0);
                const T S = _hsl(1);
                const T L = _hsl(2);

                assert(0. <= H && H <= 360.);
                assert(0. <= S && S <= 1.);
                assert(0. <= L && L <= 1.);

                const T C = ( T(1.) - std::abs( T(2.) * L - T(1.) ) ) * S;
                const T X = C * (T(1.) - std::abs(std::fmod(H / T(60.), T(2.)) - T(1.)));
                const T m = L - C * T(0.5);

                switch(size_t(H / 60.)){
                    case 0: return { C + m, X + m, m};
                    case 1: return { X + m, C + m, m};
                    case 2: return { m, C + m, X + m};
                    case 3: return { m, X + m, C + m};
                    case 4: return { X + m, m, C + m};
                    case 5: return { C + m, m, X + m};
                    default: return { 0., 0., 0.};
                }

            }; //HSL_to_RGB_deg

            //---------------------------------------------------------------------------------------//
            //                                        HSL_to_RGB_rad
            //---------------------------------------------------------------------------------------//

            /*
                h: [0, 2pi]
                s: [0, 1]
                l: [0, 1]
                rgb: [0, 1]
            */
            template <class T>
            EIGEN_STRONG_INLINE constexpr Eigen::Vector3<T> HSL_to_RGB_rad(const Eigen::Matrix<T, 3, 1>& _hsl) noexcept {
                return HSL_to_RGB_deg<T>( { _hsl(0) * T( 180. / M_PI), _hsl(1), _hsl(2) } );
            }; //HSL_to_RGB_rad

            //---------------------------------------------------------------------------------------//
            //                                        rgb_to_gs
            //---------------------------------------------------------------------------------------//

            template <class T>
            EIGEN_STRONG_INLINE constexpr Eigen::Vector3<T> rgb_to_gs(const Eigen::Matrix<T, 3, 1>& _rgb) noexcept {
                const T gs = T(0.299) * _rgb(0) + T(0.587) * _rgb(1) +  T(0.114) * _rgb(2);
                return { gs, gs, gs };
            }; //rgb_to_gs

            //---------------------------------------------------------------------------------------//
            //                                     phi_0
            //---------------------------------------------------------------------------------------//

            template <class T, template<typename> class INVARIANTS>
            EIGEN_STRONG_INLINE std::complex<T> phi_0 (
                size_t _t, 
                const Eigen::VectorX<T>& _x, 
                const INVARIANTS<T>& _inv
            ) noexcept {
                const Eigen::VectorX<std::complex<T>> xq = _x - _inv.q[_t];
                const Eigen::RowVectorX<std::complex<T>> xqt = xq.transpose();
                const std::complex<T> e1 = _inv.i_2_E_2 * xqt * _inv.P_Q_1[_t] * xq;
                const std::complex<T> e2 = _inv.i_E_2_p[_t] * xq;
                return _inv.pre[_t] * std::exp(e1 + e2);
            }; //phi_0

            //---------------------------------------------------------------------------------------//
            //                                        phi
            //---------------------------------------------------------------------------------------//

            template <class T, template<typename> class INVARIANTS>
            EIGEN_STRONG_INLINE Eigen::VectorX<std::complex<T>> phi (
                size_t _t,
                const Eigen::VectorX<T>& _x,
                const tsl::robin_map<Eigen::Index, std::complex<T>>& _phis,
                const Eigen::VectorXi& _index,
                const INVARIANTS<T>& _inv
            ) noexcept {

                using Index = Eigen::VectorXi;
                using Vector = Eigen::VectorX<std::complex<T>>;
                using Matrix = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>;

                const auto xq = _x - _inv.q[_t];

                Vector res (_inv.dimensions);

                Vector kp (_inv.dimensions);
                for (size_t j = 0; j < _inv.dimensions; ++j) {
                    if (_index(j) - 1 < 0) {
                        kp(j) = { 0., 0. };
                        continue;
                    }

                    Index k_1 = _index;
                    k_1(j) = std::max<Eigen::Index>(k_1(j) - 1, 0ll);
                    const Eigen::Index ii = Detail::index(k_1, _inv.k);
                    assert(_phis.contains(ii));
                    const std::complex<double>& p = (*_phis.find(ii)).second;
                    kp(j) = std::sqrt(_index(j)) * p;  
                }

                const Eigen::Index ii = Detail::index(_index, _inv.k);
                assert(_phis.contains(ii));
                const std::complex<double>& p = (*_phis.find(ii)).second;
                auto phi_t = _inv.Q_1[_t] * xq * p - _inv.Q_1_Q_T[_t] * kp;

                Vector phi (_inv.dimensions);
                for(size_t i = 0; i < _inv.dimensions; ++i){
                    const T sk = std::sqrt(T(_index(i)) + 1.);
                    const std::complex<T> skc = std::complex<T>(sk, 0.);
                    phi(i) = phi_t(i) / skc;
                }

                return phi;
            }; //phi

            //---------------------------------------------------------------------------------------//
            //                                  compute_cube
            //---------------------------------------------------------------------------------------//

            template <class T, template<typename> class INVARIANTS>
            EIGEN_STRONG_INLINE tsl::robin_map<Eigen::Index, std::complex<T>> compute_cube (
                size_t _t, 
                const Eigen::VectorX<T>& _x, 
                const INVARIANTS<T>& _inv
            ) {

                using Index = Eigen::VectorXi;
                using Vector = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>;
                using Matrix = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>;

                const size_t dim = _inv.dimensions;

                size_t size = _inv.k(0) + 1;
                for (size_t i = 1; i < _inv.k.size(); ++i)
                    size *= (_inv.k(i)+1);

                tsl::robin_map<Eigen::Index, std::complex<T>> phis;

                //iterate over ks
                Index index(dim);
                index.fill(0);

                bool first = true;

                while (true) {
                    for (index(dim-1) = 0; index(dim-1) <= _inv.k(dim-1); ++index(dim-1)) { 
                        if (first) {
                            first = false;
                            const auto phi0 = Detail::phi_0<T>(_t, _x, _inv);
                            phis.insert( {0, phi0} );
                            --index(dim-1);
                            continue;
                        }

                        //check if have reached the end of the shape
                        const Eigen::Index ii = Detail::index(index, _inv.k);
                        if(!_inv.k_shape.contains(ii))
                            break;

                        //compute phi for index
                        const auto phi = Detail::phi<T>(_t, _x, phis, index, _inv);

                        for (size_t d = 0; d < dim; ++d) {
                            Index ni = index;
                            ni(d) += 1;
                            const Eigen::Index ii = Detail::index(ni, _inv.k);
                            phis.insert( {ii, phi(d)} );
                            //phis[ii] = phi(d);
                        }
                    }

                    bool done = false;
                    for (Eigen::Index d = dim - 2; d >= 0; --d) {
                        index(d) += 1;
                        if (index(d) >=  _inv.k(d)) {
                            if (d == 0)
                                done = true;
                            else
                                index(d) = 0;
                        } else
                            break;
                    }
                    if (done) break;
                }

                return phis;
            } //compute

            //---------------------------------------------------------------------------------------//
            //                                   linear_combination
            //---------------------------------------------------------------------------------------//

            template <class T, template<typename> class INVARIANTS>
            std::complex<T> linear_combination(
                size_t _time_step,
                const tsl::robin_map<Eigen::Index, std::complex<T>>& _phis,
                const INVARIANTS<T>& _inv
            ) {
                std::complex<T> res (0., 0.);
                for(Eigen::Index k = 0; k < _inv.Ks.size(); ++k){ 
                    const Eigen::Index idx = Detail::index( _inv.Ks[k],  _inv.k_max);
                    assert(_phis.count(idx) != 0);
                    const std::complex<T>& p = (*_phis.find(idx)).second;
                    res +=  _inv.c_0[_time_step](k) * p;
                } 
                res *= std::exp(std::complex<double>(0., 1.) *  _inv.S(Eigen::Index(_time_step)));
                return res;
            }

        }//Detail

        //---------------------------------------------------------------------------------------//
        //                                        HyperCube
        //---------------------------------------------------------------------------------------//
        /*
            A functional that computes the functional value at a given world position and time step.
        */
        template <class T>
        class HyperCube {

            std::complex<T> value;

        public:

            template<template<typename> class INVARIANTS>
            HyperCube(size_t _time_step, const Eigen::VectorX<T>& _pos, const INVARIANTS<T>& _inv) {
                const auto phis = Detail::compute_cube(_time_step, _pos, _inv);
                value = linear_combination(_time_step, phis, _inv);
            }

            HyperCube() = delete;

            EIGEN_STRONG_INLINE std::complex<T> VAL() const noexcept {
                return value;
            }

            /**
             * @brief Returns the functions value in HSL with max value of 10.
            */
            EIGEN_STRONG_INLINE Eigen::Vector3<T> HSL(T MAXVALUE = 10.) const noexcept {
                return Detail::c_to_HSL<T>(value, MAXVALUE);
            }

            /**
             * @brief Returns the functions value in RGB
            */
            EIGEN_STRONG_INLINE Eigen::Vector3<T> RGB() const noexcept {
                return Detail::HSL_to_RGB_deg(HSL());
            }

            /**
             * @brief Returns the functions value in Grayscale
            */
            EIGEN_STRONG_INLINE Eigen::Vector3<T> GRAYSCALE() const noexcept {
                return Detail::rgb_to_gs(RGB());
            }

        };//HyperCube

    }//Hagedorn

    namespace Detail {
        
        //---------------------------------------------------------------------------------------//
        //                                   insert_sample_into
        //---------------------------------------------------------------------------------------//

        template<class T, template<typename> class INVARIANTS>
        EIGEN_STRONG_INLINE bool insert_sample_into_grid(
            Geometry::SamplerGrid<T>& _grid,
            Eigen::Index _time_step, 
            const Eigen::Vector3i& _idx,
            const Eigen::Vector3i& _cardinal,
            const Eigen::Vector3<T>& _cell_extent,
            const Eigen::VectorX<T>& _aabb_min,
            const INVARIANTS<T>& _inv
        ) {
            Eigen::VectorX<T> pos = _aabb_min; //todo: errör
            for(Eigen::Index i = 0; i < 3; ++i)
                pos(_cardinal(i)) += _idx(i) * _cell_extent(i);
            const Hagedorn::HyperCube<T> cube(_time_step, pos, _inv);
            return _grid.insert(pos, cube.VAL());
        }//insert_sample_into_grid

        //---------------------------------------------------------------------------------------//
        //                                   sample_grid
        //---------------------------------------------------------------------------------------//

        template<class T, template<typename> class INVARIANTS, class CHAR = char>
        void sample_grid(
            const std::filesystem::path& _output,
            Eigen::Index _time_step, 
            const Eigen::Vector3i& _cardinal,
            const Eigen::Vector3i& _lattice,
            const Eigen::Vector3<T>& _aabb_min,
            const Eigen::Vector3<T>& _aabb_max,
            const Eigen::VectorX<T>& _const_dims,
            const INVARIANTS<T>& _inv
        ){
            const Eigen::Vector3<T> ext = _aabb_max - _aabb_min;
            const Eigen::Vector3<T> c = { ext(_cardinal(0)) / T(_lattice(0)), ext(_cardinal(1)) / T(_lattice(1)), ext(_cardinal(2)) / T(_lattice(2)) };

            std::cout << "------------ Sampling start ------------" << std::endl;
            const auto start = std::chrono::high_resolution_clock::now();
            {
                Geometry::SamplerGrid<T> grid(_output, _cardinal, c);

                #pragma omp parallel for schedule(dynamic) collapse(3) shared(grid)
                for(Eigen::Index z = 0; z < _lattice(2); ++z){
                    for(Eigen::Index y = 0; y < _lattice(1); ++y){
                        for(Eigen::Index x = 0; x < _lattice(0); ++x){

                            const bool s = insert_sample_into_grid<T>(grid, _time_step, Eigen::Vector3i(x, y, z), _cardinal, c, _aabb_min, _inv);
                            if(!s){
                                std::cerr << "Error" << std::endl;
                            }

                        }
                    }//1119 1168 369 387
                }

            }
            const auto end = std::chrono::high_resolution_clock::now();
            std::cout << "Total Samples:\t\t" << _lattice(0) * _lattice(1) * _lattice(2) << std::endl;
            std::cout << "------------ Sampling done after " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "s ------------\n";
        }

        //---------------------------------------------------------------------------------------//
        //                                   render_grid
        //---------------------------------------------------------------------------------------//
        template<class T>
        std::vector<T> render_grid(
            const std::filesystem::path& _input,
            const Eigen::Vector3<T>& _cam_or,
            const Eigen::Vector3<T>& _cam_dir,
            const Eigen::Vector3<T>& _cam_right,
            const Eigen::Vector2i& _img_size_pixel,
            T _pixel_world_size,
            const Eigen::Vector3<T>& _aabb_min,
            const Eigen::Vector3<T>& _aabb_max,
            T step_size,
            T _hsl_max_value,
            Geometry::RenderGrid<T>& _grid
        ) {

            const auto contains = [](const Eigen::Vector3<T>& _min, const Eigen::Vector3<T>& _max, const Eigen::Vector3<T>& _pos) -> bool {
                for(Eigen::Index i = 0; i < 3; ++i)
                    if(_pos(i) < _min(i) || _pos(i) > _max(i)) return false;
                return true;
            };

            std::vector<T> buffer (_img_size_pixel(0) * _img_size_pixel(1) * 4);
            std::fill(buffer.begin(), buffer.end(), 0.);

            const Eigen::Vector3<T> cam_up = (_cam_dir.cross(_cam_right)).normalized();
            const T dS = step_size;

            const T max_dist = ((_aabb_min + (_aabb_max - _aabb_min) * 0.5) - _cam_or).norm() * 10.;

            const T dp = ((2 * M_PI) / T(36));

            //#pragma omp parallel for schedule(dynamic) shared(buffer) num_threads(4)
            for(size_t x = 0; x < _img_size_pixel(0); ++x){
                for(size_t y = 0; y < _img_size_pixel(1); ++y){

                    const Eigen::Vector2<T> pp = Eigen::Vector2<T>((2. * T(x) - T(_img_size_pixel(0))) / T(_img_size_pixel.maxCoeff()), (2. * T(y) - T(_img_size_pixel(1))) / T(_img_size_pixel.maxCoeff()));
                    const Eigen::Vector3<T> rr = _cam_or + _pixel_world_size * (pp(0) * _cam_right + pp(1) * cam_up);

                    T tt = 0.;
                    const bool hit = Hagedorn::Detail::ray_aabb_intersect(rr, _cam_dir, _aabb_min, _aabb_max, max_dist, tt);
  
                    if(!hit){
                        buffer[4 * (y * _img_size_pixel(0) + x) + 3] = 1.;
                        continue;
                    }

                    T transmission = 1.;
                    Eigen::Vector3<T> col = Eigen::Vector3<T>::Zero();
                    
                    while(true){

                        const Eigen::VectorX<T> pos = rr + _cam_dir * tt;
                        if(!contains(_aabb_min, _aabb_max, pos)) break;

                        const std::complex<T> res = _grid[pos];

                        //compute color
                        const auto hsl = Hagedorn::Detail::c_to_HSL(res, _hsl_max_value);
                        transmission *= std::exp(-hsl(2) * dS);
                        const auto rgb = Hagedorn::Detail::HSL_to_RGB_rad(hsl);
                        col += transmission * rgb * dS;
                        tt += dS;

                    }

                    buffer[4 * (y * _img_size_pixel(0) + x)] = col(0);
                    buffer[4 * (y * _img_size_pixel(0) + x) + 1] = col(1);
                    buffer[4 * (y * _img_size_pixel(0) + x) + 2] = col(2);
                    buffer[4 * (y * _img_size_pixel(0) + x) + 3] = 1.;

                    //std::cout << "[" << col(0) << ", " << col(1) << ", " << col(2) << "]" << std::endl;

                }
            }

            return buffer;
        }

        //---------------------------------------------------------------------------------------//
        //                                   normalize_to_char
        //---------------------------------------------------------------------------------------//
        template<class T>
        EIGEN_STRONG_INLINE std::vector<unsigned char> normalize_to_char(const std::vector<T>& _in) {
            using uchar = unsigned char;
            std::vector<unsigned char> out (_in.size());
            T max = -std::numeric_limits<T>::infinity();
            for(size_t i = 0; i < out.size(); ++i){
                max = std::max<T>(max, _in[i]);
            }
            for(size_t i = 0; i < out.size()-4; i+=4){
                out[i] = uchar((_in[i] / max) * 255);
                out[i+1] = uchar((_in[i+1] / max) * 255);
                out[i+2] = uchar((_in[i+2] / max) * 255);
                out[i+3] = 255;
            }
            return out;
        }

    }

    //---------------------------------------------------------------------------------------//
    //                                        Sampler
    //---------------------------------------------------------------------------------------//
    /*
        An example implementation that shows how to sample a give file from the config json.
    */
    template<class T>
    class Sampler {
    public:
        Sampler(const json& _config) {

            const auto extract = [](const json& _config) -> Eigen::VectorX<T> {
                Eigen::VectorX<T> out;
                out.resize(_config.size());
                Eigen::Index i = 0;
                for(auto it = _config.begin(); it != _config.end(); ++it){
                    out.coeffRef(i++) = T(*it);
                }
                return out;
            };

            const auto extract3 = [](const json& _config) -> Eigen::Vector3<T> {
                Eigen::Vector3<T> out;
                Eigen::Index i = 0;
                for(auto it = _config.begin(); it != _config.end(); ++it){
                    out.coeffRef(i++) = T(*it);
                }
                return out;
            };

            const auto extract3i = [](const json& _config) -> Eigen::Vector3i {
                Eigen::Vector3i out;
                Eigen::Index i = 0;
                for(auto it = _config.begin(); it != _config.end(); ++it){
                    out.coeffRef(i++) = int(*it);
                }
                return out;
            };

            //input file (must exist)
            if(!_config.contains("input")) {
                std::cerr << "[Json Error]: field 'input' missing. Aborting..." << std::endl;
                std::terminate();
            }

            const std::filesystem::path input = std::filesystem::path(_config["input"]);
            if(!std::filesystem::exists(input)){
                std::cerr << "[IO Error]: Input file does not exist. Aborting..." << std::endl;
                std::terminate();
            }

            if(!_config.contains("time step")) {
                std::cerr << "[Json Error]: field 'time step' missing. Aborting..." << std::endl;
                std::terminate();
            }

            if(!_config.contains("dims")) {
                std::cerr << "[Json Error]: field 'dims' missing. Aborting..." << std::endl;
                std::terminate();
            }

            if(!_config.contains("K")) {
                std::cerr << "[Json Error]: field 'K' missing. Aborting..." << std::endl;
                std::terminate();
            }

            if(!_config.contains("cardinal")) {
                std::cerr << "[Json Error]: field 'cardinal' missing. Aborting..." << std::endl;
                std::terminate();
            }

            if(!_config.contains("aabb min")) {
                std::cerr << "[Json Error]: field 'aabb min' missing. Aborting..." << std::endl;
                std::terminate();
            }

            if(!_config.contains("aabb max")) {
                std::cerr << "[Json Error]: field 'aabb max' missing. Aborting..." << std::endl;
                std::terminate();
            }

            if(!_config.contains("lattice")) {
                std::cerr << "[Json Error]: field 'lattice' missing. Aborting..." << std::endl;
                std::terminate();
            }

            if(!_config.contains("const dims")) {
                std::cerr << "[Json Error]: field 'const dims' missing. Aborting..." << std::endl;
                std::terminate();
            }

            const Eigen::Index time_step = _config["time step"];
            const Eigen::Index dims = _config["dims"];
            const Eigen::Index K = _config["K"];
            const Eigen::Vector3i cardinal = extract3i(_config["cardinal"]);
            const Eigen::Vector3i lattice = extract3i(_config["lattice"]);
            const Eigen::Vector3<T> aabb_min = extract3(_config["aabb min"]);
            const Eigen::Vector3<T> aabb_max = extract3(_config["aabb max"]);
            const Eigen::VectorX<T> const_dims = extract(_config["const dims"]);

            //output file (optional)
            std::filesystem::path output;
            if(!_config.contains("output")) {
                std::stringstream ss;
                ss << "out_" << input.stem().c_str() << "_" << cardinal(0) << "_" << cardinal(1) << "_" << cardinal(2);
                output = input.parent_path() / ss.str();
            } else output = std::filesystem::path(_config["output"]);

            TurboDorn::IO::Detail::File<T> file = TurboDorn::IO::Detail::load_from_file<T, size_t>(input, dims, K, TurboDorn::Hagedorn::Detail::hyperbolic_cut_shape);

            const auto inv = Hagedorn::Detail::prepare_invariants_from_file<T>(file);
            TurboDorn::Detail::sample_grid<T>(output, time_step, cardinal, lattice, aabb_min, aabb_max, const_dims, *inv);
        }
    };

    //---------------------------------------------------------------------------------------//
    //                                        Renderer
    //---------------------------------------------------------------------------------------//
    /*
        An example implementation that shows how to render a give file from the config json.
    */
    template<class T>
    class Renderer {
    public:
        Renderer(const json& _views) {

            const auto extract3 = [](const json& _config) -> Eigen::Vector3<T> {
                Eigen::Vector3<T> out;
                Eigen::Index i = 0;
                for(auto it = _config.begin(); it != _config.end(); ++it){
                    out.coeffRef(i++) = T(*it);
                }
                return out;
            };

            const auto extract2i = [](const json& _config) -> Eigen::Vector2i {
                Eigen::Vector2i out;
                Eigen::Index i = 0;
                for(auto it = _config.begin(); it != _config.end(); ++it){
                    out.coeffRef(i++) = int(*it);
                }
                return out;
            };

            size_t idx = 0;

            for(auto it = _views.begin(); it != _views.end(); ++it){

                const json& j = (*it);

                const std::filesystem::path input = std::filesystem::path(j["input"]);
                if(!std::filesystem::exists(input)){
                    std::cerr << "[IO Error]: Input does not exist. Aborting..." << std::endl;
                    std::terminate();
                }

                if(!j.contains("camera position")) {
                    std::cerr << "[Json Error]: field 'camera position' missing. Aborting..." << std::endl;
                    std::terminate();
                }

                if(!j.contains("camera forward")) {
                    std::cerr << "[Json Error]: field 'camera forward' missing. Aborting..." << std::endl;
                    std::terminate();
                }

                if(!j.contains("camera up")) {
                    std::cerr << "[Json Error]: field 'camera up' missing. Aborting..." << std::endl;
                    std::terminate();
                }

                if(!j.contains("image pixel")) {
                    std::cerr << "[Json Error]: field 'image pixel' missing. Aborting..." << std::endl;
                    std::terminate();
                }

                if(!j.contains("pixel world size")) {
                    std::cerr << "[Json Error]: field 'pixel world size' missing. Aborting..." << std::endl;
                    std::terminate();
                }

                if(!j.contains("aabb min")) {
                    std::cerr << "[Json Error]: field 'aabb min' missing. Aborting..." << std::endl;
                    std::terminate();
                }

                if(!j.contains("aabb max")) {
                    std::cerr << "[Json Error]: field 'aabb max' missing. Aborting..." << std::endl;
                    std::terminate();
                }

                if(!j.contains("step size")) {
                    std::cerr << "[Json Error]: field 'step size' missing. Aborting..." << std::endl;
                    std::terminate();
                }

                const Eigen::Vector3<T> cam_or = extract3(j["camera position"]);
                const Eigen::Vector3<T> cam_dir = extract3(j["camera forward"]);
                const Eigen::Vector3<T> cam_right = extract3(j["camera up"]);
                const Eigen::Vector2i img_size_pixel = extract2i(j["image pixel"]);
                const T pixel_world_size = j["pixel world size"];
                const Eigen::Vector3<T> aabb_min = extract3(j["aabb min"]);
                const Eigen::Vector3<T> aabb_max = extract3(j["aabb max"]);
                const T step_size = j["step size"];

                //optional
                const T hsl_max_value = j.contains("hsl max value") ? T(j["hsl max value"]) : T(10.);
                const size_t max_gridsize_byte = j.contains("max grid size byte") ? size_t(j["max grid size byte"]) : 500000000;

                //output file (optional)
                std::filesystem::path output;
                if(!j.contains("output")) {
                    std::stringstream ss;
                    ss << "out_" << idx++ << "_" << input.stem().c_str() << ".png";
                    output = input.parent_path() / ss.str();
                } else output = std::filesystem::path(j["output"]);

                Geometry::RenderGrid<T> grid(input, max_gridsize_byte);

                std::cout << "------------ Render start ------------" << std::endl;
                const auto start = std::chrono::high_resolution_clock::now();
                const std::vector<T> buffer = TurboDorn::Detail::render_grid<T>(
                    input, 
                    cam_or, 
                    cam_dir, 
                    cam_right, 
                    img_size_pixel, 
                    pixel_world_size, 
                    aabb_min, 
                    aabb_max, 
                    step_size, 
                    hsl_max_value, 
                    grid
                );
                const auto end = std::chrono::high_resolution_clock::now();
                std::cout << "------------ Render done after " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "s ------------\n";
                const std::vector<unsigned char> img = TurboDorn::Detail::normalize_to_char(buffer);

                if(lodepng::encode(output.c_str(), img.data(), img_size_pixel(0), img_size_pixel(1))){
                    std::cerr << "failed to save img" << std::endl;
                }

            }

        }
    };//Renderer

}//Hagedorn

#endif //TURBODORN_HPP
