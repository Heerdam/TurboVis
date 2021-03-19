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

}


#endif /* DATA_HPP */