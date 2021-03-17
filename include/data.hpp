#ifndef DATA_HPP
#define DATA_HPP

#include "defines.hpp"

using Triplet = std::tuple<float, uint16_t, uint16_t, uint16_t>;

enum class Gridsize {
    x2, x4, x16, x32, x64, x128,
    x256, x512, x1024, x2048, x4096
};

inline size_t GET_WIDTH(Gridsize _size){
    switch(_size){
        case Gridsize::x2: return 2;
        case Gridsize::x4: return 4;
        case Gridsize::x16: return 16;
        case Gridsize::x32: return 32;
        case Gridsize::x64: return 64;
        case Gridsize::x128: return 128;
        case Gridsize::x256: return 256;
        case Gridsize::x512: return 512;
        case Gridsize::x1024: return 1024;
        case Gridsize::x2048: return 2048;
        case Gridsize::x4096: return 4096;
        default: return 0;
    }
}

inline size_t operator*(const uint_fast32_t _index, const Gridsize _size){
    return _index * GET_WIDTH(_size);
}

inline size_t operator*(const size_t _index, const Gridsize _size){
    return _index * GET_WIDTH(_size);
}

inline size_t operator*(const Gridsize _s1, const Gridsize _s2){
    return  GET_WIDTH(_s1) *  GET_WIDTH(_s2);
}

class SparseGrid {

    std::mutex mutex;
    std::unordered_map<size_t, float> map;
    Gridsize size;

public:
    SparseGrid();
    void init(const Gridsize&);
    void downsample(const Gridsize&);
    Gridsize getSize() const;
    size_t getSizeOfMapInBytes() const;

    void insert(float, uint16_t, uint16_t, uint16_t);
    void insertBulk(const std::vector<Triplet>&);

    void insert_async(float, uint16_t, uint16_t, uint16_t);
    void insertBulk_async(const std::vector<Triplet>&);

    std::unique_ptr<float> toBuffer() const;
    std::unique_ptr<float> createBuffers(const Gridsize&, const float) const;
};

class FunctionData {

    std::mutex mutex;
    SparseGrid grid;

    static void eval(const std::function<float(int16_t, int16_t, int16_t)>, SparseGrid&, const bool, const int16_t, const int16_t, const int16_t, const int16_t, const int16_t, const int16_t);
    
public:

    void createDataFromFunction(const std::function<float(int16_t, int16_t, int16_t)>, const Gridsize, const size_t);
    void downsample(const Gridsize&);

    std::unique_ptr<float> toBuffer() const;
    size_t getDataSizeinBytes() const;
    Gridsize currentGridSize() const;

};


#endif /* DATA_HPP */