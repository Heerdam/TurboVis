#include "../include/data.hpp"

#define __AVX2__
#include <libmorton/morton.h>

using namespace Data;

SparseGrid::SparseGrid(){}

void SparseGrid::init(Gridsize _newSize) {
    size = _newSize;
    map.reserve(size*size*size);
}

void SparseGrid::downsample(Gridsize _newSize){
    //TODO
}

Gridsize SparseGrid::getSize() const {
    return size;
}

size_t SparseGrid::getSizeOfMapInBytes() const {
    return sizeof(float) * map.size();
}

void SparseGrid::insert(float _val, uint16_t _x, uint16_t _y, uint16_t _z){
    const size_t index = libmorton::morton3D_64_encode((uint_fast32_t)_x, (uint_fast32_t)_y, (uint_fast32_t)_z);
    map[index] = _val;
}

void SparseGrid::insertBulk(Iterator _begin, Iterator _end){
    for(Iterator i = _begin; i < _end; ++i){
        const auto& t = *i;
        const size_t index = libmorton::morton3D_64_encode((uint_fast32_t)std::get<1>(t), (uint_fast32_t)std::get<2>(t), (uint_fast32_t)std::get<3>(t));
        map[index] = std::get<0>(t);
    }
}

void SparseGrid::insert_async(float _val, uint16_t _x, uint16_t _y, uint16_t _z){
    std::lock_guard<std::mutex> lock(mutex);
    insert(_val, _x, _y, _z);
}

void SparseGrid::insertBulk_async(Iterator _begin, Iterator _end){
    std::lock_guard<std::mutex> lock(mutex);
    insertBulk(_begin, _end);
}

std::vector<float> SparseGrid::toBuffer() const{
    const size_t bs = static_cast<size_t>(size);
    std::vector<float> buffer(bs * bs * bs);
    for(auto it = map.begin(); it != map.end(); ++it){
        const size_t index = (*it).first;
        uint_fast32_t x, y, z;
        libmorton::morton3D_64_decode(index, x, y, z);
        buffer[z*size*size + y*size +x] = (*it).second;
    }
    return buffer;
}

std::vector<float> SparseGrid::createBuffers(float _cellwidth) const{
    const size_t bs = static_cast<size_t>(size);
    std::vector<float> v(bs * bs * bs * 3);
    for(uint16_t z = 0; z < bs; ++z){
        for(uint16_t y = 0; y < bs; ++y){
            for(uint16_t x = 0; x < bs; ++x){
                const size_t index = 3*(z*bs*bs + y*bs + x);
                v[index] = x * _cellwidth; //x
                v[index + 1] = y * _cellwidth; //y
                v[index + 2] = z * _cellwidth; //z
            }
        }
    }
    return v;
}

std::vector<float> FunctionData::createBuffers() const {
    return grid.createBuffers(25.f);
}

// [start, end + 1)
void FunctionData::eval(
    std::function<float(int16_t, int16_t, int16_t)> _func, 
    SparseGrid& _grid,
    bool _async,
    int16_t _startX, int16_t _endX, 
    int16_t _startY, int16_t _endY, 
    int16_t _startZ, int16_t _endZ)
{
    std::vector<Triplet> tmp;
    tmp.reserve((_endX - _startX)*(_endY - _startY)*(_endZ - _startZ));
    for(int16_t z = _startZ; z < _endZ; ++z){
        for(int16_t y = _startY; y < _endY; ++y){
            for(int16_t x = _startX; x < _endX; ++x){
                const float v = _func(x, y, z);
                if(std::abs(v) > std::numeric_limits<float>::epsilon()){
                    tmp.emplace_back(v, x, y, z);
                }
            }
        }
    }

    if(_async)
        _grid.insertBulk_async(tmp.begin(), tmp.end());
    else
        _grid.insertBulk(tmp.begin(), tmp.end());

}

void FunctionData::createDataFromFunction(std::function<float(int16_t, int16_t, int16_t)> _func, Gridsize _size){
    grid.init(_size);
    switch(_size){
        case Gridsize::x2:
        {        
            eval(_func, grid, false, 0, 2, 0, 2, 0, 2);
        }
        break;
        case Gridsize::x4:
        {
            eval(_func, grid, false, 0, 4, 0, 4, 0, 4);
        }
        break;
        case Gridsize::x16:
        {
            /*
            std::vector<Eigen::Triplet<float>> data;
            data.reserve(16*16*16);
            std::vector<std::future<void>> futures;
            for(int16_t z = 0; z < 4; ++z){
                for(int16_t y = 0; y < 4; ++y){
                    for(int16_t x = 0; x < 4; ++x){
                        futures.push_back(std::async(std::launch::async, &FunctionData::eval, _func, std::ref(data), x * 4, x * 4 + 4, y * 4, y * 4 + 4, z * 4, z * 4 + 4));
                    }
                }           
            }
            for(auto& f : futures){
                f.wait();
            }
            vals.setFromTriplets(data.begin(), data.end());
            vals.makeCompressed();
            */
        }
        break;
        case Gridsize::x32:
        {

        }
        break;
        case Gridsize::x64:
        {

        }
        break;
        case Gridsize::x128:
        {

        }
        break;
        case Gridsize::x256:
        {

        }
        break;
        case Gridsize::x512:
        {

        }
        break;
        case Gridsize::x1024:
        {

        }
        break;
        case Gridsize::x2048:
        {

        }
        break;
        case Gridsize::x4096:
        {

        }
        break;
    }
}

void FunctionData::downsample(Gridsize _newSize){
    grid.downsample(_newSize);
}

std::vector<float> FunctionData::toBuffer() const{
    return grid.toBuffer();
}
    
size_t FunctionData::getDataSizeinBytes() const{
    return grid.getSizeOfMapInBytes();
}
    
Gridsize FunctionData::currentGridSize() const{
    return grid.getSize();
}