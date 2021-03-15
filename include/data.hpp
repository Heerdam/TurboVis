#ifndef DATA_HPP
#define DATA_HPP

#include "defines.hpp"

class Chunk {

};

class ChunkGrid {

public:
    
    std::vector<Chunk*> getChunksForFrustum(const Frustum&) const;
};


#endif /* DATA_HPP */