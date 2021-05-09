#ifndef UTIL_HPP
#define UTIL_HPP

#include "defines.hpp"

class FilePathResolver {

    std::string shader_dir;
    std::string asset_dir;

    static FilePathResolver* instance;
    FilePathResolver(){};

public:
    static void resolve();

    [[nodiscard]] static std::string SHADERDIR() noexcept {
        return instance->shader_dir;
    }

    [[nodiscard]] static std::string ASSETDIR() noexcept {
        return instance->asset_dir;
    }

};

#endif /* UTIL_HPP */
