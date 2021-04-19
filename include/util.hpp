#ifndef UTIL_HPP
#define UTIL_HPP

#include "defines.hpp"

class FilePathResolver {

    std::string shader_dir;

    static FilePathResolver* instance;
    FilePathResolver(){};

public:
    static void resolve();

    static [[nodiscard]] std::string SHADERDIR() noexcept {
        return instance->shader_dir;
    }

};

#endif /* UTIL_HPP */
