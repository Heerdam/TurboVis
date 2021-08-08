#ifndef UTIL_HPP
#define UTIL_HPP

#include <cstring>
#include <fstream>
#include <filesystem>
#include <stdexcept>

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

inline FilePathResolver* FilePathResolver::instance = new FilePathResolver();

inline void FilePathResolver::resolve() {
    //shader
    {
        const auto s_path_1 = std::filesystem::current_path() /= "shader/";
        const auto s_path_2 = std::filesystem::current_path() /= "../shader/";
        const auto s_path_3 = std::filesystem::current_path() /= "../../shader/";

        if(std::filesystem::exists(s_path_1))
            instance->shader_dir = s_path_1.string();
        else if(std::filesystem::exists(s_path_2))
            instance->shader_dir = s_path_2.string();
        else if(std::filesystem::exists(s_path_3))
            instance->shader_dir = s_path_3.string();
        else throw std::runtime_error("folder shader not found");
    }
    //asset
    {
        const auto s_path_1 = std::filesystem::current_path() /= "assets/";
        const auto s_path_2 = std::filesystem::current_path() /= "../assets/";
        const auto s_path_3 = std::filesystem::current_path() /= "../../assets/";

        if(std::filesystem::exists(s_path_1))
            instance->asset_dir = s_path_1.string();
        else if(std::filesystem::exists(s_path_2))
            instance->asset_dir = s_path_2.string();
        else if(std::filesystem::exists(s_path_3))
            instance->asset_dir = s_path_3.string();
        else throw std::runtime_error("folder asset not found");
    }

}
#endif /* UTIL_HPP */
