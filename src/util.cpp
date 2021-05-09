
#include "../include/util.hpp"

FilePathResolver* FilePathResolver::instance = new FilePathResolver();

void FilePathResolver::resolve() {
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
        else throw new TVexception("folder shader not found");
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
        else throw new TVexception("folder asset not found");
    }

}