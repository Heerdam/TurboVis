
#include "../include/util.hpp"

FilePathResolver* FilePathResolver::instance = new FilePathResolver();

void FilePathResolver::resolve() {

    const auto path = std::filesystem::current_path() /= "shader/";
    const auto path1 = std::filesystem::current_path() /= "../shader/";
    const auto path2 = std::filesystem::current_path() /= "../../shader/";

    if(std::filesystem::exists(path)){
        instance->shader_dir = path.string();
        return;
    }

    if(std::filesystem::exists(path1)){
        instance->shader_dir = path1.string();
        return;
    }

    if(std::filesystem::exists(path2)){
        instance->shader_dir = path2.string();
        return;
    }

}