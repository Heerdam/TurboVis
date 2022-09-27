#ifndef FILEPATHRESOLVER_HPP
#define FILEPATHRESOLVER_HPP

#include <cstring>
#include <fstream>
#include <filesystem>
#include <stdexcept>

class FilePathResolver {

    std::filesystem::path example_files;

public:

    FilePathResolver() {

        const auto s_path_1 = std::filesystem::current_path() / "example_files/";
        const auto s_path_2 = std::filesystem::current_path() / "../example_files/";
        const auto s_path_3 = std::filesystem::current_path() / "../../example_files/";

        if(std::filesystem::exists(s_path_1))
            example_files = s_path_1.string();
        else if(std::filesystem::exists(s_path_2))
            example_files = s_path_2.string();
        else if(std::filesystem::exists(s_path_3))
            example_files = s_path_3.string();
        else throw std::runtime_error("folder example_files not found");

    }

    std::filesystem::path example_files_path() noexcept {
        return example_files;
    }

};

#endif /* UTIL_HPP */