
//#define DEBUG_TRACE
#include <TurboDorn/turbodorn.hpp>

/*
    //---------------------------------------------------------------------------------------//
    //                                        TURBOVIS
    //---------------------------------------------------------------------------------------//

    ******************************** How to use this example implementation *******************************
    to start the program: tv path/to/config/file
    or just tv and it assumes config.json in the same folder as the binary

    Sample Mode: Samples the provided function in hdf5 form into binary form.
    Config Json valid options for build mode:
    {
        "mode": "sample",
        "input": "path/to/input/file.hdf5",
        "output": "path/to/output/file", //optional. appends samples_[input file name]_[dim1]-[dim2]-[dim3]

        "time step": 0, //the time step of the simulation [0, 1, ...]
        "dims": 3, //total dimensions
        "K": 3,
        "cardinal": [0, 1, 2], //index of dimensions to be used
        "aabb min": [x, y, z], //lower position
        "aabb max": [x, y, z], //upper position
        "lattice": [x, y, z], //voxel resolution of the cube
        "const dims": [val, val,...], //constant values of the unused dimensions
    }

    Render Mode: Renders the provided function in binary form with the provided parameters in json to a png
    Config Json valid options for render mode:
    {
        "mode": "render",
        "input": "path/to/input/file.hdf5",
        "output": "path/to/output/file.png", //optional
        "camera position": [x, y, z], 
        "camera forward": [x, y, z], 
        "camera up": [x, y, z], 
        "pixels": [x, y],
        "pixel world size": 5, //the scaling of a pixel to world unit
        "aabb min": [x, y, z], 
        "aabb max": [x, y, z],
        "step size": 0.1 //size of a single step for the ray marcher
    }
*/

int main(int argc, char* argv[]) {

    //std::string tt = "../example_files/cam_exp.json";
    //char* targv[] = { nullptr, tt.data() };

    const auto getJsonPath = [](int _argc, char* _argv[]) {
        if(_argc == 1){
            const auto out = std::filesystem::current_path() / "config.json";
            if(!std::filesystem::exists(out) || out.extension() == "json") {
                std::cerr << "file error: no config.json in working directory found" << std::endl;
                std::terminate();
            }
            return out;
        } else {
            const auto out = std::filesystem::current_path() / _argv[1]; 
            if(!std::filesystem::exists(out) || out.extension() == "json"){
                std::cerr << "file error: no .json found at: " << _argv[1] << std::endl;
                std::terminate();
            }
            return out;
        }
    };

    const std::filesystem::path p_config = getJsonPath(2, argv);
    std::ifstream f_config(p_config);
    if(!f_config.good()){
        std::cerr << "file error: error reading provided config file" << std::endl;
        std::terminate();
    }
    const json config = json::parse(f_config);
    
    //parse mode: 0: invalid, 1: sample, 2: render
    const std::string mk = config["mode"];
    const int mode = mk == "sample" ? 1 : mk == "render" ? 2 : 0;
    if(mode == 0){
        std::cerr << "json error: key mode is missing or wrong value" << std::endl;
        std::terminate();
    }

    /*
    ***********************************************************************************************************
        If we have reached this point the config file is valid and we can proceed with sampling or rendering.
    ***********************************************************************************************************
    */
    //build
    if(mode == 1){
        TurboDorn::Sampler<double> sampler (config);
    }
    //render
    if(mode == 2){
        TurboDorn::Renderer<double> renderer (config["views"]);
    }
    return EXIT_SUCCESS;
}
