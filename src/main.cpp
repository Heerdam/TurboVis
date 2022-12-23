
#include <TurboDorn/turbodorn.hpp>

/*
    //---------------------------------------------------------------------------------------//
    //                                        TURBOVIS
    //---------------------------------------------------------------------------------------//






    ********************************** How to build TurboVis **********************************





    ******************************** How to use and run TurboVis *******************************
    to start the program: tv -c path/to/config/file
    or just tv and it assumes config.json in the same folder as the binary

    Sample Mode: Samples the provided function in hdf5 form into binary form.
    Config Json valid options for build mode:
    {
        "mode": "sample",
        "input": "path/to/input/file.hdf5",
        "output": "path/to/output/file.hdf5", //optional. appends out_[input]_[dim1]-[dim2]-[dim3].hdf

        "aabb": [[x, y, z], [x, y, z]], //upper and lower position
        "resolution": [x, y, z], //voxel resolution of the cube
        "dimensions": [0, 1, 2], //index of dimensions to be used
        "dim_val": [val, val,...], //constant values of the unused dimensions
    }

    Render Mode: Renders the provided function in binary form with the provided parameters in json to a png
    Config Json valid options for render mode:
    {
        "mode": "render",
        "input": "path/to/input/file.hdf5",
        "output": "path/to/output/file.hdf5", //optional. appends out_[input].jpg

        "pixels": [x, y],
        "camera": [[x, y, z], [x, y, z]], //dir, up
        "steps": X, //step size will be diagonale of the previously defined AABB / steps
    }

    ********************************** Example Json Sample **********************************
    {
        "mode": "build",
        "input": "path/to/input/file.hdf5",
        "output": "path/to/output/file.hdf5", //optional. appends out_[input].hdf

        "aabb": [[x, y, z], [x, y, z]], //upper and lower position
        "resolution": [x, y, z], //voxel resolution of the cube
    }
    ****************************************************************************************

    ********************************** Example Json Render *********************************
    {
        "mode": "render",
        "input": "path/to/input/file.hdf5",
        "output": "path/to/output/file.hdf5", //optional. appends out_[input].jpg

        "pixels": [x, y],
        "camera": [[x, y, z], [x, y, z]], //dir, up
        "steps": X, //step size will be diagonale of the previously defined AABB / steps
    }
    ****************************************************************************************

*/

int main(int argc, char* argv[]) {

    const auto getJsonPath = [&] {
        if(argc == 0){
            const auto out = std::filesystem::current_path() / "config.json";
            if(!std::filesystem::exists(out) || out.extension() == "json") {
                std::cerr << "file error: no config.json in working directory found" << std::endl;
                std::terminate();
            }
            return out;
        } else {
            const auto out = std::filesystem::current_path() / argv[0]; 
            if(!std::filesystem::exists(out) || out.extension() == "json"){
                std::cerr << "file error: no .json found at: " << argv[0] << std::endl;
                std::terminate();
            }
            return out;
        }
    };

    /*
    ***********************************************************************************************************
        We are doing only preliminary validation of the config file here and if a valid input file exists.
        All the specific extraction and validation of the config file is done in the build and render modes.
    ***********************************************************************************************************
    */

    const std::filesystem::path p_config = getJsonPath();
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

    //input file (must exist)
    const std::filesystem::path input = config["input"];
    if(!std::filesystem::exists(input)){
        std::cerr << "json error: key mode is missing or wrong value" << std::endl;
        std::terminate();
    }

    //output file (optional)
    std::filesystem::path output = config["output"];
    if(output.empty()) {
        std::stringstream ss;
        ss << input.parent_path() << "/out_" << input.filename();
        output = ss.str();
    }

    /*
    ***********************************************************************************************************
        If we have reached this point the config file is valid and we can proceed with sampling or rendering.
    ***********************************************************************************************************
    */

    //build
    if(mode == 1){

        TurboDorn::Sampler sampler (config);

        if(sampler.success()) EXIT_SUCCESS;
        else EXIT_FAILURE;
    }
    
    //render
    if(mode == 2){
        TurboDorn::Renderer renderer (config);
        
        if(renderer.success()) EXIT_SUCCESS;
        else EXIT_FAILURE;
    }

    return EXIT_SUCCESS;

}













/*

enum class File : size_t {
    none = 0,
    simulation_results = 1,
    simulation_results_phi000 = 2,
    simulation_results_phi100 = 3,
    simulation_results_phi121 = 4,
    simulation_results_phi412 = 5,
};
    if(argc == 1){
        std::cout << "no valid input" << std::endl;
        return 0;
    }

    std::vector<std::string_view> inputs;
    for(size_t i = 0; i < argc; ++i)
        inputs.emplace_back(argv[i]);

    //-v [file] [grid size] [path]
    if(inputs[1] == "-v"){

        File file = File::none;
        size_t gs = 1000;
        std::filesystem::path path_in = std::filesystem::current_path();
        std::filesystem::path path_out = std::filesystem::current_path();

        //parse file
        File file = File::none;
        if(inputs[1] == "phi") file = File::simulation_results;
        else if(inputs[1] == "phi000") file = File::simulation_results_phi000;
        else if(inputs[1] == "phi100") file = File::simulation_results_phi100;
        else if(inputs[1] == "phi121") file = File::simulation_results_phi121;
        else if(inputs[1] == "phi412") file = File::simulation_results_phi412;
        else {
            std::cout << "no valid input" << std::endl;
            return 0;
        }

        //grid size
        if(inputs.size() < 3){
            std::cout << "no valid input" << std::endl;
            return 0;
        } else {
            gs = std::stoul(inputs[3].data());
        }

        //path
        if(inputs.size() >= 4){
            path_out = inputs[3];
        }

        TurboDorn::Voxeliser voxel;
        voxel.compute(path_in, path_out);

    } else if(inputs[1] == "-r"){

    } else {
        std::cout << "no valid input" << std::endl;
        return 0;
    }

}
*/
/*

int main() {

    FilePathResolver::resolve();

    

    struct Work{
        std::string name;
        size_t width = 1500, height = 1500;
        File file;
        //camera
        Vec3 position, target;
        Vec3 upAxis;
        Eigen::Matrix<double, -1, 1> lower, upper;
        double scale = 1.;
        double maxDist = 35.;
        size_t steps = 1000;
        size_t curT = 0;
        double MAX = 10.f;
    };

    const auto getFile = [](File _file){
        switch(_file){
            case File::simulation_results: return IO::simulation_results();
            case File::simulation_results_phi000: return IO::simulation_results_phi000();
            case File::simulation_results_phi100: return IO::simulation_results_phi100();
            case File::simulation_results_phi121: return IO::simulation_results_phi121();
            case File::simulation_results_phi412: return IO::simulation_results_phi412();
            default: return std::optional<IO::File<double>>();
        }
    };

    std::queue<Work> todo;
    if constexpr(true){
        Work w;
        w.width = 1500;
        w.height = 1500;
        w.steps = 1500;
        w.name = "sm_000_z_0.png";
        w.file = File::simulation_results_phi000;
        w.position = Vec3(0.f, 0.f, -25.f);
        w.target = Vec3(0.f, 0.f, 0.f);
        w.upAxis = Vec3(0.f, 1.f, 0.f);
        w.lower.resize(3);
        w.upper.resize(3);
        w.lower(0) = -100.;
        w.lower(1) = -100.;
        w.lower(2) = -24.;
        w.upper(0) = 100.;
        w.upper(1) = 100.;
        w.upper(2) = 100.;
        w.scale = 4.; 
        w.maxDist = 30.;
        w.curT = 0;

        todo.push(std::move(w));
    }    
    if constexpr(true){
        Work w;
        w.width = 1500;
        w.height = 1500;
        w.steps = 1500;
        w.name = "sm_000_z_1.png";
        w.file = File::simulation_results_phi000;
        w.position = Vec3(0.f, 0.f, -25.f);
        w.target = Vec3(0.f, 0.f, 0.f);
        w.upAxis = Vec3(0.f, 1.f, 0.f);
        w.lower.resize(3);
        w.upper.resize(3);
        w.lower(0) = -100.;
        w.lower(1) = -100.;
        w.lower(2) = -24.;
        w.upper(0) = 100.;
        w.upper(1) = 100.;
        w.upper(2) = 100.;
        w.scale = 4.; 
        w.maxDist = 30.;
        w.curT = 1;

        todo.push(std::move(w));
    }
    if constexpr(true){
        Work w;
        w.width = 1500;
        w.height = 1500;
        w.steps = 1500;
        w.name = "sm_000_z_2.png";
        w.file = File::simulation_results_phi000;
        w.position = Vec3(1.f, 0.f, -25.f);
        w.target = Vec3(1.f, 0.f, 0.f);
        w.upAxis = Vec3(0.f, 1.f, 0.f);
        w.lower.resize(3);
        w.upper.resize(3);
        w.lower(0) = -100.;
        w.lower(1) = -100.;
        w.lower(2) = -24.;
        w.upper(0) = 100.;
        w.upper(1) = 100.;
        w.upper(2) = 100.;
        w.scale = 4.; 
        w.maxDist = 30.;
        w.curT = 2;

        todo.push(std::move(w));
    }
    if constexpr(true){
        Work w;
        w.width = 1500;
        w.height = 1500;
        w.steps = 1500;
        w.name = "sm_000_z_3.png";
        w.file = File::simulation_results_phi000;
        w.position = Vec3(1.f, 0.f, -25.f);
        w.target = Vec3(1.f, 0.f, 0.f);
        w.upAxis = Vec3(0.f, 1.f, 0.f);
        w.lower.resize(3);
        w.upper.resize(3);
        w.lower(0) = -100.;
        w.lower(1) = -100.;
        w.lower(2) = -24.;
        w.upper(0) = 100.;
        w.upper(1) = 100.;
        w.upper(2) = 100.;
        w.scale = 4.; 
        w.maxDist = 30.;
        w.curT = 3;

        todo.push(std::move(w));
    }

    const auto start = std::chrono::high_resolution_clock::now();
    while(!todo.empty()){
        const Work w = todo.front();
        todo.pop();

        GL::Camera cam(int64_t(w.width), int64_t(w.height), glm::radians(55.f), 0.01f, 100.f);
        cam.position = w.position;
        cam.target = w.target;
        cam.upAxis = w.upAxis;
        cam.update();

        std::vector<double> buffer;
        buffer.resize(w.width * w.height * 4);

        //render
        const auto file = getFile(w.file);
        if(!file.has_value()){
            std::cout << "could not load file: " << w.name << std::endl;
            continue;
        }
        const auto inv = Math::Hagedorn::computeInvariants(file.value());

        AsyncRenderer renderer;
        renderer.start(64, w.width, w.height, [&](size_t _threat_index, size_t _x, size_t _y){

            using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;

            const Eigen::Matrix<double, 2, 1> bounds = Eigen::Matrix<double, 2, 1>(double(w.width), double(w.height));
            const Eigen::Matrix<double, 2, 1> uv = Eigen::Matrix<double, 2, 1>(double(_x), double(_y));
            const Eigen::Matrix<double, 2, 1> pp = (2. * uv - bounds) / bounds(1);  //eye ray

            const Vec3 r_o_glm = cam.position + w.scale*(pp(0) * cam.right + pp(1) * cam.up );
            const Vec3 r_d_glm = cam.dir;

            Vector r_o (3);
            r_o(0) = r_o_glm.x;
            r_o(1) = r_o_glm.y;
            r_o(2) = r_o_glm.z;

            Vector r_d (3);
            r_d(0) = r_d_glm.x;
            r_d(1) = r_d_glm.y;
            r_d(2) = r_d_glm.z;

            const double dS = w.maxDist / double(w.steps);

            Eigen::Matrix<double, 3, 1> col;
            col.setZero();

            double t = 0.;
            double transmission = 1.;
            
            //intersect bounding box for early out
            const bool hit = Math::intersect(r_o, r_d, w.lower, w.upper, w.maxDist, t); 
            t +=dS;
            if(!hit){   

                col(0) = col(1) = col(2) = 0.;

            } else if(w.steps == 1){   

                col(0) = 1.;
                transmission = 0.;
                
            } else {
    
                for (size_t s = 0; s < w.steps; ++s) {
                    const Vector pos = r_o + t * r_d;
                    t += dS;

                    if(t > w.maxDist || 
                        pos(0) < w.lower(0) || 
                        pos(1) < w.lower(1) ||
                        pos(2) < w.lower(2) ||
                        pos(0) > w.upper(0) ||
                        pos(1) > w.upper(1) ||
                        pos(2) > w.upper(2)
                    )
                        break;

                    //calculate basis function
                    const std::unordered_map<Eigen::Index, std::complex<double>> phis = Math::Hagedorn::compute(
                        w.curT,
                        pos,
                        inv
                    );

                    std::complex<double> res (0., 0.);
                    for(Eigen::Index k = 0; k < file.value().Ks.size(); ++k){ 
                        const Eigen::Index idx = Math::Hagedorn::Detail::index(file.value().Ks[k], file.value().k_max);
                        assert(phis.contains(idx));
                        const std::complex<double>& p = (*phis.find(idx)).second;
                        res += file.value().c_0[w.curT](k) * p;
                    } 
                    res *= std::exp(std::complex<double>(0., 1.) * file.value().S(Eigen::Index(w.curT)));

                    //compute color
                    const auto hsl = Math::c_to_HSL(double(w.MAX), res);
                    transmission *= std::exp(-hsl(2) * dS);
                    const auto rgb = Math::HSL_to_RGB_rad(hsl);
                    col += transmission * rgb * dS;

                    if(renderer.isShutdown() || renderer.isRestart(_threat_index))
                        return;
                }
            }        

            const size_t idx = _y * cam.width + _x;
            buffer[4 * idx] = col(0);//uint16_t(std::round(col(0) * 255.));
            buffer[4 * idx + 1] = col(1);//uint16_t(std::round(col(1) * 255.));
            buffer[4 * idx + 2] = col(2);//uint16_t(std::round(col(2) * 255.));;
            buffer[4 * idx + 3] = (1. - transmission);//uint16_t(std::round((1. - transmission) * 255.));

        });

        const auto fstart = std::chrono::high_resolution_clock::now();
        std::cout << "started file: " << w.name << std::endl;
        std::cout << "Progress: [";
        std::cout.flush();
        double lastProg = 0.;
        while(true){
            if(renderer.getProgress() >= 1.) break;
            const double p = renderer.getProgress();
            if((p - lastProg) >= 0.1){
                std::cout << "x";
                std::cout.flush();
                lastProg = p;
            } 
            std::this_thread::sleep_for(5s);
        }
        std::cout << "]" << std::endl;
        renderer.stop();
        
        const double ftime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - fstart).count();
        std::cout << "finished file: " << w.name << " [" << ftime / (60.) << "mins]" << std::endl;

        //save file
        std::vector<unsigned char> out (buffer.size());
        double maxVal = 0.f;
        for(size_t i = 0; i < buffer.size(); ++i){
            maxVal = std::max(maxVal, buffer[i]);
        }
        for(size_t i = 0; i < buffer.size(); ++i){
            out[i] = uint16_t(std::round((buffer[i]/maxVal) * 255.));
        }
        const uint32_t error = lodepng::encode(w.name, out, uint32_t(w.width), uint32_t(w.height));
        if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
        else std::cout << "file successfully saved [" << w.name << "]" << std::endl;
    }

    const double time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
    std::cout << "ALl finished after " << time / (60.) << "mins" << std::endl;
    std::cout << "all done. goodbye" << std::endl;

}

*/