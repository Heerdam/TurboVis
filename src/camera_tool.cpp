
#include <TurboDorn/turbodorn.hpp>

int main() {

    constexpr size_t its = 36;

    json file;
    file["mode"] = "render";
    file["views"] = json::array();

    const Eigen::Vector3d center = Eigen::Vector3d::Zero();
    const double dPi = (2. * M_PI) / double(its);

    for(size_t i = 0; i < its; ++i){

        json t;
        t["input"] = "../example_files/out_simulation_results_phi_0_1_2";
        {
            std::stringstream ss;
            ss << "../example_files/out_" << i << ".png";
            t["output"] = ss.str();
        }
        t["pixel world size"] = 5;
        t["step size"] = 0.1;

        t["image pixel"] = json::array();
        t["image pixel"].push_back(500);
        t["image pixel"].push_back(500);

        t["aabb min"] = json::array();
        t["aabb min"].push_back(-10);
        t["aabb min"].push_back(-10);
        t["aabb min"].push_back(-10);

        t["aabb max"] = json::array();
        t["aabb max"].push_back(10);
        t["aabb max"].push_back(10);
        t["aabb max"].push_back(10);

        const double ax = std::cos(i * dPi);
        const double ay = std::sin(i * dPi);

        const Eigen::Vector3d dir (ax, ay, 0.);
        const Eigen::Vector3d cpos = center + dir * 15.;

        t["camera position"] = json::array();
        t["camera position"].push_back(cpos(0));
        t["camera position"].push_back(cpos(1));
        t["camera position"].push_back(cpos(2));

        t["camera forward"] = json::array();
        t["camera forward"].push_back(-dir(0));
        t["camera forward"].push_back(-dir(1));
        t["camera forward"].push_back(-dir(2));
        
        t["camera up"] = json::array();
        t["camera up"].push_back(0);
        t["camera up"].push_back(0);
        t["camera up"].push_back(1);

        file["views"].push_back(t);

    }

    std::ofstream out ("cam.json");
    out << file.dump(4);

    return EXIT_SUCCESS;

}
