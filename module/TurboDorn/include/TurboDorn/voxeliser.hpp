#ifndef VOXELISER_HPP
#define VOXELISER_HPP

#include <iostream>
#include <filesystem>
#include <vector>
#include <future>
#include <exception>

#include <hdf5.h>
#include <omp.h>

#include "turbodorn.hpp"
#include "hdf5.hpp"

namespace TurboDorn {

    namespace Detail {

    }//Detail


    class Voxeliser {


        Eigen::VectorXcd chunk(Eigen::VectorXi _x, Eigen::VectorXi _e) {

        }


    public:

        void compute(size_t _f_in, size_t _voxels, const std::filesystem::path& _hdf5_file_out) {

            //parse hdf5 file
            IO::Detail::File file
            switch(_f_in){
                case 0:
                {
                    throw std::runtime_error("file does not exist");
                }
                case 1:
                {
                    file = IO::simulation_results();
                    break;
                }
                case 2:
                {
                    file = IO::simulation_results_phi000();
                    break;
                }
                case 3:
                {
                    file = IO::simulation_results_phi100();
                    break;
                }
                case 4:
                {
                    file = IO::simulation_results_phi121();
                    break;
                }
                case 5:
                {
                    file = IO::simulation_results_phi412();
                    break;
                }
                default: throw std::runtime_error("file does not exist");
            }

            //precompute invariants
            const Detail::Invariants inv = Detail::computeInvariants(file);

            //compute chunks and write them to hdf5 file
            for(size_t x = 0; x < _voxels; ++x){
                for(size_t x = 0; x < _voxels; ++x){
                    for(size_t x = 0; x < _voxels; ++x){

                        
                        const Eigen::VectorXcd&& res = chunk();
                
                    }
                }
            }

        }

    };//Voxeliser

}//TurboDorn

#endif //VOXELISER_HPP