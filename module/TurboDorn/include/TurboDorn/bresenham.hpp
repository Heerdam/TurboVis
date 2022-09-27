#ifndef BRESENHAM_HPP
#define BRESENHAM_HPP

#include <cassert>

#include <Eigen/Eigen>
#include <tsl/robin_map.h>

namespace TurboDorn{
    namespace Detail {

    } //Detail

    class Bresenham {

    public:
        Bresenham(const Eigen::VectorXd& _p1, const Eigen::VectorXd& _p2, const double _m) {
            assert(_p1.rows() == _p2.rows());
            
            const size_t dims = _p1.rows();

            const Eigen::Matrix<T, -1, 1> dp = _m * (_p2 - _p1);

            Eigen::Matrix<I, -1, 1> in (dp.rows());

            Eigen::Matrix<I, -1, 1> eps (dp.rows());
            



        }//Bresenham

        void Bresenham3D(int x1, int y1, int z1, const int x2, const int y2, const int z2, WorldMap *output, int symbol){
    
            int i, dx, dy, dz, l, m, n, x_inc, y_inc, z_inc, err_1, err_2, dx2, dy2, dz2;
            int point[3];
            
            point[0] = x1;
            point[1] = y1;
            point[2] = z1;

            dx = x2 - x1;
            dy = y2 - y1;
            dz = z2 - z1;
            
            x_inc = (dx < 0) ? -1 : 1;
            l = abs(dx);
            y_inc = (dy < 0) ? -1 : 1;
            m = abs(dy);
            z_inc = (dz < 0) ? -1 : 1;
            n = abs(dz);
            dx2 = l << 1;
            dy2 = m << 1;
            dz2 = n << 1;
            
            if ((l >= m) && (l >= n)) {
                err_1 = dy2 - l;
                err_2 = dz2 - l;
                for (i = 0; i < l; i++) {
                    output->getTileAt(point[0], point[1], point[2])->setSymbol(symbol);
                    if (err_1 > 0) {
                        point[1] += y_inc;
                        err_1 -= dx2;
                    }
                    if (err_2 > 0) {
                        point[2] += z_inc;
                        err_2 -= dx2;
                    }
                    err_1 += dy2;
                    err_2 += dz2;
                    point[0] += x_inc;
                }
            } else if ((m >= l) && (m >= n)) {
                err_1 = dx2 - m;
                err_2 = dz2 - m;
                for (i = 0; i < m; i++) {
                    output->getTileAt(point[0], point[1], point[2])->setSymbol(symbol);
                    if (err_1 > 0) {
                        point[0] += x_inc;
                        err_1 -= dy2;
                    }
                    if (err_2 > 0) {
                        point[2] += z_inc;
                        err_2 -= dy2;
                    }
                    err_1 += dx2;
                    err_2 += dz2;
                    point[1] += y_inc;
                }
            } else {
                err_1 = dy2 - n;
                err_2 = dx2 - n;
                for (i = 0; i < n; i++) {
                    output->getTileAt(point[0], point[1], point[2])->setSymbol(symbol);
                    if (err_1 > 0) {
                        point[1] += y_inc;
                        err_1 -= dz2;
                    }
                    if (err_2 > 0) {
                        point[0] += x_inc;
                        err_2 -= dz2;
                    }
                    err_1 += dy2;
                    err_2 += dx2;
                    point[2] += z_inc;
                }
            }
            output->getTileAt(point[0], point[1], point[2])->setSymbol(symbol);
        }

    };
}

#endif //BRESENHAM_HPP
