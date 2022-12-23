#ifndef TURBODORN_HPP
#define TURBODORN_HPP

#include <complex>
#include <numeric>
#include <vector>
#include <filesystem>
#include <fstream>
#include <cassert>

#include <Eigen/Eigen>
#include <tsl/robin_map.h>

#include <highfive/H5File.hpp>
#include <highfive/h5easy_bits/H5Easy_Eigen.hpp>

#include <nlohmann/json.hpp>

using namespace nlohmann;
using namespace std::chrono_literals;

namespace TurboDorn {

    namespace Griderator {

        using Vector3 = Eigen::Vector3d;
        using Point3 = Eigen::Vector3d;
        using Vector3int32 = Eigen::Vector3i;
        using Point3int32 = Eigen::Vector3i;

        using int32 = int32_t;

        namespace Detail {

            inline Vector3 operator/ (const Vector3& lkVector, const Vector3& rkVector) {
                return Vector3(lkVector.x() / rkVector.x(), lkVector.y() / rkVector.y(), lkVector.z() / rkVector.z());
            }

            template <typename T> 
            int sign(T _val) {
                return (T(0) < _val) - (_val < T(0));
            }

            template<typename T>
            T square(T _v) {
                return _v * _v;
            }

            class AABox {
                /** NaN if empty */
                Point3 lo;

                /** NaN if empty */
                Point3 hi;

            public:
                /** Creates the empty bounds, i.e., an empty set of points. */
                AABox() : lo(NAN, NAN, NAN), hi(NAN, NAN, NAN) {}

                AABox(Point3 _lo, Point3 _hi) : lo(std::move(_lo)), hi(std::move(_hi)) {}

                /**
                 Constructs a zero-volume AABox at v.
                */
                explicit AABox(const Point3& v) {
                    lo = hi = v;
                }

                /** Returns not-a-number if empty */
                const Point3& low() const {
                    return lo;
                }

                /** Returns not-a-number if empty */
                const Point3& high() const {
                    return hi;
                }

                /** Returns the centroid of the box (NaN if empty) */
                Point3 center() const {
                    return (lo + hi) * 0.5;
                }

                Point3 corner(int index) const {
                    // default constructor inits all components to 0
                    Vector3 v;

                    switch (index) {
                    case 0:
                        v.x() = lo.x();
                        v.y() = lo.y();
                        v.z() = hi.z();
                        break;

                    case 1:
                        v.x() = hi.x();
                        v.y() = lo.y();
                        v.z() = hi.z();
                        break;

                    case 2:
                        v.x() = hi.x();
                        v.y() = hi.y();
                        v.z() = hi.z();
                        break;

                    case 3:
                        v.x() = lo.x();
                        v.y() = hi.y();
                        v.z() = hi.z();
                        break;

                    case 4:
                        v.x() = lo.x();
                        v.y() = lo.y();
                        v.z() = lo.z();
                        break;

                    case 5:
                        v.x() = hi.x();
                        v.y() = lo.y();
                        v.z() = lo.z();
                        break;

                    case 6:
                        v.x() = hi.x();
                        v.y() = hi.y();
                        v.z() = lo.z();
                        break;

                    case 7:
                        v.x() = lo.x();
                        v.y() = hi.y();
                        v.z() = lo.z();
                        break;

                    default:
                        assert(false && "Invalid corner index");
                        break;
                    }
                    return v;
                }

                bool isEmpty() const {
                    if(std::isnan(lo[0])) return true;
                    if(std::isnan(lo[1])) return true;
                    if(std::isnan(lo[2])) return true;
                }

                /** Distance from corner(0) to the next corner along axis a. */
                float extent(int a) const {
                    if (isEmpty()) return 0.0f;
                    assert(a < 3);
                    return hi[a] - lo[a];
                }

                Vector3 extent() const {
                    if (isEmpty()) return Vector3::Zero();
                    return hi - lo;
                }

            };

            class Ray {
                Point3 m_origin;

                double m_minDistance;

                /** Unit length */
                Vector3 m_direction;

                double m_maxDistance;

            public:
                const Point3& origin() const {
                    return m_origin;
                }
                
                /** Unit direction vector. */
                const Vector3& direction() const {
                    return m_direction;
                }

                /** \param direction Assumed to have unit length */
                void set(const Point3& origin, const Vector3& direction, double minDistance = 0.0f, double maxDistance = INFINITY);

                Ray() {
                    set(Point3::Zero(), Vector3::UnitX());
                }

                /** \param direction Assumed to have unit length */
                Ray(const Point3& origin, const Vector3& direction, double minDistance = 0.0f, double maxDistance = INFINITY) {
                    set(origin, direction, minDistance, maxDistance);
                }

                /**
                 Creates a Ray from a origin and a (nonzero) unit direction.
                */
                static Ray fromOriginAndDirection(const Point3& point, const Vector3& direction, double minDistance = 0.0f, double maxDistance = INFINITY) {
                    return Ray(point, direction, minDistance, maxDistance);
                }

            };

            /**
                @brief Calculates intersection of a ray and a static Axis-Aligned Box (AABox).

                @note Avoids the sqrt from collisionTimeForMovingPointFixedAABox; 
                early-out branches and operations optimized for Intel Core2 architecture.
                
                @param invDir      1/dir
                @param location    Location of collision. [Post Condition]
                @param inside      Does the ray originate inside the box? [Post Condition]

                @return True if the ray hits the box
            */
            bool inter_rayAABox(
                const Ray& ray, 
                const Vector3& invDir, 
                const AABox& box, 
                const Vector3& boxCenter,
                double boundingRadiusSquared,
                Vector3& location,
                bool& inside) {
                
                //assert(fabs(ray.direction().squaredLength() - 1.0f) < 0.01f, format("Length = %f", ray.direction().length()));
                {
                    // Pre-emptive partial bounding sphere test
                    const Vector3 L(boxCenter - ray.origin());
                    double d = L.dot(ray.direction());

                    double L2 = L.dot(L);
                    double D2 = d * d;
                    double M2 = L2 - D2;

                    if (((d < 0) && (L2 > boundingRadiusSquared)) || (M2 > boundingRadiusSquared)) {
                        inside = false;
                        return false;
                    }
                    // Passing here does not mean that the ray hits the bounding sphere;
                    // we would still have to perform more expensive tests to determine
                    // that.
                }

                inside = true;
                const Vector3& MinB = box.low();
                const Vector3& MaxB = box.high();
                Vector3 MaxT(-1.0f, -1.0f, -1.0f);

                // Find candidate planes.
                for (int i = 0; i < 3; ++i) {
                    if (ray.origin()[i] < MinB[i]) {
                        location[i]    = MinB[i];
                        inside      = false;
                        
                        // Calculate T distances to candidate planes
                        if (ray.direction()[i] != 0.) {
                            MaxT[i] = (MinB[i] - ray.origin()[i]) * invDir[i];
                        }
                    } else if (ray.origin()[i] > MaxB[i]) {
                        location[i]    = MaxB[i];
                        inside      = false;

                        // Calculate T distances to candidate planes
                        if (ray.direction()[i] != 0.) {
                            MaxT[i] = (MaxB[i] - ray.origin()[i]) * invDir[i];
                        }
                    }
                }

                if (inside) {
                    // Ray origin inside bounding box
                    location = ray.origin();
                    return true;
                }
                
                // Get largest of the maxT's for final choice of intersection
                int WhichPlane = 0;
                if (MaxT[1] > MaxT[WhichPlane]) {
                    WhichPlane = 1;
                }

                if (MaxT[2] > MaxT[WhichPlane]) {
                    WhichPlane = 2;
                }

                // Check final candidate actually inside box
                if (MaxT[WhichPlane] < 0.0f) {
                    // Miss the box
                    return false;
                }

                for (int i = 0; i < 3; ++i) {
                    if (i != WhichPlane) {
                        location[i] = ray.origin()[i] + MaxT[WhichPlane] * ray.direction()[i];
                        if ((location[i] < MinB[i]) ||
                            (location[i] > MaxB[i])) {
                            // On this plane we're outside the box extents, so
                            // we miss the box
                            return false;
                        }
                    }
                }
                
                return true;
            }

        }  // namespace Detail

        using namespace Detail;

        /**
        Computes conservative line raster/voxelization across a grid for
        use in walking a grid spatial data structure or or voxel scene
        searching for intersections.  At each iteration, the iterator
        steps exactly one cell in exactly one dimension.

        Example of this iterator applied to ray-primitive intersection in a
        grid:

        \code
        bool firstRayIntersection(const Ray& ray, Value*& value, float& distance) const {

        for (RayGridIterator it(ray, cellSize); inBounds(it.index); ++it) {
            // Search for an intersection within this grid cell
            const Cell& c = cell(it.index);
            float maxdistance = min(distance, t.tExit);
            if (c.firstRayIntersection(ray, value, maxdistance)) {
                distance = maxdistance;
                return true;
            }
        }
        }
        \endcode

        \sa CollisionDetection, PointHashGrid, Ray, Intersect

        */
        class RayGridIterator {
            
        protected:
            /** Extent of the grid in each dimension, in grid cell units.*/
            Vector3int32 m_numCells;

            /** Current grid cell m_index */
            Vector3int32 m_index;

            /** Sign of the direction that the ray moves along each axis; +/-1
             or 0 */
            Vector3int32 m_step;

            /** Size of one cell in units of t along each axis. */
            Vector3 m_tDelta;

            /** Distance along the ray of the first intersection with the
             current cell (i.e., that given by m_index). Zero for the cell that contains the ray origin. */
            double m_enterDistance;

            /** Distance along the ray to the intersection with the next grid
             cell.  enterDistance and m_exitDistance can be used to bracket ray
            ray-primitive intersection tests within a cell.
            */
            Vector3 m_exitDistance;

            /** The axis along which the ray entered the cell; X = 0, Y = 1, Z = 2.  This is always set to X for the cell that contains the ray origin. */
            int m_enterAxis;

            /** The original ray */
            Ray m_ray;

            /** Size of each cell along each axis */
            Vector3 m_cellSize;

            /** True if index() refers to a valid cell inside the grid.  This
                is usually employed as the loop termination condition.*/
            bool m_insideGrid;

            /** The value that the index will take on along each boundary when
                it just leaves the grid. */
            Vector3int32 m_boundaryIndex;

            /** True if this cell contains the ray origin */
            bool m_containsRayOrigin;

        public:
            /** \copydoc m_ray */
            const Ray& ray() const {
                return m_ray;
            }

            /** \copydoc m_numCells */
            Vector3int32 numCells() const {
                return m_numCells;
            }

            /** \copydoc m_enterAxis */
            int enterAxis() const {
                return m_enterAxis;
            }

            /** Outward-facing normal to the current grid cell along the
             partition just entered.  Initially zero. */
            Vector3int32 enterNormal() const {
                Vector3int32 normal(0, 0, 0);
                normal[m_enterAxis] = -m_step[m_enterAxis];
                return normal;
            }

            /** \copydoc m_cellSize */
            const Vector3& cellSize() const {
                return m_cellSize;
            }

            /** Location where the ray entered the current grid cell */
            Point3 enterPoint() const {
                return m_ray.origin() + m_enterDistance * m_ray.direction();
            }

            /** Location where the ray exits the current grid cell */
            Point3 exitPoint() const {
                return m_ray.origin() + m_exitDistance.minCoeff() * m_ray.direction();
            }

            /** \copydoc m_enterDistance */
            double enterDistance() const {
                return m_enterDistance;
            }

            /** Distance from the ray origin to the exit point in this cell. */
            double exitDistance() const {
                return m_exitDistance.minCoeff();
            }

            /** \copydoc m_step */
            const Vector3int32& step() const {
                return m_step;
            }

            /** \copydoc m_index */
            const Vector3int32& index() const {
                return m_index;
            }

            /** \copydoc m_tDelta */
            const Vector3& tDelta() const {
                return m_tDelta;
            }

            /** \copydoc m_insideGrid */
            bool insideGrid() const {
                return m_insideGrid;
            }

            /** \copydoc m_containsRayOrigin */
            bool containsRayOrigin() const {
                return m_containsRayOrigin;
            }

            /**
                \brief Initialize the iterator to the first grid cell hit by
                the ray and precompute traversal variables.

                The grid is assumed to have a corner at (0,0,0) and extend
                along the canonical axes.  For intersections grids transformed
                by a rigid body transformation, first transform the ray into
                the grid's object space with CFrame::rayToObjectSpace.

                If the ray never passes through the grid, insideGrid() will
                be false immediately after intialization.

                If using for 2D iteration, set <code>numCells.z = 1</code> and
                <code>ray.origin().z = 0.5</code>

                \param cellSize The extent of one cell

                \param minBoundsLocation The location of the lowest corner of grid cell minBoundsCellIndex along each axis.
                This translates the grid relative to the ray's coordinate frame.

                \param minBoundsCellIndex The index of the first grid cell.  This allows
                operation with grids defined on negative indices.  This translates all
                grid indices.
            */
            RayGridIterator(Ray ray,
                                        const Vector3int32& numCells,
                                        const Vector3& cellSize,
                                        const Point3& gridOrigin,
                                        const Point3int32& gridOriginIndex) : m_numCells(numCells),
                                                                            m_enterDistance(0.0f),
                                                                            m_enterAxis(0),
                                                                            m_ray(ray),
                                                                            m_cellSize(cellSize),
                                                                            m_insideGrid(true),
                                                                            m_containsRayOrigin(true) {
                if (!gridOrigin.isZero()) {
                    // Change to the grid's reference frame
                    ray = Ray::fromOriginAndDirection(ray.origin() - gridOrigin, ray.direction());
                }

                //////////////////////////////////////////////////////////////////////
                // See if the ray begins inside the box
                const AABox gridBounds(Vector3::Zero(), Vector3(numCells(0) * cellSize(0), numCells(1) * cellSize(1), numCells(2) * cellSize(2)));

                bool startsOutside = false;
                bool inside = false;
                Point3 startLocation = ray.origin();

                const bool passesThroughGrid = inter_rayAABox(
                    ray, 
                    Vector3::Ones() / ray.direction(),
                    gridBounds, 
                    gridBounds.center(),
                    square(gridBounds.extent().norm() * 0.5f),
                    startLocation,
                    inside
                );

                if (!inside) {
                    if (passesThroughGrid) {
                        // Back up slightly so that we immediately hit the
                        // start location.  The precision here is tricky--if
                        // the ray strikes at a very glancing angle, we need
                        // to move a large distance along the ray to enter the
                        // grid.  If the ray strikes head on, we only need to
                        // move a small distance.
                        m_enterDistance = (ray.origin() - startLocation).norm() - 0.0001f;
                        startLocation = ray.origin() + ray.direction() * m_enterDistance;
                        startsOutside = true;
                    } else {
                        // The ray never hits the grid
                        m_insideGrid = false;
                    }
                }

                //////////////////////////////////////////////////////////////////////
                // Find the per-iteration variables
                for (int a = 0; a < 3; ++a) {
                    m_index[a] = (int32)floor(startLocation[a] / cellSize[a]);
                    m_tDelta[a] = cellSize[a] / fabs(ray.direction()[a]);

                    m_step[a] = (int32)sign(ray.direction()[a]);

                    // Distance to the edge fo the cell along the ray direction
                    float d = startLocation[a] - m_index[a] * cellSize[a];
                    if (m_step[a] > 0) {
                        // Measure from the other edge
                        d = cellSize[a] - d;

                        // Exit on the high side
                        m_boundaryIndex[a] = m_numCells[a];
                    } else {
                        // Exit on the low side (or never)
                        m_boundaryIndex[a] = -1;
                    }
                    assert(d >= 0 && d <= cellSize[a]);

                    if (ray.direction()[a] != 0) {
                        m_exitDistance[a] = d / fabs(ray.direction()[a]) + m_enterDistance;
                    } else {
                        // Ray is parallel to this partition axis.
                        // Avoid dividing by zero, which could be NaN if d == 0
                        m_exitDistance[a] = (float)INFINITY;
                    }
                }

                if (!gridOriginIndex.isZero()) {
                    // Offset the grid coordinates
                    m_boundaryIndex += gridOriginIndex;
                    m_index += gridOriginIndex;
                }

                if (startsOutside) {
                    // Let the increment operator bring us into the first cell
                    // so that the starting axis is initialized correctly.
                    ++(*this);
                }
            }

            /** Increment the iterator, moving to the next grid cell */
            RayGridIterator& operator++() {
                // Find the axis of the closest partition along the ray
                if (m_exitDistance.x() < m_exitDistance.y()) {
                    if (m_exitDistance.x() < m_exitDistance.z()) {
                        m_enterAxis = 0;
                    } else {
                        m_enterAxis = 2;
                    }
                } else if (m_exitDistance.y() < m_exitDistance.z()) {
                    m_enterAxis = 1;
                } else {
                    m_enterAxis = 2;
                }

                m_enterDistance = m_exitDistance[m_enterAxis];
                m_index[m_enterAxis] += m_step[m_enterAxis];
                m_exitDistance[m_enterAxis] += m_tDelta[m_enterAxis];

                // If the index just hit the boundary exit, we have
                // permanently exited the grid.
                m_insideGrid = m_insideGrid &&
                            (m_index[m_enterAxis] != m_boundaryIndex[m_enterAxis]);

                m_containsRayOrigin = false;

                return *this;
            }
        };

    }//Griderator

    namespace Hagedorn::Detail {
        std::vector<Eigen::VectorXi> hyperbolic_cut_shape(size_t _dim, size_t _K) noexcept;
        Eigen::Index index(const Eigen::VectorXi& _i, const Eigen::VectorXi& _e ) noexcept;
    }

    namespace IO {

        namespace Detail {

            //---------------------------------------------------------------------------------------//
            //                                        FilePathResolver
            //---------------------------------------------------------------------------------------//
            class FilePathResolver {

                std::filesystem::path path_asset;
            
            public:

                FilePathResolver() {
                    //asset
                    {
                        const auto s_path_1 = std::filesystem::current_path() /= "assets/";
                        const auto s_path_2 = std::filesystem::current_path() /= "../assets/";
                        const auto s_path_3 = std::filesystem::current_path() /= "../../assets/";

                        if(std::filesystem::exists(s_path_1))
                            path_asset = s_path_1;
                        else if(std::filesystem::exists(s_path_2))
                            path_asset = s_path_2;
                        else if(std::filesystem::exists(s_path_3))
                            path_asset = s_path_3;
                        else throw std::runtime_error("folder asset not found");
                    }
                }

                std::filesystem::path operator()() const {
                    return path_asset;
                }

            };//FilePathResolver

            //---------------------------------------------------------------------------------------//
            //                                        File
            //---------------------------------------------------------------------------------------//

            template<class T>
            struct File {
                size_t dimensions;
                size_t timesteps;
                size_t K;
                T epsilon = 1.;  
                Eigen::VectorXcd S;       
                std::vector<Eigen::VectorX<std::complex<T>>> c_0;
                std::vector<Eigen::VectorX<std::complex<T>>> p, q;
                std::vector<Eigen::MatrixX<std::complex<T>>> P, Q;

                //shapefunction
                std::vector<Eigen::VectorXi> Ks; 
                Eigen::VectorXi k_max; //max k in d direction
                tsl::robin_map<Eigen::Index, bool> b_Ks;
            };//File

            //---------------------------------------------------------------------------------------//
            //                                   load_from_file
            //---------------------------------------------------------------------------------------//

            template<class T, class I, class CUTSHAPE>
            File<T> load_from_file(const std::filesystem::path& _path, I _dims, I _K, CUTSHAPE&& _cs) {

                using Vector = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>;
                using Matrix = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;

                if(!std::filesystem::exists(_path))
                    throw std::runtime_error("file does not exist");

                File<T> out;
                out.dimensions = _dims;
                out.K = _K;

                try {

                    //file
                    HighFive::File file(_path.string(), HighFive::File::ReadOnly);

                    //timesteps
                    {
                        const Eigen::Matrix<T, Eigen::Dynamic, 1> steps = H5Easy::load<Eigen::Matrix<T, Eigen::Dynamic, 1>>(file, "datablock_0/wavepacket/timegrid");
                        out.timesteps = steps.rows();
                    }

                    //S
                    {
                        const H5Easy::DataSet dataset = file.getDataSet("datablock_0/wavepacket/Pi/S");
                        std::vector<std::vector<std::vector<std::complex<T>>>> S;
                        dataset.read(S);
                        out.S.resize(S.size());
                        for(size_t i = 0; i < S.size(); ++i){
                            out.S(i) = S[i][0][0];
                        }
                    }
                
                    //c_0
                    {
                        const Matrix tmp = H5Easy::load<Matrix>(file, "datablock_0/wavepacket/coefficients/c_0");
                        for(size_t t = 0; t < out.timesteps; ++t){
                            Vector c_0(tmp.cols());
                            for(size_t d = 0; d < tmp.cols(); ++d)
                                c_0(d) = tmp(t, d);
                            out.c_0.push_back(std::move(c_0));
                        }
                    }
                
                    //p
                    {
                        const H5Easy::DataSet dataset = file.getDataSet("datablock_0/wavepacket/Pi/p");
                        std::vector<std::vector<std::vector<std::complex<T>>>> p;
                        dataset.read(p);
                        for(size_t t = 0; t < p.size(); ++t){
                            Vector temp(out.dimensions);
                            for(size_t y = 0; y < p[t].size(); ++y){
                                temp(y) = p[t][y][0];
                            }
                            out.p.push_back(std::move(temp));
                        }
                    }
                
                    //q
                    {
                        const H5Easy::DataSet dataset = file.getDataSet("datablock_0/wavepacket/Pi/q");
                        std::vector<std::vector<std::vector<std::complex<T>>>> q;
                        dataset.read(q);
                        for(size_t t = 0; t < q.size(); ++t){
                            Vector temp(out.dimensions);
                            for(size_t y = 0; y < q[t].size(); ++y){
                                temp(y) = q[t][y][0];
                            }
                            out.q.push_back(std::move(temp));
                        }
                    }
            
                    //P
                    {
                        const H5Easy::DataSet dataset = file.getDataSet("datablock_0/wavepacket/Pi/P");
                        std::vector<std::vector<std::vector<std::complex<T>>>> P;
                        dataset.read(P);
                        for(size_t t = 0; t < P.size(); ++t){
                            Matrix temp(out.dimensions, out.dimensions);
                            for(size_t x = 0; x < out.dimensions; ++x){
                                for(size_t y = 0; y < out.dimensions; ++y)
                                    temp(x, y) = P[t][x][y];
                            }
                            out.P.push_back(std::move(temp));
                        }
                    }

                    //Q
                    {
                        const H5Easy::DataSet dataset = file.getDataSet("datablock_0/wavepacket/Pi/Q");
                        std::vector<std::vector<std::vector<std::complex<T>>>> Q;
                        dataset.read(Q);
                        for(size_t t = 0; t < Q.size(); ++t){
                            Matrix temp(out.dimensions, out.dimensions);
                            for(size_t x = 0; x < out.dimensions; ++x){
                                for(size_t y = 0; y < out.dimensions; ++y)
                                    temp(x, y) = Q[t][x][y];
                            }
                            //std::cout << temp << std::endl;
                            out.Q.push_back(std::move(temp));
                        }
                    }
                    

                } catch(const std::exception& _e){
                    throw _e;
                }

                // -------------- cut shape --------------
                {
                    out.Ks = _cs(out.dimensions, out.K);

                    out.k_max = Eigen::VectorXi(out.dimensions);
                    out.k_max.setZero();

                    for(size_t i = 0; i < out.dimensions; ++i){
                        for(const auto& k : out.Ks){
                            out.k_max(i) = std::max(out.k_max(i), k(i));
                        }
                    }
            
                    for(const auto& c : out.Ks){
                        const Eigen::Index ii = ::TurboDorn::Hagedorn::Detail::index(c, out.k_max);
                        //std::cout << c(0) << ", " << c(1) << ", " << c(2) << " | " << out.k_max(0) << ", " << out.k_max(1) << ", " << out.k_max(2) << " | " << ii << std::endl;
                        out.b_Ks.insert( {ii, true} );
                    }
                }

                return out ;
            }; //loadFromFile

        }

        Detail::File<double> simulation_results() {
            try {
                const auto path = Detail::FilePathResolver()() /= "simulation_results.hdf5";
                return Detail::load_from_file<double> (path, 3, 4, ::TurboDorn::Hagedorn::Detail::hyperbolic_cut_shape);
            } catch(const std::exception& _e){
                throw _e;
            }
        }//simulation_results

        Detail::File<double> simulation_results_phi000() {
            try {
                const auto path = Detail::FilePathResolver()() /= "simulation_results_phi000.hdf5";
                return Detail::load_from_file<double>(path, 3, 1, ::TurboDorn::Hagedorn::Detail::hyperbolic_cut_shape);
            } catch(const std::exception& _e){
                throw _e;
            }
        }//simulation_results_phi000

        Detail::File<double> simulation_results_phi100() {
            try {
                const auto path = Detail::FilePathResolver()() /= "simulation_results_phi100.hdf5";
                return Detail::load_from_file<double>(path, 3, 2, ::TurboDorn::Hagedorn::Detail::hyperbolic_cut_shape);
            } catch(const std::exception& _e){
                throw _e;
            }
        }//simulation_results_phi100

        Detail::File<double> simulation_results_phi121() {
            try {
                const auto path = Detail::FilePathResolver()() /= "simulation_results_phi121.hdf5";
                return Detail::load_from_file<double>(path, 3, 12, ::TurboDorn::Hagedorn::Detail::hyperbolic_cut_shape);
            } catch(const std::exception& _e){
                throw _e;
            }
        }//simulation_results_phi121

        Detail::File<double> simulation_results_phi412() {
            try {
                const auto path = Detail::FilePathResolver()() /= "simulation_results_phi412.hdf5";
                return Detail::load_from_file<double>(path, 3, 30, ::TurboDorn::Hagedorn::Detail::hyperbolic_cut_shape);
            } catch(const std::exception& _e){
                throw _e;
            }
        }//simulation_results_phi412

    }//IO

    namespace Hagedorn {

        namespace Detail {

            //---------------------------------------------------------------------------------------//
            //                                        Invariants
            //---------------------------------------------------------------------------------------//

            template<class T>
            struct Invariants {
                size_t dimensions;
                Eigen::VectorXi k; //k extends
                tsl::robin_map<Eigen::Index, bool> k_shape; //lookup for shape
                std::vector<Eigen::VectorX<std::complex<T>>> p;
                std::vector<Eigen::VectorX<std::complex<T>>> q;
                //phi_0
                std::vector<std::complex<T>> pre;
                std::complex<T> i_2_E_2;
                std::vector<Eigen::MatrixX<std::complex<T>>> P_Q_1;
                std::vector<Eigen::RowVectorX<std::complex<T>>> i_E_2_p;
                //phi
                //sqrt*Q-1
                std::vector<Eigen::MatrixX<std::complex<T>>> Q_1;
                //Q-1*QT
                std::vector<Eigen::MatrixX<std::complex<T>>> Q_1_Q_T;
                //
                std::vector<Eigen::VectorX<std::complex<T>>> c_0;
                std::vector<Eigen::VectorXi> Ks; 
                Eigen::VectorXi k_max; //max k in d direction
                Eigen::VectorXcd S;  
            };//Invariants

            //---------------------------------------------------------------------------------------//
            //                             prepare_invariants_from_file
            //---------------------------------------------------------------------------------------//

            template<class T, template<typename> class FILE>
            EIGEN_STRONG_INLINE std::unique_ptr<Invariants<T>> prepare_invariants_from_file(const FILE<T>& _file) {
                std::unique_ptr<Invariants<T>> out = std::make_unique<Invariants<T>>();
                out->dimensions = _file.dimensions;
                out->k = _file.k_max;
                out->i_2_E_2 = std::complex<T>(0., 1.) / (2. * _file.epsilon * _file.epsilon);
                out->p = _file.p;
                out->q = _file.q;
                out->k_shape = _file.b_Ks;

                for(size_t t = 0; t < _file.timesteps; ++t){
                    //phi0
                    out->pre.push_back( std::pow(T(M_PI) * _file.epsilon * _file.epsilon, - T(_file.dimensions) / 4.) * std::pow(_file.Q[t].determinant(), -0.5) );
                    out->P_Q_1.push_back( _file.P[t] * _file.Q[t].inverse() );

                    out->i_E_2_p.push_back( (std::complex<T>(0., 1.) / _file.epsilon * _file.epsilon) * _file.p[t].transpose() );
                    //phi
                    out->Q_1.push_back( std::sqrt(2. / (_file.epsilon * _file.epsilon)) * _file.Q[t].inverse() );
                    out->Q_1_Q_T.push_back( _file.Q[t].inverse() * _file.Q[t].conjugate() );
                }
                
                // -------------- copy some things --------------
                {
                    out->c_0 = _file.c_0;
                    out->Ks = _file.Ks;
                    out->k_max = _file.k_max;
                    out->S = _file.S;
                }
                return out;
            }

            //---------------------------------------------------------------------------------------//
            //                                            index
            //---------------------------------------------------------------------------------------//

            /*
            i: index
            e: extents, # units
            */
            EIGEN_STRONG_INLINE Eigen::Index index(const Eigen::VectorXi& _i, const Eigen::VectorXi& _e ) noexcept {
                assert(_i.size() == _e.size());
                Eigen::Index out = _i(0);
                for (Eigen::Index k = 1; k < _i.size(); ++k) {
                    out *= _e(k);
                    out += _i(k);
                }
                return out;
            }; //index

            //---------------------------------------------------------------------------------------//
            //                             hyperbolic_cut_shape
            //---------------------------------------------------------------------------------------//

            EIGEN_STRONG_INLINE std::vector<Eigen::VectorXi> hyperbolic_cut_shape(size_t _dim, size_t _K) noexcept {
                std::vector<Eigen::VectorXi> out;

                Eigen::VectorXi index(_dim);
                index.setZero();

                while (true) {
                    for (index(_dim - 1) = 0; index(_dim - 1) <= _K; ++index(_dim - 1)) {
                        Eigen::Index p = 1;
                        for (Eigen::Index d = 0; d < _dim; ++d)
                            p *= (1 + index(d));
                        if (p <= _K)
                            out.push_back(index);
                    }

                    bool done = false;
                    for (Eigen::Index d = _dim - 2; d >= 0; --d) {
                        index(d) += 1;

                        if (index(d) >= _K) {
                            if (d == 0)
                                done = true;
                            else
                                index(d) = 0;
                        } else
                            break;
                    }
                    if (done) break;
                }

                return out;
            }; //hyperbolic_cut_shape

            //---------------------------------------------------------------------------------------//
            //                                            
            //---------------------------------------------------------------------------------------//

            //struct StreamingOctreeOptions {

            //};

            /*
                voxel size: size of a voxel in world coords
                depth: how many subdivides the octree has
                chunk size: symmetric 2^n extends -> size: 4*2*(2^n)^3 bytes (n = 32 (262'144 bytes) or 64 (2'097'152 bytes))

                only lowest level has data. rest only boolmaps

                format:
                double                              voxel size
                size_t                              depth
                size_t                              chunk size
                n bits with 1 byte/subdivision      ((8^(n+1))-1)/7-1 nodes
                chunks...                           chunk_size^3 bytes for every none-empty chunk

            */
        /*
            template<class T = double>
            class Chunkerator {

                iVec extends;

                bVec bmap;
                cVec temp_chunk;

                //iterator state
                const bool isEndIterator = false;

                //file meta information
                std::filesystem::path path;
                dVec camPos, camPlaneNorm, halfExtCamPlane;


            public:

                using iterator_category = std::forward_iterator_tag;

                /*
                    path - path to the file
                    p - position of the camera plane in world coords
                    n - normalized plane normal of the camera plane
                    e - half extends of the camera plane in world units
                */
                //Chunkerator(std::filesystem::path&& _path, dVec&& _p, dVec&& _n, dVec&& _e) : 
                    //   path(std::move(_path)), camPos(std::move(_p)),  camPlaneNorm(std::move(_n)), halfExtCamPlane(std::move(_e)){}

                //creates an end chunkerator
                //Chunkerator() : isEndIterator(true) {};

                //void intitialise() {
                //    if(!isEndIterator) throw std::runtime_error("Can't initialize an end chunkerator!");
                    
                //}

                //void operator*() {

                //}

                //void operator++() {

                //}

                //void operator++(int) {

                //}

                //bool operator==(const Chunkerator<T>& _other) const noexcept {
                //    if(isEndIterator && _other.isEndIterator) return true;
                //}

            //};//Chunkerator

            //---------------------------------------------------------------------------------------//
            //                                    ray_aabb_intersect
            //---------------------------------------------------------------------------------------//

            template<class T>
            EIGEN_STRONG_INLINE bool ray_aabb_intersect(const Eigen::Vector3<T>& _r_o, const Eigen::Vector3<T>& _r_d, const Eigen::Vector3<T>& _low, const Eigen::Vector3<T>& _high, T _tmax, T& _t) noexcept {
                _t = -std::numeric_limits<T>::infinity();
                for (size_t i = 0; i < _r_o.size(); ++i) {
                    if (std::abs(_r_d[i]) < std::numeric_limits<T>::epsilon()){
                        if (_r_o[i] < _low[i] || _r_o[i] > _high[i]) 
                            return false;
                    } else {
                        const T ood = 1. / _r_d[i];
                        T t1 = (_low[i] - _r_o[i]) * ood;
                        T t2 = (_high[i] - _r_o[i]) * ood;
                        if (t1 > t2) std::swap(t1, t2);
                        _t = std::max(_t, t1);
                        _tmax = std::min(_tmax, t2);
                        if (_t > _tmax) return false;
                    }
                }
                return true;
            }//intersect

            //---------------------------------------------------------------------------------------//
            //                                        c_to_HSL
            //---------------------------------------------------------------------------------------//

            template <class T, T MAX>
            EIGEN_STRONG_INLINE constexpr Eigen::Vector3<T> c_to_HSL(const std::complex<T>& _c) noexcept {
                const T H = std::clamp(std::abs(std::fmod(std::arg(_c), 2. * M_PI)), 0., 2. * M_PI);
                const T S = 1.;
                const T L = std::clamp(std::abs(MAX * std::atan(std::abs(_c)) / (0.5 * M_PI)), 0., 1.);
                return { H, S, L };
            }//c_to_HSL

            //---------------------------------------------------------------------------------------//
            //                                        HSL_to_RGB_deg
            //---------------------------------------------------------------------------------------//

            /*.
                h: [0, 360]
                s: [0, 1]
                l: [0, 1]
                rgb: [0, 1]
            */
            template <class T>
            EIGEN_STRONG_INLINE constexpr Eigen::Vector3<T> HSL_to_RGB_deg(const Eigen::Matrix<T, 3, 1>& _hsl) noexcept {

                const T H = _hsl(0);
                const T S = _hsl(1);
                const T L = _hsl(2);

                assert(0. <= H && H <= 360.);
                assert(0. <= S && S <= 1.);
                assert(0. <= L && L <= 1.);

                const T C = ( T(1.) - std::abs( T(2.) * L - T(1.) ) ) * S;
                const T X = C * (T(1.) - std::abs(std::fmod(H / T(60.), T(2.)) - T(1.)));
                const T m = L - C * T(0.5);

                switch(size_t(H / 60.)){
                    case 0: return { C + m, X + m, m};
                    case 1: return { X + m, C + m, m};
                    case 2: return { m, C + m, X + m};
                    case 3: return { m, X + m, C + m};
                    case 4: return { X + m, m, C + m};
                    case 5: return { C + m, m, X + m};
                    default: return { 0., 0., 0.};
                }

            }; //HSL_to_RGB_deg

            //---------------------------------------------------------------------------------------//
            //                                        HSL_to_RGB_rad
            //---------------------------------------------------------------------------------------//

            /*
                h: [0, 2pi]
                s: [0, 1]
                l: [0, 1]
                rgb: [0, 1]
            */
            template <class T>
            EIGEN_STRONG_INLINE constexpr Eigen::Vector3<T> HSL_to_RGB_rad(const Eigen::Matrix<T, 3, 1>& _hsl) noexcept {
                return HSL_to_RGB_deg<T>( { _hsl(0) * T( 180. / M_PI), _hsl(1), _hsl(2) } );
            }; //HSL_to_RGB_rad

            //---------------------------------------------------------------------------------------//
            //                                        rgb_to_gs
            //---------------------------------------------------------------------------------------//

            template <class T>
            EIGEN_STRONG_INLINE constexpr Eigen::Vector3<T> rgb_to_gs(const Eigen::Matrix<T, 3, 1>& _rgb) noexcept {
                const T gs = T(0.299) * _rgb(0) + T(0.587) * _rgb(1) +  T(0.114) * _rgb(2);
                return { gs, gs, gs };
            }; //rgb_to_gs

            //---------------------------------------------------------------------------------------//
            //                                     phi_0
            //---------------------------------------------------------------------------------------//

            template <class T, template<typename> class INVARIANTS>
            EIGEN_STRONG_INLINE std::complex<T> phi_0 (
                size_t _t, 
                const Eigen::VectorX<T>& _x, 
                const INVARIANTS<T>& _inv
            ) noexcept {
                const Eigen::VectorX<std::complex<T>> xq = _x - _inv.q[_t];
                const Eigen::VectorX<std::complex<T>> xqt = xq.transpose();
                const std::complex<T> e1 = _inv.i_2_E_2 * xqt * _inv.P_Q_1[_t] * xq;
                const std::complex<T> e2 = _inv.i_E_2_p[_t] * xq;
                return _inv.pre[_t] * std::exp(e1 + e2);
            }; //phi_0

            //---------------------------------------------------------------------------------------//
            //                                        phi
            //---------------------------------------------------------------------------------------//

            template <class T, template<typename> class INVARIANTS>
            EIGEN_STRONG_INLINE Eigen::VectorX<std::complex<T>> phi (
                size_t _t,
                const Eigen::VectorX<T>& _x,
                const tsl::robin_map<Eigen::Index, std::complex<T>>& _phis,
                const Eigen::VectorX<T>& _index,
                const INVARIANTS<T>& _inv
            ) noexcept {

                using Index = Eigen::VectorX<Eigen::Index>;
                using Vector = Eigen::VectorX<std::complex<T>>;
                using Matrix = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>;

                const auto xq = _x - _inv.q[_t];

                Vector res (_inv.dimensions);

                Vector kp (_inv.dimensions);
                for (size_t j = 0; j < _inv.dimensions; ++j) {
                    if (_index(j) - 1 < 0) {
                        kp(j) = { 0., 0. };
                        continue;
                    }

                    Index k_1 = _index;
                    k_1(j) = std::max<Eigen::Index>(k_1(j) - 1, 0ll);
                    const Eigen::Index ii = Detail::index(k_1, _inv.k);
                    assert(_phis.contains(ii));
                    const std::complex<double>& p = (*_phis.find(ii)).second;
                    kp(j) = std::sqrt(_index(j)) * p;  
                }

                const Eigen::Index ii = Detail::index(_index, _inv.k);
                assert(_phis.contains(ii));
                const std::complex<double>& p = (*_phis.find(ii)).second;
                auto phi_t = _inv.Q_1[_t] * xq * p - _inv.Q_1_Q_T[_t] * kp;

                Vector phi (_inv.dimensions);
                for(size_t i = 0; i < _inv.dimensions; ++i){
                    const T sk = std::sqrt(T(_index(i)) + 1.);
                    const std::complex<T> skc = std::complex<T>(sk, 0.);
                    phi(i) = phi_t(i) / skc;
                }

                return phi;
            }; //phi

            //---------------------------------------------------------------------------------------//
            //                                  compute_cube
            //---------------------------------------------------------------------------------------//

            template <class T, template<typename> class INVARIANTS>
            EIGEN_STRONG_INLINE tsl::robin_map<Eigen::Index, std::complex<T>> compute_cube (
                size_t _t, 
                const Eigen::VectorX<T>& _x, 
                const INVARIANTS<T>& _inv
            ) {

                using Index = Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>;
                using Vector = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>;
                using Matrix = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>;

                const size_t dim = _inv.dimensions;

                size_t size = _inv.k(0) + 1;
                for (size_t i = 1; i < _inv.k.size(); ++i)
                    size *= (_inv.k(i)+1);

                tsl::robin_map<Eigen::Index, std::complex<T>> phis;

                //iterate over ks
                Index index(dim);
                index.fill(0);

                bool first = true;

                while (true) {
                    for (index(dim-1) = 0; index(dim-1) <= _inv.k(dim-1); ++index(dim-1)) { 
                        if (first) {
                            first = false;
                            const auto phi0 = Detail::phi_0(_t, _x, _inv);
                            phis.insert( {0, phi0} );
                            --index(dim-1);
                            continue;
                        }

                        //check if have reached the end of the shape
                        const Eigen::Index ii = Detail::index(index, _inv.k);
                        if(!_inv.k_shape.contains(ii))
                            break;

                        //compute phi for index
                        const auto phi = Detail::phi(_t, _x, phis, index, _inv);

                        for (size_t d = 0; d < dim; ++d) {
                            Index ni = index;
                            ni(d) += 1;
                            const Eigen::Index ii = Detail::index(ni, _inv.k);
                            phis.insert( {ii, phi(d)} );
                            //phis[ii] = phi(d);
                        }
                    }

                    bool done = false;
                    for (Eigen::Index d = dim - 2; d >= 0; --d) {
                        index(d) += 1;
                        if (index(d) >=  _inv.k(d)) {
                            if (d == 0)
                                done = true;
                            else
                                index(d) = 0;
                        } else
                            break;
                    }
                    if (done) break;
                }

                return phis;
            } //compute

            //---------------------------------------------------------------------------------------//
            //                                   linear_combination
            //---------------------------------------------------------------------------------------//

            template <class T, template<typename> class INVARIANTS>
            std::complex<T> linear_combination(
                size_t _time_step,
                const tsl::robin_map<Eigen::Index, std::complex<T>>& _phis,
                const INVARIANTS<T>& _inv
            ) {
                std::complex<T> res (0., 0.);
                for(Eigen::Index k = 0; k < _inv.Ks.size(); ++k){ 
                    const Eigen::Index idx = Detail::index( _inv.Ks[k],  _inv.k_max);
                    assert(_phis.count(idx) != 0);
                    const std::complex<T>& p = (*_phis.find(idx)).second;
                    res +=  _inv.c_0[_time_step](k) * p;
                } 
                res *= std::exp(std::complex<double>(0., 1.) *  _inv.S(Eigen::Index(_time_step)));
            }

        }//Detail

        //---------------------------------------------------------------------------------------//
        //                                        HyperCube
        //---------------------------------------------------------------------------------------//

        template <class T>
        class HyperCube {

            std::complex<T> value;

        public:

            template<template<typename> class INVARIANTS>
            HyperCube(size_t _time_step, const Eigen::Vector3<T>& _pos, const INVARIANTS<T>& _inv) {
                const auto phis = Detail::compute_cube(_time_step, _pos, _inv);
                value = linear_combination(_time_step, phis, _inv);
            }

            /**
             * @brief Returns the functions value in HSL with max value of 10.
            */
            template<T MAXVALUE = 10.>
            EIGEN_STRONG_INLINE Eigen::Vector3<T> HSL() const noexcept {
                return Detail::c_to_HSL<T, MAXVALUE>(value);
            }

            /**
             * @brief Returns the functions value in RGB
            */
            EIGEN_STRONG_INLINE Eigen::Vector3<T> RGB() const noexcept {
                return Detail::HSL_to_RGB_deg(HSL());
            }

            /**
             * @brief Returns the functions value in Grayscale
            */
            EIGEN_STRONG_INLINE Eigen::Vector3<T> GRAYSCALE() const noexcept {
                return Detail::rgb_to_gs(RGB());
            }

        };//HyperCube

    }//Hagedorn

    //---------------------------------------------------------------------------------------//
    //                                        Sampler
    //---------------------------------------------------------------------------------------//

    class Sampler {

    public:

        Sampler(const json& _config) {


        }

        bool success() {
            return true; 
        }

    };

    //---------------------------------------------------------------------------------------//
    //                                        Renderer
    //---------------------------------------------------------------------------------------//

    class Renderer {

    public:
        Renderer(const json& _config) {

        }

        bool success() {
            return true; 
        }

    };//Renderer

}//Hagedorn

#endif //TURBODORN_HPP
