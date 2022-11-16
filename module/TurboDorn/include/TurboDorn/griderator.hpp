
#ifndef G3D_RayGridIterator_h
#define G3D_RayGridIterator_h

#include <Eigen/Eigen>

#include <cassert>
#include <cmath>

namespace TurboDorn {
    
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
            const AABox gridBounds(Vector3::Zero(), Vector3(numCells) * cellSize);

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

}  // namespace TurboDorn

#endif