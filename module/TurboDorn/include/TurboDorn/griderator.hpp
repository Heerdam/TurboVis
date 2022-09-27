
#ifndef G3D_RayGridIterator_h
#define G3D_RayGridIterator_h

#include <Eigen/Eigen>

namespace TurboDorn {

    namespace Detail {
        
        struct Ray {
            Eigen::Vector3d o, d;
        };

        struct AABB {
            Eigen::Vector3d min, max;
        };

        struct Intersection {
            std::pair<Eigen::Vector3d, Eigen::Vector3d> p;
            std::pair<double, double> t;
        };

        std::optional<Intersection> intersect_ray_aabb(const Ray& _r, const AABB& _b) {
            double tmin = -INFINITY, tmax = INFINITY;
            for (int i = 0; i < 3; i++) {
                const double t1 = (_b.min(i) - _r.o(i)) * _r.d(i);
                const double t2 = (_b.max(i) - _r.o(i)) * _r.d(i);
                tmin = std::max(tmin, std::min(t1, t2));
                tmax = std::min(tmax, std::max(t1, t2));
            }
            Intersection ret;
            if(tmax < tmin) return std::nullopt;     
            ret.t.first = tmin;
            ret.t.second = tmax;
            ret.p.first = _r.o + _r.d * tmin;
            ret.p.second = _r.o + _r.d * tmax;
            return ret;
        };

    } //Detail

    class Griderator {

        using iterator_category = std::forward_iterator_tag;

        /* Extent of the grid in each dimension, in grid cell units.*/
        Eigen::Vector3i m_numCells;

        /* Current grid cell m_index */
        Eigen::Vector3i m_index;

        /* Sign of the direction that the ray moves along each axis; +/-1 or 0 */
        Eigen::Vector3i m_step;

        /* Size of one cell in units of t along each axis. */
        Eigen::Vector3d m_tDelta;

        /* Distance along the ray of the first intersection with the
        current cell (i.e., that given by m_index). Zero for the cell that contains the ray origin. */
        float m_enterDistance;

        /* Distance along the ray to the intersection with the next grid
        cell. enterDistance and m_exitDistance can be used to bracket ray
        ray-primitive intersection tests within a cell.
        */
        Eigen::Vector3d m_exitDistance;

        /* The axis along which the ray entered the cell; X = 0, Y = 1, Z = 2.  This is always set to X for the cell that contains the ray origin. */
        int m_enterAxis;

        /* The original ray */
        Detail::Ray m_ray;

        /* Size of each cell along each axis */
        Eigen::Vector3d m_cellSize;

        /* True if index() refers to a valid cell inside the grid.  This
            is usually employed as the loop termination condition.*/
        bool m_insideGrid;

        /* The value that the index will take on along each boundary when
            it just leaves the grid. */
        Eigen::Vector3i m_boundaryIndex;

        /* True if this cell contains the ray origin */
        bool m_containsRayOrigin;

    public:

        const Detail::Ray& ray() const {
            return m_ray;
        }

        Eigen::Vector3i numCells() const {
            return m_numCells;
        }

        int enterAxis() const {
            return m_enterAxis;
        }

        /* Outward-facing normal to the current grid cell along the
        partition just entered. Initially zero. */
        Eigen::Vector3i enterNormal() const {
            Eigen::Vector3i normal(0, 0, 0);
            normal(m_enterAxis) = -m_step(m_enterAxis);
            return normal;
        }

        const Eigen::Vector3d& cellSize() const {
            return m_cellSize;
        }

        /* Location where the ray entered the current grid cell */
        Eigen::Vector3d enterPoint() const {
            return m_ray.o + m_enterDistance * m_ray.d;
        }

        /* Location where the ray exits the current grid cell */
        Eigen::Vector3d exitPoint() const {
            return m_ray.o + m_exitDistance.minCoeff() * m_ray.d;
        }

        float enterDistance() const {
            return m_enterDistance;
        }

        /* Distance from the ray origin to the exit point in this cell. */
        float exitDistance() const {
            return m_exitDistance.minCoeff();
        }

        const Eigen::Vector3i& step() const {
            return m_step;
        }

        const Eigen::Vector3i& index() const {
            return m_index;
        }

        const Eigen::Vector3d& tDelta() const {
            return m_tDelta;
        }

        bool insideGrid() const {
            return m_insideGrid;
        }

        bool containsRayOrigin() const {
            return m_containsRayOrigin;
        }

        /*
            adapted from https://sourceforge.net/p/g3d/code/HEAD/tree/G3D10/G3D-base.lib/include/G3D-base/RayGridIterator.h

            Computes conservative line raster/voxelization across a grid for
            use in walking a grid spatial data structure or or voxel scene
            searching for intersections.  At each iteration, the iterator
            steps exactly one cell in exactly one dimension.

            Example of this iterator applied to ray-primitive intersection in a
            grid:

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

            Initialize the iterator to the first grid cell hit by
            the ray and precompute traversal variables.

            The grid is assumed to have a corner at (0,0,0) and extend
            along the canonical axes. 

            If the ray never passes through the grid, insideGrid() will
            be false immediately after intialization.

            If using for 2D iteration, set <code>numCells.z = 1 and ray.origin().z = 0.5

            \param cellSize The extent of one cell

            \param _min_bounds_cell_index The location of the lowest corner of grid cell minBoundsCellIndex along each axis.
            This translates the grid relative to the ray's coordinate frame.

            \param _min_bounds_cell_index The index of the first grid cell.  This allows
            operation with grids defined on negative indices.  This translates all
            grid indices.
        */

        /*
        RayGridIterator::RayGridIterator
            (Ray                    ray, 
            const Vector3int32&    numCells, 
            const Vector3&         cellSize,
            const Point3&          gridOrigin,
            const Point3int32&     gridOriginIndex)
        */
        Griderator( Detail::Ray _r,
                    Eigen::Vector3i _num_cells,
                    Eigen::Vector3d _cell_size = Eigen::Vector3d::Ones(),
                    Eigen::Vector3d _min_bounds_location = Eigen::Vector3d::Zero(),
                    Eigen::Vector3i _min_bounds_cell_index = Eigen::Vector3i::Zero()) : 
                        m_numCells(std::move(_num_cells)),
                        m_enterDistance(0.f),
                        m_enterAxis(0),
                        m_ray({_r.o - _min_bounds_location, std::move(_r.d)}),
                        m_cellSize(std::move(_cell_size)),
                        m_insideGrid(true),
                        m_containsRayOrigin(true) 
        {

            //////////////////////////////////////////////////////////////////////
            // See if the ray begins inside the box
            const Detail::AABB gridBounds { Eigen::Vector3d::Zero(), Eigen::Vector3d(_num_cells) * _cell_size };

            bool startsOutside = false;
            bool inside = false;
            Eigen::Vector3d startLocation = _r.o;

            const auto hit = Detail::intersect_ray_aabb(m_ray, { _min_bounds_location, _min_bounds_location + m_cellSize * m_numCells});
            const bool passesThroughGrid = hit.has_value();
            
            /*CollisionDetection::rayAABox (
                                                ray, 
                                                Eigen::Vector3d(1., 1., 1.) / ray.d,
                                                gridBounds, 
                                                gridBounds.center(),
                                                square(gridBounds.extent().length() * 0.5),
                                                startLocation,
                                                inside
                                            );*/
                
            if (!inside) {
                if (passesThroughGrid) {
                    // Back up slightly so that we immediately hit the
                    // start location.  The precision here is tricky--if
                    // the ray strikes at a very glancing angle, we need
                    // to move a large distance along the ray to enter the
                    // grid.  If the ray strikes head on, we only need to
                    // move a small distance.
                    m_enterDistance = (m_ray.o - startLocation).length() - 0.0001;
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
                debugAssert(d >= 0 && d <= cellSize[a]);

                if (ray.direction()[a] != 0) {
                    m_exitDistance[a] = d / fabs(ray.direction()[a]) + m_enterDistance;
                } else {
                    // Ray is parallel to this partition axis.
                    // Avoid dividing by zero, which could be NaN if d == 0
                    m_exitDistance[a] = (double)inf();
                }
            }

            if (gridOriginIndex.nonZero()) {
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

        /* Increment the iterator, moving to the next grid cell */
        Griderator& operator++() {
            // Find the axis of the closest partition along the ray
            if (m_exitDistance.x < m_exitDistance.y) {
                if (m_exitDistance.x < m_exitDistance.z) 
                    m_enterAxis = 0;
                else m_enterAxis = 2;
            } else if (m_exitDistance.y < m_exitDistance.z) 
                m_enterAxis = 1; 
            else m_enterAxis = 2;

            m_enterDistance = m_exitDistance(m_enterAxis);
            m_index(m_enterAxis) += m_step(m_enterAxis);
            m_exitDistance(m_enterAxis) += m_tDelta(m_enterAxis);

            // If the index just hit the boundary exit, we have
            // permanently exited the grid.
            m_insideGrid = m_insideGrid && (m_index(m_enterAxis) != m_boundaryIndex(m_enterAxis));
            m_containsRayOrigin = false;
            return *this;
        }

    };

}//TurboDorn

#endif