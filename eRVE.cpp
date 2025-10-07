#include <vector>
#include <array>
#include <functional>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cassert>
#include <random>
#include <algorithm>
#include <fstream>
#include <string>
#include <numeric>

struct Subdomain {
    int row, col;
    std::array<double, 4> bounds; // xmin, xmax, ymin, ymax

    // Each neighbor: {nrow, ncol, wrap_x, wrap_y, dx, dy, is_wrapped}
    std::vector<std::array<double, 7>> neighbors;
};

std::vector<std::vector<Subdomain>> createSubdomains(double length, double breadth, int nx, int ny) {
    double dx = length / nx;
    double dy = breadth / ny;

    std::vector<std::vector<Subdomain>> grid(ny, std::vector<Subdomain>(nx));

    for (int row = 0; row < ny; ++row) {
        for (int col = 0; col < nx; ++col) {
            Subdomain& sd = grid[row][col];
            sd.row = row;
            sd.col = col;
            sd.bounds = {col * dx, (col + 1) * dx, row * dy, (row + 1) * dy};

            const std::vector<std::array<int, 2>> directions = {
                {-1, -1}, {-1, 0}, {-1, 1},
                { 0, -1}, {0, 0},  { 0, 1},
                { 1, -1}, { 1, 0}, { 1, 1}
            };

            for (const auto& d : directions) {
                int raw_nrow = row + d[0];
                int raw_ncol = col + d[1];

                int nrow = (raw_nrow + ny) % ny;
                int ncol = (raw_ncol + nx) % nx;

                int wrap_x = 0, wrap_y = 0;
                bool is_wrapped = false;

                if (raw_ncol < 0)         { wrap_x = -1; is_wrapped = true; }
                else if (raw_ncol >= nx)  { wrap_x = +1; is_wrapped = true; }

                if (raw_nrow < 0)         { wrap_y = -1; is_wrapped = true; }
                else if (raw_nrow >= ny)  { wrap_y = +1; is_wrapped = true; }

                // Bring neighbor into current frame: negate wrap translation
                double dx_wrap = -wrap_x * length;
                double dy_wrap = -wrap_y * breadth;

                sd.neighbors.push_back({
                    (double)nrow, (double)ncol,
                    (double)wrap_x, (double)wrap_y,
                    -1*dx_wrap, -1*dy_wrap,
                    is_wrapped ? 1.0 : 0.0
                });
            }
        }
    }

    return grid;
}

double computePolygonArea(const std::vector<std::array<double, 2>>& vertices) {
    double A = 0.0;
    size_t n = vertices.size();
    for (size_t i = 0; i < n; ++i) {
        const auto& p = vertices[i];
        const auto& q = vertices[(i + 1) % n];
        A += (p[0] * q[1] - q[0] * p[1]);
    }
    return std::abs(A) * 0.5;
}

enum class ShapeType {
    Circle,
    Polygon
};


class Shape {
public:
    Shape(ShapeType type): type(type) {}

    bool moved = false;

    ShapeType getType() const { return type; }
    std::array<double, 2> getCentroid() const { return centroid; }
    std::vector<std::array<double, 2>> getVertices() const { return vertices; }
    double getRadius() const { return radius; }
    double getArea() const { return area; }
    double lmin;
    int overlap_count;
    int stagnant_moves = 0;
    std::array<double, 4> getBoundingBox() const {
        if (type == ShapeType::Circle) {
            double x = centroid[0];
            double y = centroid[1];
            return {x - radius, x + radius, y - radius, y + radius};
        } else {
            double xmin = std::numeric_limits<double>::max();
            double xmax = std::numeric_limits<double>::lowest();
            double ymin = std::numeric_limits<double>::max();
            double ymax = std::numeric_limits<double>::lowest();

            for (const auto& v : vertices) {
                xmin = std::min(xmin, v[0]);
                xmax = std::max(xmax, v[0]);
                ymin = std::min(ymin, v[1]);
                ymax = std::max(ymax, v[1]);
            }
            return {xmin, xmax, ymin, ymax};
        }
    }

    void setCentroid(const std::array<double, 2>& c) { centroid = c; }
    void setRadius(double r) { radius = r; }

    void setVertices(const std::vector<std::array<double, 2>>& v) {
        vertices = v;
        if (type == ShapeType::Polygon) {
            area = computePolygonArea(vertices);
        }
    }

    // Assign subdomain location
    void assignSubdomain(int row, int col, const std::vector<std::vector<Subdomain>>& grid) {
        subdomain_row = row;
        subdomain_col = col;
        subdomain_bounds = grid[row][col].bounds;
    }


    std::pair<int, int> getSubdomainIndex() const {
        return {subdomain_row, subdomain_col};
    }

    // Custom setter for circle type
    void setCircle(double r, const std::array<double, 2>& c) {
        radius = r;
        centroid = c;
        area = M_PI * r * r;
    }

    // Motion
    void setVelocity(double magnitude, double angle_deg) {
        velocity_mag = magnitude;
        velocity_angle_deg = angle_deg;
    }

    void setRotation(double rotation_per_move_deg) {
        rot_velocity_deg = rotation_per_move_deg;
    }

    // Getters
    double getVelocityMagnitude() const { return velocity_mag; }
    double getVelocityAngle() const { return velocity_angle_deg; }
    double getRotationAngle() const { return rot_velocity_deg; }

    // movement step
    void move() {
        double angle_rad = velocity_angle_deg * M_PI / 180.0;
        double dx_unit = std::cos(angle_rad);
        double dy_unit = std::sin(angle_rad);

        double distance_remaining = velocity_mag;

        const double sub_xmin = subdomain_bounds[0];
        const double sub_xmax = subdomain_bounds[1];
        const double sub_ymin = subdomain_bounds[2];
        const double sub_ymax = subdomain_bounds[3];

        double cx = centroid[0];
        double cy = centroid[1];

        while (distance_remaining > 1e-8) {
            // Compute distances to walls in the direction of motion
            double tx = (dx_unit > 0) ? (sub_xmax - cx) / dx_unit :
                        (dx_unit < 0) ? (sub_xmin - cx) / dx_unit :
                        std::numeric_limits<double>::infinity();

            double ty = (dy_unit > 0) ? (sub_ymax - cy) / dy_unit :
                        (dy_unit < 0) ? (sub_ymin - cy) / dy_unit :
                        std::numeric_limits<double>::infinity();

            double t_min = std::min(tx, ty);

            // If we can travel full remaining distance, do so
            if (t_min >= distance_remaining) {
                cx += dx_unit * distance_remaining;
                cy += dy_unit * distance_remaining;
                distance_remaining = 0;
            } else {
                // Move to boundary and reflect
                cx += dx_unit * t_min;
                cy += dy_unit * t_min;
                distance_remaining -= t_min;

                // Reflect appropriate component
                if (t_min == tx) dx_unit *= -1.0;
                if (t_min == ty) dy_unit *= -1.0;
            }
        }

        // Update centroid
        std::array<double, 2> old_centroid = centroid;
        centroid[0] = cx;
        centroid[1] = cy;

        // Update direction
        velocity_angle_deg = std::atan2(dy_unit, dx_unit) * 180.0 / M_PI;

        // Move vertices if polygon
        if (type == ShapeType::Polygon) {
            std::array<double, 2> delta = {cx - old_centroid[0], cy - old_centroid[1]};
            for (auto& v : vertices) {
                v[0] += delta[0];
                v[1] += delta[1];
            }
        }
    }

        void rotatePolygon(double angle_deg) {
        double angle_rad = angle_deg * M_PI / 180.0;
        double cx = centroid[0];
        double cy = centroid[1];

        for (auto& v : vertices) {
            double x = v[0] - cx;
            double y = v[1] - cy;
            v[0] = cx + x * std::cos(angle_rad) - y * std::sin(angle_rad);
            v[1] = cy + x * std::sin(angle_rad) + y * std::cos(angle_rad);
        }
    }


    void moveToSubdomainRandomly(const Subdomain& sd) {
        static std::mt19937 rng(std::random_device{}());

        std::uniform_real_distribution<double> dist_x(sd.bounds[0], sd.bounds[1]);
        std::uniform_real_distribution<double> dist_y(sd.bounds[2], sd.bounds[3]);

        double new_x = dist_x(rng);
        double new_y = dist_y(rng);

        std::array<double, 2> delta = {new_x - centroid[0], new_y - centroid[1]};
        centroid = {new_x, new_y};

        if (type == ShapeType::Polygon) {
            for (auto& v : vertices) {
                v[0] += delta[0];
                v[1] += delta[1];
            }
        }

        subdomain_row = sd.row;
        subdomain_col = sd.col;
    }

    bool isOutOfBounds() const {
        return !boundary_func(centroid);
    }

private:
    ShapeType type;

    // Motion state
    double rot_velocity_deg = 0.0;
    double velocity_mag = 0.0;
    double velocity_angle_deg = 0.0;

    std::array<double, 2> centroid{};
    std::vector<std::array<double, 2>> vertices;
    double radius = 0.0;
    double area = 0.0;

    // double min_distance = 0.0;
    double vicinity_range = 0.0;
    std::function<bool(const std::array<double, 2>&)> boundary_func;

    int subdomain_row = -1;
    int subdomain_col = -1;
    std::array<double, 4> subdomain_bounds{}; // xmin, xmax, ymin, ymax

};




// === Helper function: Euclidean distance ===
double euclideanDistance(const std::array<double, 2>& p1,
                         const std::array<double, 2>& p2) {
    double dx = p1[0] - p2[0];
    double dy = p1[1] - p2[1];
    return std::sqrt(dx * dx + dy * dy);
}

// Point inside polygon using ray-casting
bool pointInPolygon(const std::array<double, 2>& point,
                    const std::vector<std::array<double, 2>>& polygon) {
    int count = 0;
    size_t n = polygon.size();
    for (size_t i = 0; i < n; ++i) {
        auto& a = polygon[i];
        auto& b = polygon[(i + 1) % n];

        if ((a[1] > point[1]) != (b[1] > point[1])) {
            double x = (b[0] - a[0]) * (point[1] - a[1]) / (b[1] - a[1]) + a[0];
            if (point[0] < x) count++;
        }
    }
    return count % 2 == 1;
}

// Check segment-segment intersection
bool segmentsIntersect(const std::array<double, 2>& p1, const std::array<double, 2>& q1,
                       const std::array<double, 2>& p2, const std::array<double, 2>& q2) {
    auto orientation = [](const std::array<double, 2>& a,
                          const std::array<double, 2>& b,
                          const std::array<double, 2>& c) {
        double val = (b[1] - a[1]) * (c[0] - b[0]) -
                     (b[0] - a[0]) * (c[1] - b[1]);
        if (std::abs(val) < 1e-10) return 0;
        return (val > 0) ? 1 : 2;
    };

    auto onSegment = [](const std::array<double, 2>& p,
                        const std::array<double, 2>& q,
                        const std::array<double, 2>& r) {
        return std::min(p[0], r[0]) <= q[0] && q[0] <= std::max(p[0], r[0]) &&
               std::min(p[1], r[1]) <= q[1] && q[1] <= std::max(p[1], r[1]);
    };

    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    if (o1 != o2 && o3 != o4) return true;
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;
    return false;
}

// Point-segment distance
double pointToSegmentDistance(const std::array<double, 2>& p,
                              const std::array<double, 2>& a,
                              const std::array<double, 2>& b) {
    double dx = b[0] - a[0], dy = b[1] - a[1];
    double t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy);
    t = std::max(0.0, std::min(1.0, t));
    double proj_x = a[0] + t * dx;
    double proj_y = a[1] + t * dy;
    return std::hypot(p[0] - proj_x, p[1] - proj_y);
}

bool BoundingBoxesIntersect(const Shape& a, const Shape& b) {
    auto boxA = a.getBoundingBox();
    auto boxB = b.getBoundingBox();
    double threshold = std::min(a.lmin,b.lmin);

    // Expand both boxes by threshold/2 in all directions
    boxA[0] -= threshold / 2;
    boxA[1] += threshold / 2;
    boxA[2] -= threshold / 2;
    boxA[3] += threshold / 2;

    boxB[0] -= threshold / 2;
    boxB[1] += threshold / 2;
    boxB[2] -= threshold / 2;
    boxB[3] += threshold / 2;

    // Check if they overlap
    double x_overlap = std::max(0.0, std::min(boxA[1], boxB[1]) - std::max(boxA[0], boxB[0]));
    double y_overlap = std::max(0.0, std::min(boxA[3], boxB[3]) - std::max(boxA[2], boxB[2]));

    return x_overlap > 0.0 && y_overlap > 0.0;
}

double boundingBoxOverlapArea(const std::array<double, 4>& a,
                              const std::array<double, 4>& b) {
    double x_overlap = std::max(0.0, std::min(a[1], b[1]) - std::max(a[0], b[0]));
    double y_overlap = std::max(0.0, std::min(a[3], b[3]) - std::max(a[2], b[2]));
    return x_overlap * y_overlap;
}

bool isBoundingBoxInsideAnother(const Shape& inner, const Shape& outer) {
    const auto& box_in = inner.getBoundingBox();  // [xmin, xmax, ymin, ymax]
    const auto& box_out = outer.getBoundingBox();

    bool x_inside = box_in[0] >= box_out[0] && box_in[1] <= box_out[1];
    bool y_inside = box_in[2] >= box_out[2] && box_in[3] <= box_out[3];

    return x_inside && y_inside;
}

// Circle–Circle distance
double distance_circle_circle(const Shape& a, const Shape& b) {
    if (!BoundingBoxesIntersect(a,b)) return 1e12; 
    double center_dist = euclideanDistance(a.getCentroid(), b.getCentroid());
    double boundary_dist = center_dist - a.getRadius() - b.getRadius();
    return boundary_dist;  // Return 0 if overlapping
}

double distance_polygon_polygon(const Shape& a, const Shape& b) {
    if (!BoundingBoxesIntersect(a,b)) return 1e12; 
    const auto& pa = a.getVertices();
    const auto& pb = b.getVertices();

    // Intersection or containment check
    for (size_t i = 0; i < pa.size(); ++i) {
        auto p1 = pa[i];
        auto q1 = pa[(i + 1) % pa.size()];
        for (size_t j = 0; j < pb.size(); ++j) {
            auto p2 = pb[j];
            auto q2 = pb[(j + 1) % pb.size()];
            if (segmentsIntersect(p1, q1, p2, q2)) return 0.0;
        }
    }

    for (const auto& v : pa) if (pointInPolygon(v, pb)) return 0.0;
    for (const auto& v : pb) if (pointInPolygon(v, pa)) return 0.0;

    // Compute min distance from each vertex to each edge
    double min_dist = std::numeric_limits<double>::max();
    for (const auto& v : pa)
        for (size_t j = 0; j < pb.size(); ++j)
            min_dist = std::min(min_dist, pointToSegmentDistance(v, pb[j], pb[(j + 1) % pb.size()]));

    for (const auto& v : pb)
        for (size_t i = 0; i < pa.size(); ++i)
            min_dist = std::min(min_dist, pointToSegmentDistance(v, pa[i], pa[(i + 1) % pa.size()]));

    return min_dist;
}

double distance_circle_polygon(const Shape& circle, const Shape& poly) {
    if (!BoundingBoxesIntersect(circle,poly)) return 1e12; 
    const auto& verts = poly.getVertices();
    const auto& c = circle.getCentroid();
    double r = circle.getRadius();

    for (size_t i = 0; i < verts.size(); ++i) {
        if (pointInPolygon(c, verts)) return 0.0;
    }

    for (size_t i = 0; i < verts.size(); ++i) {
        double d = pointToSegmentDistance(c, verts[i], verts[(i + 1) % verts.size()]);
        if (d <= r) return 0.0;
    }

    // No overlap
    double min_dist = std::numeric_limits<double>::max();
    for (size_t i = 0; i < verts.size(); i++) {
        double d = pointToSegmentDistance(c, verts[i], verts[(i + 1) % verts.size()]);
        min_dist = std::min(min_dist, d);
    }

    return std::max(0.0, min_dist - r);
}


//Master dispatcher
double distanceBetweenShapes(const Shape& a, const Shape& b,
                             const std::vector<std::vector<Subdomain>>& grid) {
    auto [rowA, colA] = a.getSubdomainIndex();
    auto [rowB, colB] = b.getSubdomainIndex();

    const Subdomain& sdA = grid[rowA][colA];

    // Step 2: Check if B is a neighbor of A
    bool found = false;
    double dx = 0.0, dy = 0.0;

    for (const auto& nb : sdA.neighbors) {
        int nrow = static_cast<int>(nb[0]);
        int ncol = static_cast<int>(nb[1]);
        if (nrow == rowB && ncol == colB) {
            dx = nb[4];  // translation x
            dy = nb[5];  // translation y
            found = true;
            break;
        }
    }

    if (!found) {
        return 1e9;  // Not a neighbor — return large value
    }

    // Step 4: Translate shape B into A's frame
    Shape b_translated = b;  // Copy
    std::array<double, 2> new_center = b.getCentroid();
    new_center[0] += dx;
    new_center[1] += dy;
    b_translated.setCentroid(new_center);

    if (b.getType() == ShapeType::Polygon) {
        auto verts = b.getVertices();
        for (auto& v : verts) {
            v[0] += dx;
            v[1] += dy;
        }
        b_translated.setVertices(verts);
    }

    if (a.getType() == ShapeType::Circle && b.getType() == ShapeType::Circle)
        return distance_circle_circle(a, b_translated);
    else if (a.getType() == ShapeType::Polygon && b.getType() == ShapeType::Polygon)
        return distance_polygon_polygon(a, b_translated);
    else if (a.getType() == ShapeType::Circle && b.getType() == ShapeType::Polygon)
        return distance_circle_polygon(a, b_translated);
    else if (a.getType() == ShapeType::Polygon && b.getType() == ShapeType::Circle)
        return distance_circle_polygon(b_translated, a); // flip
    else{
        std::cout<<"distance calculation failed\n";
        return -1;
    }
        
}


double evaluateShape(Shape& target,
                     const std::vector<Shape>& all_shapes,
                     const std::vector<std::vector<Subdomain>>& grid,
                     const double& min_distance1,
                     const double& min_step,
                     const double& max_step) {
    double total_depth = 0.0;
    target.overlap_count = 0;

    double min_distance = min_distance1;

    auto [rowA, colA] = target.getSubdomainIndex();
    const Subdomain& sdA = grid[rowA][colA];

    std::array<double, 2> resultant = {0.0, 0.0};

        double local_depth = 0.0;
        double local_rx = 0.0;
        double local_ry = 0.0;

        for (size_t i = 0; i < all_shapes.size(); ++i) {
            const auto& other = all_shapes[i];
            if (&target == &other) continue;

            auto [rowB, colB] = other.getSubdomainIndex();

            bool found = false;
            double dx = 0.0, dy = 0.0;

            for (const auto& nb : sdA.neighbors) {
                int nrow = static_cast<int>(nb[0]);
                int ncol = static_cast<int>(nb[1]);
                if (nrow == rowB && ncol == colB) {
                    dx = nb[4];
                    dy = nb[5];
                    found = true;
                    break;
                }
            }

            if (!found) continue;
            if(other.lmin==0||target.lmin==0) min_distance = 0;

            auto boxA = target.getBoundingBox();
            boxA[0] -= min_distance / 2;
            boxA[1] += min_distance / 2;
            boxA[2] -= min_distance / 2;
            boxA[3] += min_distance / 2;

            auto boxB = other.getBoundingBox();
            boxB[0] += dx - min_distance / 2;
            boxB[1] += dx + min_distance / 2;
            boxB[2] += dy - min_distance / 2;
            boxB[3] += dy + min_distance / 2;

            double intersection = boundingBoxOverlapArea(boxA, boxB);
            
            if (intersection > 0.0 && distanceBetweenShapes(target, other, grid) < min_distance) {
                // double dist = distanceBetweenShapes(target, other, grid);
                local_depth += intersection;

                std::array<double, 2> pA = target.getCentroid();
                std::array<double, 2> pB = other.getCentroid();
                pB[0] += dx;
                pB[1] += dy;
                
                // local_depth += dist;

                local_rx += (-pB[0] + pA[0]) * intersection;
                local_ry += (-pB[1] + pA[1]) * intersection;
                total_depth += local_depth;
                resultant[0] += local_rx;
                resultant[1] += local_ry;
                target.overlap_count++;
                if(isBoundingBoxInsideAnother(target,other)){
                    const auto& c1 = target.getCentroid();
                    const auto& c2 = other.getCentroid();
                    double dx = c1[0] - c2[0];
                    double dy = c1[1] - c2[1];
                    double centroid_dist = std::sqrt(dx * dx + dy * dy);

                    total_depth *= (1+centroid_dist*0.00001);

                }
                
                
            }
        }

        
    // Update direction
    double mag = target.getVelocityMagnitude();
    if (total_depth > 0.0 && (resultant[0] != 0.0 /*|| resultant[1] != 0.0*/) && target.stagnant_moves<=10) {
        double angle = std::atan2(resultant[1], resultant[0]) * 180.0 / M_PI;
        target.setVelocity(mag, angle);
    } else {
        static thread_local std::mt19937 rng(std::random_device{}());
        std::uniform_real_distribution<double> angle_dist(0.0, 360.0);
        target.setVelocity(mag, angle_dist(rng));
    }
    return total_depth;
}





std::vector<Shape> poolOfShapes(
    int num_circles,
    int num_polygons,
    int num_squares,
    int num_ellipses,
    int num_stars,
    double circle_radius,
    double polygon_radius,
    double square_side,
    std::pair<double, double> ellipse_axes,  // (a, b)
    double star_radius,
    double lmin) {
    std::vector<Shape> shapes;
    std::mt19937 rng(std::random_device{}());

    std::uniform_real_distribution<double> angle_dist(0.0, 360.0);     // degrees
    std::uniform_real_distribution<double> speed_dist(0.1, 1);    // adjust as needed

    auto assign_random_velocity = [&](Shape& s) {
        double angle = angle_dist(rng);
        double speed = speed_dist(rng);
        s.setVelocity(speed, angle);
        s.lmin = lmin;
    };

    // Circles
    for (int i = 0; i < num_circles; ++i) {
        Shape s(ShapeType::Circle);
        s.setCircle(circle_radius, {0.0, 0.0});
        assign_random_velocity(s);
        shapes.push_back(s);
    }

    // Hexagons
    for (int i = 0; i < num_polygons; ++i) {
        Shape s(ShapeType::Polygon);
        std::vector<std::array<double, 2>> verts;
        int n = 6;
        for (int j = 0; j < n; ++j) {
            double angle = 2 * M_PI * j / n;
            verts.push_back({polygon_radius * std::cos(angle), polygon_radius * std::sin(angle)});
        }
        s.setCentroid({0.0, 0.0});
        s.setVertices(verts);
        s.rotatePolygon(angle_dist(rng));
        assign_random_velocity(s);
        shapes.push_back(s);
    }

    // Squares
    for (int i = 0; i < num_squares; ++i) {
        Shape s(ShapeType::Polygon);
        double a = square_side / 2.0;
        std::vector<std::array<double, 2>> verts = {
            {-a, -a}, {a, -a}, {a, a}, {-a, a}
        };
        s.setCentroid({0.0, 0.0});
        s.setVertices(verts);
        s.rotatePolygon(angle_dist(rng));
        assign_random_velocity(s);
        shapes.push_back(s);
    }

    // Ellipses (approximated)
    for (int i = 0; i < num_ellipses; ++i) {
        Shape s(ShapeType::Polygon);
        std::vector<std::array<double, 2>> verts;
        int n = 40;
        double a = ellipse_axes.first;
        double b = ellipse_axes.second;
        for (int j = 0; j < n; ++j) {
            double angle = 2 * M_PI * j / n;
            verts.push_back({a * std::cos(angle), b * std::sin(angle)});
        }
        s.setCentroid({0.0, 0.0});
        s.setVertices(verts);
        s.rotatePolygon(angle_dist(rng));
        assign_random_velocity(s);
        shapes.push_back(s);
    }

    // 5-pointed stars
    for (int i = 0; i < num_stars; ++i) {
        Shape s(ShapeType::Polygon);
        std::vector<std::array<double, 2>> verts;
        int n = 5;
        double r_outer = star_radius;
        double r_inner = r_outer * 0.5;
        for (int j = 0; j < 2 * n; ++j) {
            double angle = M_PI * j / n;
            double r = (j % 2 == 0) ? r_outer : r_inner;
            verts.push_back({r * std::cos(angle), r * std::sin(angle)});
        }
        s.setCentroid({0.0, 0.0});
        s.setVertices(verts);
        s.rotatePolygon(angle_dist(rng));
        assign_random_velocity(s);
        shapes.push_back(s);
    }

    return shapes;
}

std::vector<int> histogram_based_sampling(int N, int total,
                                          int minvf, int maxvf,
                                          double avgvf, double confid,
                                          double x_percent, double y_percent, double z_percent) {
    std::mt19937 gen(std::random_device{}());
    const int num_bins = 1000;
    std::vector<int> bin_counts(num_bins, 0);
    std::vector<int> result;

    // Step 1: Build bin edges
    double bin_width = static_cast<double>(maxvf - minvf) / num_bins;
    std::vector<double> bin_edges(num_bins + 1);
    for (int i = 0; i <= num_bins; ++i)
        bin_edges[i] = minvf + i * bin_width;

    // Step 2: Identify bin regions
    int low_end = static_cast<int>((avgvf - confid - minvf) / bin_width);
    int high_start = static_cast<int>((avgvf + confid - minvf) / bin_width);
    low_end = std::clamp(low_end, 0, num_bins);
    high_start = std::clamp(high_start, 0, num_bins);

    int n = low_end;
    int o = high_start - low_end;
    int m = num_bins - high_start;

    // Step 3: Compute target number of samples per region
    int n_count = static_cast<int>(std::round((x_percent / 100.0) * N));
    int o_count = static_cast<int>(std::round((z_percent / 100.0) * N));
    int m_count = N - n_count - o_count;  // Remaining goes to m region

    // Step 4: Randomly sample values from each bin region
    auto sample_from_bins = [&](int bin_start, int bin_end, int count) {
        std::uniform_real_distribution<> dist(0.0, 1.0);
        for (int i = 0; i < count; ++i) {
            int bin_idx = std::uniform_int_distribution<>(bin_start, bin_end - 1)(gen);
            double value = bin_edges[bin_idx] + dist(gen) * bin_width;
            result.push_back(static_cast<int>(std::round(value)));
        }
    };

    sample_from_bins(0, low_end, n_count);
    sample_from_bins(low_end, high_start, o_count);
    sample_from_bins(high_start, num_bins, m_count);

    // Step 5: Adjust result to match total sum
    int current_sum = std::accumulate(result.begin(), result.end(), 0);
    int diff = total - current_sum;
    while (diff != 0) {
        std::uniform_int_distribution<> pick(0, N - 1);
        int i = pick(gen);
        if (diff > 0 && result[i] < maxvf) {
            result[i]++;
            diff--;
        } else if (diff < 0 && result[i] > minvf) {
            result[i]--;
            diff++;
        }
    }

    return result;
}

void assignShapesToSubdomains(
    std::vector<Shape>& shapes,
    const std::vector<std::vector<Subdomain>>& grid,
    int min_shapes = 0,
    int max_shapes = 41,
    double confid = 5.0,
    double x_percent = 15.0,
    double y_percent = 60.0,
    double z_percent = 25.0
)
 {
    int ny = grid.size();
    int nx = grid[0].size();
    int num_subdomains = nx * ny;
    int total_shapes = shapes.size();
    double avgvf = static_cast<double>(total_shapes) / num_subdomains;


    // Call histogram-based sampling to assign shape count per subdomain
    std::vector<int> shape_counts = histogram_based_sampling(
        num_subdomains, total_shapes,
        min_shapes, max_shapes,
        avgvf, confid,
        x_percent, y_percent, z_percent
    );
    int assigned = std::accumulate(shape_counts.begin(), shape_counts.end(), 0);
    if (assigned < shapes.size()) {
        std::cerr << "ERROR: Not enough subdomain capacity for all shapes!\n";
        std::cerr << "Shapes: " << shapes.size() << ", total assigned: " << assigned << "\n";
        std::exit(EXIT_FAILURE);
    }
    int t=0;
    for (size_t i = 0; i < shape_counts.size(); ++i) {
        std::cout << "Subdomain " << i << ": " << shape_counts[i] << "\n";
        t=t+shape_counts[i];
    }
    // --- Correct shape count mismatch ---
    int current_sum = std::accumulate(shape_counts.begin(), shape_counts.end(), 0);
    int diff = total_shapes - current_sum;

    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<> pick(0, shape_counts.size() - 1);

    while (diff != 0) {
        int i = pick(gen);
        if (diff > 0 && shape_counts[i] < max_shapes) {
            shape_counts[i]++;
            diff--;
        } else if (diff < 0 && shape_counts[i] > min_shapes) {
            shape_counts[i]--;
            diff++;
        }
    }


    // === Assign shapes according to shape_counts ===
    std::random_device rd;
    std::mt19937 rng(rd());
    std::shuffle(shapes.begin(), shapes.end(), rng);
    std::shuffle(shape_counts.begin(), shape_counts.end(), rng);

    int shape_idx = 0;
    std::vector<Shape> placed_shapes;

    for (int r = 0; r < ny; ++r) {
        for (int c = 0; c < nx; ++c) {
            int idx = r * nx + c;
            int count = shape_counts[idx];
            const Subdomain& sd = grid[r][c];

            for (int i = 0; i < count && shape_idx < total_shapes; ++i) {
                Shape& shape = shapes[shape_idx++];
                shape.moveToSubdomainRandomly(sd);
                shape.assignSubdomain(r, c, grid);
            }
        }
    }


    // Optional: Compute and print CV
    double mean = static_cast<double>(total_shapes) / num_subdomains;
    double sq_sum = 0.0;
    for (int n : shape_counts) sq_sum += (n - mean) * (n - mean);
    double stddev = std::sqrt(sq_sum / num_subdomains);
    double cv = stddev / mean;


    std::cout << "Finished assignment with histogram-based distribution (total,shapes = " << t << ","<<total_shapes<<")\n";
}

void assignShapesUniformly(
    std::vector<Shape>& shapes,
    const std::vector<std::vector<Subdomain>>& grid
) {
    int ny = grid.size();
    int nx = grid[0].size();
    int num_subdomains = nx * ny;
    int total_shapes = shapes.size();

    // --- Determine uniform count per subdomain ---
    int base_count = total_shapes / num_subdomains;
    int remainder = total_shapes % num_subdomains;

    std::vector<int> shape_counts(num_subdomains, base_count);
    
    // Distribute remaining shapes randomly across subdomains
    std::mt19937 gen(std::random_device{}());
    std::vector<int> indices(num_subdomains);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), gen);
    for (int i = 0; i < remainder; ++i)
        shape_counts[indices[i]]++;

    std::shuffle(shapes.begin(), shapes.end(), gen);

    // === Assign shapes according to shape_counts ===
    int shape_idx = 0;
    for (int r = 0; r < ny; ++r) {
        for (int c = 0; c < nx; ++c) {
            int idx = r * nx + c;
            int count = shape_counts[idx];
            const Subdomain& sd = grid[r][c];

            for (int i = 0; i < count && shape_idx < total_shapes; ++i) {
                Shape& shape = shapes[shape_idx++];
                shape.moveToSubdomainRandomly(sd);
                shape.assignSubdomain(r, c, grid);
            }
        }
    }

    std::cout << "Uniform assignment done. Avg = " << (double)total_shapes / num_subdomains << ", Max = "
              << *std::max_element(shape_counts.begin(), shape_counts.end()) << ", Min = "
              << *std::min_element(shape_counts.begin(), shape_counts.end()) << "\n";
}

void untangleShapes(std::vector<Shape>& shapes,
                    const std::vector<std::vector<Subdomain>>& grid,
                    // int max_iterations = 1000,
                    double tolerance = 1e-5, //min dist
                    double min_step = 1e-3,
                    double max_step = 0.05) {

    int tot_iter = 0;
    int comp_iter = 0;
    for (int iter = 0; iter < 1e20; ++iter) {
        double total_depth = 0.0;
        int successful_moves = 0;
        double global_depth = 0;
        tot_iter = tot_iter + iter;
        for (auto& shape : shapes) {
            double original_depth = evaluateShape(shape, shapes, grid, tolerance, min_step, max_step);
            total_depth += original_depth;
            comp_iter += iter;
            std::mt19937 rng(std::random_device{}());
            std::uniform_real_distribution<double> angle(-5, 5);

            // Backup state
            auto old_centroid = shape.getCentroid();
            auto old_vertices = shape.getVertices();
            double old_velocity = shape.getVelocityMagnitude();

            // Attempt move
            shape.move();
            // shape.rotatePolygon(angle(rng)); 

            // Compute new overlap
            double new_depth = evaluateShape(shape, shapes, grid, tolerance, min_step, max_step);

            if (new_depth <= original_depth || new_depth == 0 ) {
                // Accept move and increase velocity
                shape.setVelocity(std::min(max_step,old_velocity * 1.1), shape.getVelocityAngle());
                successful_moves++;
                shape.moved = true;
                shape.stagnant_moves = 0;
                global_depth = global_depth + new_depth;
            } else {
                // Revert move and reduce velocity
                shape.setCentroid(old_centroid);
                shape.setVertices(old_vertices);
                shape.moved = false;
                shape.stagnant_moves++;
                shape.setVelocity(std::max(min_step,old_velocity * 0.9), shape.getVelocityAngle());
                global_depth = global_depth + original_depth;
            }
        }
        if (global_depth == 0) break;

        std::cout << "Iteration " << iter + 1
                  << " | Total overlap depth: " << total_depth
                  << " | Successful moves: " << successful_moves << "\n";

    }
    std::cout << "Completed in "<< tot_iter << " iterations.\n";
    std::cout << "Completed in "<< comp_iter << " comprehensive iterations.\n";

}

void saveConfiguration(const std::vector<Shape>& shapes,
                       const std::string& filename = "shape_configuration.txt") {
    std::ofstream fout(filename);
    if (!fout) {
        std::cerr << "Error opening file: " << filename << "\n";
        return;
    }

    fout << std::fixed << std::setprecision(6);  // Format float precision
    fout << "Shape Configuration:\n";

    for (size_t i = 0; i < shapes.size(); ++i) {
        const auto& s = shapes[i];
        auto [r, c] = s.getSubdomainIndex();
        auto bb = s.getBoundingBox();
        auto center = s.getCentroid();

        fout << center[0] << ", " << center[1] << "," << s.getRadius();

        if (s.getType() == ShapeType::Polygon) {
            const auto& verts = s.getVertices();
            fout << " | Vertices: ";
            for (const auto& v : verts) {
                fout << "(" << v[0] << "," << v[1] << ") ";
            }
        }

        fout << "\n";
    }

    fout.close();
    std::cout << "Configuration saved to " << filename << "\n";
}

std::array<int, 2> getBoundaryCrossingDomain(const Shape& shape,
                                                const double& Lx, 
                                                const double& Ly) {

    auto bb = shape.getBoundingBox();  // [xmin, xmax, ymin, ymax]

    std::array<int, 2> direction = {0, 0};

    if (bb[0] < 0) direction[0] += 1;    // left
    if (bb[1] > Lx) direction[0] += -1;   // right
    if (bb[2] < 0) direction[1] += 1;    // bottom
    if (bb[3] > Ly) direction[1] += -1;   // top

    return direction;
}

void saveConfigurationWithPeriodic(const std::vector<Shape>& original_shapes,
                                //    const std::vector<std::vector<Subdomain>>& grid,
                                   double length,
                                   double breadth,
                                   const std::string& filename = "shape_configuration.txt") {
    std::vector<Shape> all_shapes = original_shapes;

    for (const auto& shape : original_shapes) {
        std::array<int, 2> shift = getBoundaryCrossingDomain(shape, length, breadth);

        // Generate shifted copies based on direction
        std::vector<std::array<int, 2>> shifts_to_add;

        if (shift[0] != 0 && shift[1] != 0) {
            // corner crossing (e.g., top-right)
            shifts_to_add.push_back({shift[0], 0});
            shifts_to_add.push_back({0, shift[1]});
            shifts_to_add.push_back(shift);
        } else if (shift[0] != 0 || shift[1] != 0) {
            // edge crossing
            shifts_to_add.push_back(shift);
        }

        for (const auto& s : shifts_to_add) {
            Shape copy = shape;

            std::array<double, 2> c = copy.getCentroid();
            c[0] += s[0] * length;
            c[1] += s[1] * breadth;
            copy.setCentroid(c);

            if (copy.getType() == ShapeType::Polygon) {
                auto verts = copy.getVertices();
                for (auto& v : verts) {
                    v[0] += s[0] * length;
                    v[1] += s[1] * breadth;
                }
                copy.setVertices(verts);
            }

            all_shapes.push_back(copy);
        }
    }

    // Now call the original save function on expanded shape list
    saveConfiguration(all_shapes, filename);
}


int main() {
    // Domain and grid setup
    double Lx = 0.031330, Ly = 0.03133;
    double vf = 0.60;
    int nx = 4, ny = 4;
    double ratio = 50; // Lx/radius
    double radius = 0.0025; Lx/ratio;
    double min_distance = radius * 0.07;
    int number = 30;// vf*Lx*Ly/M_PI/radius/radius;
    // radius = std::sqrt(vf*Lx*Ly/M_PI/(static_cast<double>(number)));
    double min_step = min_distance*0.01;
    double max_step = Lx/nx; //square setting
    double void_vf = 0.01;
    double void_radius = radius * 0.2;
    int void_number = void_vf*Lx*Ly/M_PI/void_radius/void_radius;
    void_radius = std::sqrt(void_vf*Lx*Ly/M_PI/(static_cast<double>(void_number)));
    int min_shapes = 0;
    int max_shapes = 45;
    double min_void_fraction = 0;
    double max_void_fraction = 0.05;
    int min_void_number = 0;
    int max_void_number = void_number/nx/ny*2;
    

    // Create subdomain grid
    auto grid = createSubdomains(Lx, Ly, nx, ny);

    // Define boundary function
    auto boundary = [=](const std::array<double, 2>& pt) {
        return pt[0] >= 0.0 && pt[0] <= Lx && pt[1] >= 0.0 && pt[1] <= Ly;
    };

    // Create a pool of shapes
    std::vector<Shape> shapes = poolOfShapes(
        number,  // circles
        100*0,  // polygons
        400*0,  // squares
        100*0,  // ellipses
        0,  // stars
        radius, 0.1, radius*std::sqrt(M_PI),           // circle, polygon, square side
        {0.15, 0.1},                // ellipse radii
        0.1,                       // star radius
        min_distance
    );
    // assignShapesToSubdomains(shapes, grid, min_shapes, max_shapes);
    assignShapesUniformly(shapes , grid);
    std::cout << "Calling untangleShapes on " << shapes.size() << " shapes\n";

    

    // std::vector<Shape> void_shapes = poolOfShapes(
    //     void_number,  // circles
    //     100*0,  // polygons
    //     200*0,  // squares
    //     100*0,  // ellipses
    //     0,  // stars
    //     void_radius, 0.1, 0.15,           // circle, polygon, square sizes
    //     {0.15, 0.1},                // ellipse radii
    //     0.1,                       // star radius
    //     min_distance * 0.000001
    // );
    //     assignShapesToSubdomains(void_shapes, grid, min_void_number, max_void_number);  

    // std::cout<<"void_shapes:"<<void_shapes.size()<<"," << "min_void_number:"<<min_void_number<<","<< "max_void_number:"<<max_void_number<<"\n";
    
    // for (size_t i = 0; i < shapes.size(); ++i) {
    //     auto [r, c] = shapes[i].getSubdomainIndex();
    //     std::cout << "Shape " << i << " -> Subdomain (" << r << ", " << c << ")\n";
    // }
    
    std::vector<Shape> all_shapes;
    // all_shapes.reserve(shapes.size() + void_shapes.size());
    // all_shapes.insert(all_shapes.end(), shapes.begin(), shapes.end());
    // all_shapes.insert(all_shapes.end(), void_shapes.begin(), void_shapes.end());
    all_shapes.reserve(shapes.size());
    all_shapes.insert(all_shapes.end(), shapes.begin(), shapes.end());


    saveConfiguration(all_shapes, "initial.txt");
    
    untangleShapes(all_shapes,
                    grid,
                    // 10000,
                    min_distance,
                    min_step,
                    max_step);
    
    saveConfigurationWithPeriodic(all_shapes, Lx, Ly, "final.txt");
    return 0;
}