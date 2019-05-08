#include "pybind11/pybind11.h"
#include <utility>
#include "cartesian_kernel.h"
#include "solve_lp.h"
#include "pybind11/numpy.h"
#include <eigen3/Eigen/Eigen>
#include "polytopes.h"
#include <iostream>
#include "polytope_generators.h"
#include <vector>
#include "samplers.h"
#include <memory>

typedef float                    NT;
typedef Cartesian<NT>             Kernel;
typedef typename Kernel::Point    Point;
typedef HPolytope<Point> Polytope;
typedef boost::mt19937            RNGType;

namespace py=pybind11;

PYBIND11_MODULE(volesti, m) {
    py::class_<Point>(m, "Point", py::buffer_protocol())
        .def(py::init([](py::array_t<NT> b) {
            py::buffer_info info = b.request();
            Point* p = new Point(info.shape[0]);
            p->set_dimension(info.shape[0]);
            NT* buffer_data = static_cast<NT*>(info.ptr);
            for (int i=0; i<p->dimension(); i++) {
                p->set_coord(i, buffer_data[i]);
            }
            return p;
        }))
        .def_buffer([](Point &p) -> py::buffer_info {
            return py::buffer_info(
                p.data(),
                sizeof(NT),
                py::format_descriptor<NT>::format(),
                1,
                { p.dimension() },
                { sizeof(NT) }
            );
        })
        .def("__add__", [](Point &p, Point &p2) {
            return p+p2;
        }, py::is_operator());

    py::class_<Polytope>(m, "Polytope")
        .def(py::init([](py::array_t<NT> A, py::array_t<NT> b) {
            py::buffer_info A_info = A.request();
            py::buffer_info b_info = b.request();

            Eigen::Map<Polytope::MT> A_map((NT*) A_info.ptr, A_info.shape[0], A_info.shape[1]);
            Eigen::Map<Polytope::VT> b_map((NT*) b_info.ptr, b_info.shape[0], 1);
            Polytope* p = new Polytope(A_map, b_map);
            return p;
        }))
        .def("print", &Polytope::print)
        .def("create_point_representation", [](Polytope &p, py::array_t<NT> internalPoint, int num_bits) {
            py::buffer_info ip_info = internalPoint.request();
            p.create_point_representation((NT*) ip_info.ptr, num_bits);
        }, py::arg(), py::arg("num_bits") = 0)
        .def("contains_point", [](Polytope &p, py::array_t<NT> point) -> bool {
            py::buffer_info point_info = point.request();
            return p.contains_point((NT* ) point_info.ptr);
        })
        .def("is_in", [](Polytope &p, py::array_t<NT> point) -> bool {
            py::buffer_info point_info = point.request();
            Point the_point(point_info.shape[0]);
            the_point.set_dimension(point_info.shape[0]);
            NT* buffer_data = static_cast<NT*>(point_info.ptr);
            for (int i=0; i<the_point.dimension(); i++) {
                the_point.set_coord(i, buffer_data[i]);
            }
            return p.is_in(the_point)==-1;
        })
        .def("intersect", [](Polytope &p, py::array_t<NT> source, py::array_t<NT> direction, int facet_idx) -> py::array_t<NT> {
            py::buffer_info source_i = source.request();
            py::buffer_info ray_i = direction.request();

            Point source_p(p.dimension(), (NT* ) source_i.ptr);
            Point ray_p(p.dimension(), (NT* ) ray_i.ptr);

            Point intersection = p.intersect_ray_hyperplane(source_p, ray_p, facet_idx);

            return py::array_t<NT>(
                py::buffer_info(
                    intersection.data(),
                    sizeof(NT),
                    py::format_descriptor<NT>::format(),
                    1,
                    { p.dimension() },
                    { sizeof(NT) }
                )
            );
        })
        .def("boundary", [](Polytope &p, py::array_t<NT> source, py::array_t<NT> direction, float epsilon, int maxSteps) -> py::array_t<NT> {
            py::buffer_info source_i = source.request();
            py::buffer_info ray_i = direction.request();

            Point source_p(p.dimension(), (NT* ) source_i.ptr);
            Point ray_p(p.dimension(), (NT* ) ray_i.ptr);

            Point intersection = p.compute_boundary_intersection(source_p, ray_p, epsilon, maxSteps);

            return py::array_t<NT>(
                py::buffer_info(
                    intersection.data(),
                    sizeof(NT),
                    py::format_descriptor<NT>::format(),
                    1,
                    { p.dimension() },
                    { sizeof(NT) }
                )
            );
        })
        .def("line_intersect", [](Polytope &p, py::array_t<NT> source, py::array_t<NT> direction) -> py::array_t<NT> {
            py::buffer_info source_i = source.request();
            py::buffer_info ray_i = direction.request();

            Point source_p(p.dimension(), (NT* ) source_i.ptr);
            Point ray_p(p.dimension(), (NT* ) ray_i.ptr);

            std::pair<Point, Point> intersection = p.line_intersect_as_p(source_p, ray_p);
            intersection.first.get_coeffs().insert(intersection.first.get_coeffs().end(), 
                    intersection.second.get_coeffs().begin(), intersection.second.get_coeffs().end());

            return py::array_t<NT>(
                py::buffer_info(
                    intersection.first.data(),
                    sizeof(NT),
                    py::format_descriptor<NT>::format(),
                    2,
                    { 2, (int) p.dimension() },
                    { sizeof(NT) * p.dimension(), sizeof(NT) }
                )
            );
        }).
        def("sample_boundary", [](Polytope &p, py::array_t<NT> internal_point,
                                   int rnum, unsigned int walk_len) -> py::array_t<NT> {
            py::buffer_info internal_point_i = internal_point.request();
            Point internal_point_p(p.dimension(), (NT* ) internal_point_i.ptr);
            std::vector<Point> rand_points;

            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            RNGType rng(seed);
            boost::normal_distribution<> rdist(0,1);
            boost::random::uniform_real_distribution<>(urdist);
            boost::random::uniform_real_distribution<> urdist1(-1,1);

            vars<NT, RNGType> var(0, p.dimension(), 0, 1, 0, 0, 0, 0.0, 0, 0, rng,
            urdist, urdist1, 0, false, false, false, false, false, false, false, true, false);

            boundary_rand_point_generator(p, internal_point_p, rnum, walk_len, rand_points, var);

            NT* return_data = new NT[rnum*p.dimension()];
            for(int i = 0; i < rnum; ++i) {
                for (int j=0; j<p.dimension(); ++j) {
                    return_data[i*p.dimension()+j] = rand_points[i][j];
                }
            }

            return py::array_t<NT>(
                py::buffer_info(
                    return_data,
                    sizeof(NT),
                    py::format_descriptor<NT>::format(),
                    2,
                    { rnum, (int) p.dimension() },
                    { sizeof(NT) * p.dimension(), sizeof(NT) }
                )
            );
        });
}
