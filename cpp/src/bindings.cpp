#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "placeopt.hpp"

namespace py = pybind11;
using namespace placeopt;

PYBIND11_MODULE(_core, m) {
    m.doc() = "PlaceOpt C++ core: simulated-annealing and analytical standard-cell placement";
    py::class_<Pin>(m, "Pin")
        .def(py::init([](int owner, double dx, double dy) { return Pin{owner, dx, dy}; }),
             py::arg("owner"), py::arg("dx") = 0.0, py::arg("dy") = 0.0)
        .def_readwrite("owner", &Pin::owner).def_readwrite("dx", &Pin::dx).def_readwrite("dy", &Pin::dy);
    py::class_<Netlist>(m, "Netlist")
        .def(py::init<>())
        .def_readwrite("rows", &Netlist::rows).def_readwrite("sites", &Netlist::sites)
        .def_readwrite("site_w", &Netlist::site_w).def_readwrite("row_h", &Netlist::row_h)
        .def_readwrite("width", &Netlist::width).def_readwrite("term_x", &Netlist::term_x)
        .def_readwrite("term_y", &Netlist::term_y).def_readwrite("nets", &Netlist::nets)
        .def("build", &Netlist::build).def("num_cells", &Netlist::num_cells);
    py::class_<Placement>(m, "Placement")
        .def(py::init<>())
        .def(py::init([](std::vector<int> r, std::vector<int> s) { return Placement{r, s}; }))
        .def_readwrite("row", &Placement::row).def_readwrite("site", &Placement::site);
    py::class_<AnnealParams>(m, "AnnealParams")
        .def(py::init<>())
        .def_readwrite("seed", &AnnealParams::seed).def_readwrite("moves_per_cell", &AnnealParams::moves_per_cell)
        .def_readwrite("init_temp", &AnnealParams::init_temp).def_readwrite("exit_ratio", &AnnealParams::exit_ratio)
        .def_readwrite("max_moves", &AnnealParams::max_moves)
        .def_readwrite("start_window_full", &AnnealParams::start_window_full)
        .def_readwrite("init_window", &AnnealParams::init_window).def_readwrite("snapshots", &AnnealParams::snapshots);
    py::class_<TempStep>(m, "TempStep")
        .def_readonly("temp", &TempStep::temp).def_readonly("cost", &TempStep::cost)
        .def_readonly("accept_rate", &TempStep::accept_rate).def_readonly("window", &TempStep::window)
        .def_readonly("moves", &TempStep::moves);
    py::class_<AnnealResult>(m, "AnnealResult")
        .def_readonly("placement", &AnnealResult::placement).def_readonly("hpwl", &AnnealResult::hpwl)
        .def_readonly("moves", &AnnealResult::moves).def_readonly("seconds", &AnnealResult::seconds)
        .def_readonly("trace", &AnnealResult::trace).def_readonly("snapshots", &AnnealResult::snapshots);
    py::class_<GlobalResult>(m, "GlobalResult")
        .def_readonly("x", &GlobalResult::x).def_readonly("y", &GlobalResult::y)
        .def_readonly("overflow", &GlobalResult::overflow).def_readonly("seconds", &GlobalResult::seconds);
    m.def("total_hpwl", &total_hpwl);
    m.def("is_legal", [](const Netlist& nl, const Placement& p) { std::string why; bool ok = is_legal(nl, p, &why); return py::make_tuple(ok, why); });
    m.def("random_placement", &random_placement);
    m.def("anneal", &anneal, py::call_guard<py::gil_scoped_release>());
    m.def("global_place", &global_place, py::arg("nl"), py::arg("iterations") = 30, py::arg("target_density") = 1.0,
          py::call_guard<py::gil_scoped_release>());
    m.def("legalize", &legalize);
}
