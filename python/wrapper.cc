#include <pybind11/pybind11.h>
#include <iostream>
#include "Bend.h"
#include "Sext.h"
#include "Quad.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

PYBIND11_MODULE(pyGBCuda, m)
{
    py::class_<Bend>(m, "Bend")
        .def(py::init<const std::string &,
                      double const,
                      double const,
                      double const,
                      double const,
                      double const,
                      double const,
                      double const,
                      double const,
                      double const,
                      double const,
                      double const,
                      double const>())
        .def("Update", &Bend::Update)
        .def("SetSympass", &Bend::SetSympass)
        .def("TransMatrix", &Bend::TransMatrix);

    py::class_<Sext>(m, "Sext")
        .def(py::init<std::string const &,
                      double const,
                      double const,
                      int const,
                      double const,
                      double const,
                      double const>())
        .def("Update", &Sext::Update)
        .def("SetSympass", &Sext::SetSympass);

    py::class_<Quad>(m, "Quad")
        .def(py::init<std::string const &,
                      double const,
                      double const,
                      int const,
                      double const,
                      double const,
                      double const>())
        .def("run_qsympass4", &Quad::run_qsympass4)
        .def("Update", &Quad::Update)
        .def("SetSympass", &Quad::SetSympass)
        .def("PrintTM", &Quad::PrintTM)
        .def("TransMatrix", &Quad::TransMatrix)
        .def("DoDxDyTilt", &Quad::DoDxDyTilt)
        .def("Name", &Quad::Name);
}
