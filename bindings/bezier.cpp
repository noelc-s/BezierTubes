// pybind11_wrapper.cpp
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "../inc/Bezier.h"
#include "../inc/Types.h"

namespace py = pybind11;

PYBIND11_MODULE(bezier, m) {
  py::class_<Bezier>(m,"Bezier")
  .def(py::init<int, int, scalar_t>()) // constructor
  .def("R_matrix", &Bezier::R_matrix)
  .def("S_matrix", &Bezier::S_matrix)
  .def("D_matrix", &Bezier::D_matrix)
  .def("H_matrix", &Bezier::H_matrix)
  .def("R_n", &Bezier::R_n)
  .def("T", py::overload_cast<scalar_t>(&Bezier::T))
  .def("T", py::overload_cast<vector_t>(&Bezier::T))
  .def("b", py::overload_cast<scalar_t, matrix_t>(&Bezier::b))
  .def("b", py::overload_cast<vector_t, matrix_t>(&Bezier::b))
  .def("db", py::overload_cast<scalar_t,int, matrix_t>(&Bezier::db))
  .def("db", py::overload_cast<vector_t,int, matrix_t>(&Bezier::db))
  .def_readwrite("R", &Bezier::R)
  .def_readwrite("S", &Bezier::S)
  .def_readwrite("D", &Bezier::D)
  .def_readwrite("H", &Bezier::H)
  .def_readwrite("order", &Bezier::order)
  .def_readwrite("gamma", &Bezier::gamma)
  .def_readwrite("tau", &Bezier::tau)
 ;
}


