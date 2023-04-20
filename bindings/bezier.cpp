// pybind11_wrapper.cpp
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "../inc/Bezier.h"
#include "../inc/Types.h"

PYBIND11_MODULE(bezier, m) {
  pybind11::class_<Bezier>(m,"Bezier")
  .def(pybind11::init<int, int, scalar_t>()) // constructor
  .def("R_matrix", &Bezier::R_matrix)
  .def("S_matrix", &Bezier::S_matrix)
  .def("D_matrix", &Bezier::D_matrix)
  .def("H_matrix", &Bezier::H_matrix)
  .def("R_n", &Bezier::R_n)
  .def_readwrite("R", &Bezier::R)
  .def_readwrite("S", &Bezier::S)
  .def_readwrite("D", &Bezier::D)
  .def_readwrite("H", &Bezier::H)
  .def_readwrite("order", &Bezier::order)
  .def_readwrite("gamma", &Bezier::gamma)
  .def_readwrite("T", &Bezier::T)
 ;

 m.def("add",&add);
}


