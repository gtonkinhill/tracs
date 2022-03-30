#include <pybind11/stl.h>
#include "pairsnp.hpp"

namespace py = pybind11;

PYBIND11_MODULE(MTRAN, m) {
  m.doc() = "Meta Transmission Clustering";

  m.def("pairsnp", &pairsnp, py::return_value_policy::take_ownership,
        "Run pairsnp", py::arg("fasta"), py::arg("n_threads"), py::arg("dist"),
        py::arg("knn"));

}
