#include <pybind11/stl.h>
#include "pairsnp.hpp"
#include "transcluster.hpp"
#include "dmultinomial.hpp"

namespace py = pybind11;

PYBIND11_MODULE(TRACM, m)
{
  m.doc() = "Meta Transmission Clustering";

  m.def("pairsnp", &pairsnp, py::return_value_policy::take_ownership,
        "Run pairsnp", py::arg("fasta"), py::arg("n_threads"), py::arg("dist"));

  m.def("lprob_k_given_N", &lprob_k_given_N, py::return_value_policy::take_ownership,
        "Probability of K intermediate hosts given N SNPs",
        py::arg("N"), py::arg("k"), py::arg("delta"), py::arg("lamb"), py::arg("beta"), py::arg("lgamma"));

  m.def("trans_dist", &trans_dist, py::return_value_policy::take_ownership,
        "Calculate transmission probabilities and expectation for a vector of distances",
        py::arg("snpdiff"), py::arg("datediff"), py::arg("lamb"), py::arg("beta"), py::arg("threshold_Ek"));

  m.def("calculate_posteriors", &calculate_posteriors, py::return_value_policy::take_ownership,
        "Calculate posterior count estimates",
        py::arg("counts"), py::arg("alphas"), py::arg("keep"), py::arg("threshold"));
}
