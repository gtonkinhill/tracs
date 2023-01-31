#include <pybind11/numpy.h>
#include <math.h>
#include <numeric>   // std::iota
#include <algorithm> // std::sort, std::stable_sort

namespace py = pybind11;

py::array_t<double> calculate_posteriors(py::array_t<double> counts, std::vector<double> alphas,
                                         bool keep, double expected)
{

    // total alpha
    std::sort(alphas.begin(), alphas.end(), std::greater<double>());
    double a0 = std::accumulate(alphas.begin(), alphas.end(), 0.0);
    double a_min = alphas[0] / a0;

    py::buffer_info buf1 = counts.request();

    // Apply resources
    py::array_t<double> posterior = py::array_t<double>(buf1.size);
    // Resize to 2d array
    posterior.resize({buf1.shape[0], buf1.shape[1]});

    py::buffer_info buf_posterior = posterior.request();

    // Pointer reads and writes numpy.ndarray
    double *ptr1 = (double *)buf1.ptr;
    double *ptr_result = (double *)buf_posterior.ptr;

    std::vector<size_t> idx(alphas.size());
    std::vector<double> row(buf1.shape[1]);

    for (int i = 0; i < buf1.shape[0]; i++)
    {
        double denom = 0;
        std::iota(idx.begin(), idx.end(), 0);

        for (int j = 0; j < buf1.shape[1]; j++)
        {
            denom += ptr1[i * buf1.shape[1] + j];
            row[j] = ptr1[i * buf1.shape[1] + j];
        }

        // argsort counts
        std::stable_sort(idx.begin(), idx.end(),
                         [&row](size_t i1, size_t i2)
                         { return row[i1] > row[i2]; });

        size_t alpha_index = 0;

        for (int j = 0; j < buf1.shape[1]; j++)
        {
            if (denom <= 0)
            {
                ptr_result[i * buf1.shape[1] + j] = a_min;
            }
            else
            {
                ptr_result[i * buf1.shape[1] + idx[j]] = (row[idx[j]] + alphas[alpha_index]) / (denom + a0);

                if ((j < (buf1.shape[1] - 1)) && (row[idx[j]] != row[idx[j + 1]]))
                {
                    alpha_index += 1;
                }
            }
        }

        // keep alleles even if the posterior is below the threshold else set to 0.
        for (size_t j = 0; j < buf1.shape[1]; j++)
        {
            if ((ptr_result[i * buf1.shape[1] + j] <= expected))
            {
                if (keep && (row[j] > 0))
                {
                    ptr_result[i * buf1.shape[1] + j] = expected;
                }
                else
                {
                    ptr_result[i * buf1.shape[1] + j] = 0.0;
                }
            }
        }
    }

    return posterior;
}