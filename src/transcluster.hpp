// #include <pybind11/numpy.h>
#include <math.h>
#include <map>
#include <tuple>

namespace std
{
    namespace
    {

        // Code from boost
        // Reciprocal of the golden ratio helps spread entropy
        //     and handles duplicates.
        // See Mike Seymour in magic-numbers-in-boosthash-combine:
        //     http://stackoverflow.com/questions/4948780

        template <class T>
        inline void hash_combine(std::size_t &seed, T const &v)
        {
            seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }

        // Recursive template code derived from Matthieu M.
        template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
        struct HashValueImpl
        {
            static void apply(size_t &seed, Tuple const &tuple)
            {
                HashValueImpl<Tuple, Index - 1>::apply(seed, tuple);
                hash_combine(seed, std::get<Index>(tuple));
            }
        };

        template <class Tuple>
        struct HashValueImpl<Tuple, 0>
        {
            static void apply(size_t &seed, Tuple const &tuple)
            {
                hash_combine(seed, std::get<0>(tuple));
            }
        };
    }

    template <typename... TT>
    struct hash<std::tuple<TT...>>
    {
        size_t
        operator()(std::tuple<TT...> const &tt) const
        {
            size_t seed = 0;
            HashValueImpl<std::tuple<TT...>>::apply(seed, tt);
            return seed;
        }
    };
}

/* Computes ㏒ₑ(𝑒ˣ + 𝑒ʸ) in a safe and accurate way.
 *
 * For example, `log(exp(1e3) + exp(-INFINITY))` will likely overflow,
 * while `logaddexp(1e3, -INFINITY)` will return `1e3`.
 */
inline static double logaddexpd(double x, double y)
{
    double const tmp = x - y;

    if (x == y)
        return x + M_LN2;

    if (tmp > 0)
        return x + log1p(exp(-tmp));
    else if (tmp <= 0)
        return y + log1p(exp(tmp));

    return tmp;
}

double lprob_k_given_N(size_t N, size_t k, double delta, double lamb, double beta, const std::vector<double> &lgamma)
{

    double lprob;

    if (delta > 0)
    {
        lprob = ((N + 1) * log(lamb) - delta * (lamb + beta) + k * log(beta) - lgamma[k + 1]);

        // ugly poisson cdf
        double pois_cdf = -INFINITY;
        for (size_t i = 0; i <= N; i++)
        {
            pois_cdf = logaddexpd(i * log(lamb * delta) - lgamma[i + 1], pois_cdf);
        }
        pois_cdf -= lamb * delta;
        lprob -= pois_cdf;

        double integral = -INFINITY;
        for (size_t i = 0; i <= N + k; i++)
        {
            integral = logaddexpd(
                lgamma[N + k + 1] - lgamma[i + 1] - lgamma[N + k - i + 1] + (N + k - i) * log(delta) + lgamma[i + 1] - (i + 1) * log(lamb + beta),
                integral);
        }

        integral -= lgamma[N + 1];
        lprob += integral;
    }
    else
    {
        lprob = ((N + 1) * log(lamb) + k * log(beta) + lgamma[N + k + 1] - lgamma[N + 1] - lgamma[k + 1] - (N + k + 1) * log(lamb + beta));
    }

    return (lprob);
}

double expected_k(int N, double delta, double lamb, double beta, int max_k,
                  const std::vector<double> &lgamma,
                  std::unordered_map<std::tuple<int, int, double>, double> &kN_map)
{

    double lprob = -INFINITY;
    double lkN;
    std::tuple<int, int, double> key;

    for (int k = 1; k <= max_k; k++)
    {
        key = std::make_tuple(N, k, delta);

        if (kN_map.count(key))
        {
            lkN = kN_map[key];
        }
        else
        {
            lkN = lprob_k_given_N(N, k, delta, lamb, beta, lgamma);
            kN_map[key] = lkN;
        }

        lprob = logaddexpd(lprob, lkN + log(k));
    }

    return (exp(lprob));
}

inline std::tuple<std::vector<double>, std::vector<double>>
trans_dist(const std::vector<int> &snpdiff, const std::vector<double> &datediff, double lamb, double beta)
{

    // chache results
    std::unordered_map<std::tuple<int, double>, double> eK_map;
    std::unordered_map<std::tuple<int, int, double>, double> kN_map;

    // results vectors
    std::vector<double> eK(snpdiff.size());
    std::vector<double> p0(snpdiff.size());

    // precalculate lgamma
    std::vector<double> lg;
    lg.reserve(10000);
    for (double i = 0; i < 10000; i++)
    {
        lg.push_back(std::lgamma(i));
    }

    std::tuple<int, double> key;
    std::tuple<int, int, double> keyB;

    for (size_t i = 0; i < snpdiff.size(); i++)
    {
        key = std::make_tuple(snpdiff[i], datediff[i]);
        if (eK_map.count(key))
        {
            eK[i] = eK_map[key];
        }
        else
        {
            eK[i] = expected_k(snpdiff[i], datediff[i], lamb, beta, 100, lg, kN_map);
            eK_map[key] = eK[i];
        }

        keyB = std::make_tuple(snpdiff[i], 0, datediff[i]);
        if (kN_map.count(keyB))
        {
            p0[i] = kN_map[keyB];
        }
        else
        {
            p0[i] = lprob_k_given_N(snpdiff[i], 0, datediff[i], lamb, beta, lg);
            kN_map[keyB] = p0[i];
        }
    }

    return std::make_tuple(p0, eK);
}