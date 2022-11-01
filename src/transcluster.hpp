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


// inline static double logsubexpd(double x, double y)
// {
//     if (x == y)
//         return -INFINITY;

//     if (y>x)
//         return NAN;

//     return x + log1p(-exp(-(x - y)));

// }

inline std::tuple<double, double>
lprob_k_given_N(size_t N, size_t k, double delta, double lamb, double beta, const std::vector<double> &lgamma)
{

    double lprob;
    double lhs;

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
        lhs = lprob;
        lprob += integral;
    }
    else
    {
        lprob = ((N + 1) * log(lamb) + k * log(beta) + lgamma[N + k + 1] - lgamma[N + 1] - lgamma[k + 1] - (N + k + 1) * log(lamb + beta));
        lhs = lprob;
    }

    return (std::make_tuple(lprob, lhs));
}


double upper_bound_E(const std::vector<double> &lgamma, double delta, double lamb, double beta, size_t N)
{

    double diff;

    // ugly poisson cdf
    double pois_cdf = -INFINITY;
    for (size_t i = 0; i <= N; i++)
    {
        pois_cdf = logaddexpd(i * log(lamb * delta) - lgamma[i + 1], pois_cdf);
    }

    diff = exp(log(beta) + delta*lamb + log(N+1) - (log(lamb) + pois_cdf));

    return (diff);
}


double expected_k(int N, double delta, double lamb, double beta, double threshold_Ek,
                  const std::vector<double> &lgamma,
                  std::unordered_map<std::tuple<int, int, double>, std::tuple<double, double>> &kN_map)
{

    double lprob = -INFINITY;
    double elprob = -INFINITY;
    double lkN, upper_bound, diff_bound;
    int k;
    std::tuple<int, int, double> key;


    // Calculate expected value of E(K) to a given error margin.
    k=1;
    upper_bound = upper_bound_E(lgamma, delta, lamb, beta, N);
    diff_bound = threshold_Ek+1;
    while ((diff_bound > threshold_Ek) && (k<10000))
    {
        // std::cout << exp(logsubexpd(upper_bound, lprob)) << " : " << k << std::endl;
        if ((k%100)==0){
        printf("upper_bound: %f\n", upper_bound);
        printf("exp(lprob): %f\n", exp(lprob));
        printf("exp(elprob): %f\n", exp(elprob));
        printf("diff_bound: %f\n", diff_bound);
        printf("k: %i, N: %i, delta: %f\n", k, N, delta);
        }

        key = std::make_tuple(N, k, delta);

        if (!kN_map.count(key))
        {
            kN_map[key] = lprob_k_given_N(N, k, delta, lamb, beta, lgamma);
        }

        lprob = logaddexpd(lprob, std::get<0>(kN_map[key]) + log(k));
        

        // Find upper bound on E(K)
        elprob = logaddexpd(elprob, std::get<1>(kN_map[key]) + log(k) + delta*(lamb+beta) - (N+k+1)*log(lamb+beta));
        diff_bound = upper_bound - exp(elprob);

        k++;
    }

    return (exp(lprob));
}

inline std::tuple<std::vector<double>, std::vector<double>>
trans_dist(const std::vector<int> &snpdiff, const std::vector<double> &datediff, double lamb, double beta, double threshold_Ek=1e-6)
{

    // chache results
    std::unordered_map<std::tuple<int, double>, double> eK_map;
    std::unordered_map<std::tuple<int, int, double>, std::tuple<double, double>> kN_map;

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
            eK[i] = expected_k(snpdiff[i], datediff[i], lamb, beta, threshold_Ek, lg, kN_map);
            eK_map[key] = eK[i];
        }

        keyB = std::make_tuple(snpdiff[i], 0, datediff[i]);

        if (!kN_map.count(keyB))
        {
            kN_map[keyB] = lprob_k_given_N(snpdiff[i], 0, datediff[i], lamb, beta, lg);
        }
        p0[i] = std::get<0>(kN_map[keyB]);

    }

    return std::make_tuple(p0, eK);
}