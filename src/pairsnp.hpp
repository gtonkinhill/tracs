#pragma once

#include <chrono>
#include <iostream>
#include <stdio.h>
#include <string>
#include <zlib.h>
#include <map>

#include "kseq.h"

#include <boost/dynamic_bitset.hpp>
#include <boost/math/distributions/binomial.hpp>

template <typename T>
std::vector<T> combine_vectors(const std::vector<std::vector<T>> &vec,
                               const size_t len)
{
  std::vector<T> all(len);
  auto all_it = all.begin();
  for (size_t i = 0; i < vec.size(); ++i)
  {
    std::copy(vec[i].cbegin(), vec[i].cend(), all_it);
    all_it += vec[i].size();
  }
  return all;
}

double cached_binomial_cdf(int n, double p, int k)
{
    static std::map<std::tuple<int, double, int>, double> cache;

    auto key = std::make_tuple(n, p, k);

    auto it = cache.find(key);
    if (it != cache.end()) {
        // Return cached value if it exists
        return it->second;
    } else {
        // Calculate CDF and cache result
        boost::math::binomial_distribution<> dist(n, p);
        double cdf = boost::math::cdf(dist, k);
        cache[key] = cdf;
        return cdf;
    }
}

KSEQ_INIT(gzFile, gzread)

std::pair<size_t, size_t> load_seqs(std::string fasta,
                                    std::vector<boost::dynamic_bitset<>> &A_snps,
                                    std::vector<boost::dynamic_bitset<>> &C_snps,
                                    std::vector<boost::dynamic_bitset<>> &G_snps,
                                    std::vector<boost::dynamic_bitset<>> &T_snps,
                                    std::vector<std::string> &seq_names)
{

  int l;
  size_t seq_length = 0;
  size_t count = 0;
  gzFile fp;

  fp = gzopen(fasta.c_str(), "r");
  kseq_t *seq = kseq_init(fp);

  while (true)
  {
    l = kseq_read(seq);

    if (l == -1) // end of file
      break;
    if (l == -2)
    {
      throw std::runtime_error("Error reading FASTA!");
    }
    if (l == -3)
    {
      throw std::runtime_error("Error reading FASTA!");
    }

    // check sequence length
    if ((count > 0) && (seq->seq.l != seq_length))
    {
      throw std::runtime_error(
          "Error reading FASTA, variable sequence lengths!");
    }
    seq_length = seq->seq.l;

    seq_names.push_back(seq->name.s);
    boost::dynamic_bitset<> As(seq_length);
    boost::dynamic_bitset<> Cs(seq_length);
    boost::dynamic_bitset<> Gs(seq_length);
    boost::dynamic_bitset<> Ts(seq_length);

    for (size_t j = 0; j < seq_length; j++)
    {

      seq->seq.s[j] = std::toupper(seq->seq.s[j]);

      switch (seq->seq.s[j])
      {
      case 'A':
        As[j] = 1;
        break;
      case 'C':
        Cs[j] = 1;
        break;
      case 'G':
        Gs[j] = 1;
        break;
      case 'T':
        Ts[j] = 1;
        break;

      // M = A or C
      case 'M':
        As[j] = 1;
        Cs[j] = 1;
        break;

      // R = A or G
      case 'R':
        As[j] = 1;
        Gs[j] = 1;
        break;

      // W = A or T
      case 'W':
        As[j] = 1;
        Ts[j] = 1;
        break;

      // S = C or G
      case 'S':
        Cs[j] = 1;
        Gs[j] = 1;
        break;

      // Y = C or T
      case 'Y':
        Cs[j] = 1;
        Ts[j] = 1;
        break;

      // K = G or T
      case 'K':
        Gs[j] = 1;
        Ts[j] = 1;
        break;

      // V = A,C or G
      case 'V':
        As[j] = 1;
        Cs[j] = 1;
        Gs[j] = 1;
        break;

      // H = A,C or T
      case 'H':
        As[j] = 1;
        Cs[j] = 1;
        Ts[j] = 1;
        break;

      // D = A,G or T
      case 'D':
        As[j] = 1;
        Gs[j] = 1;
        Ts[j] = 1;
        break;

      // B = C,G or T
      case 'B':
        Cs[j] = 1;
        Gs[j] = 1;
        Ts[j] = 1;
        break;

      // N = A,C,G or T
      default:
        As[j] = 1;
        Cs[j] = 1;
        Gs[j] = 1;
        Ts[j] = 1;
        break;
      }
    }
    A_snps.push_back(As);
    C_snps.push_back(Cs);
    G_snps.push_back(Gs);
    T_snps.push_back(Ts);

    count++;
  }
  kseq_destroy(seq);
  gzclose(fp);

  return std::make_pair(count, seq_length);
}


std::pair<size_t, size_t>  range_count(const boost::dynamic_bitset<> &bit_set, size_t start, size_t end) {
    if (start > end || end > bit_set.size()) {
        throw std::out_of_range("Invalid range specified.");
    }

    size_t count = 0;
    size_t pos = bit_set.find_first();
    size_t left = 0;
    size_t length;

    while (pos != boost::dynamic_bitset<>::npos && pos < end) {
        if (pos >= start) {
            if (left == 0) {
                left = pos;
            }
            count++;
            length = pos-left+1;
        }
        pos = bit_set.find_next(pos);
    }

    return std::make_pair(length, count);
}


size_t filter_recomb(boost::dynamic_bitset<> res) {
  // Invert to look at SNPs
  res.flip();

  // Determine window size
  double d = static_cast<double>(res.count());
  int aln_length = res.size();

  double p = d / aln_length;
  double p_value_threshold = 0.05 / d;

  // Determine window size
  int window_size_half = static_cast<int>(1.0 / p / 2.0 + 1);
  window_size_half = std::min(window_size_half, 5000);
  window_size_half = std::max(window_size_half, 50);

  // std::cout << "window_size_half: " << window_size_half << std::endl;

  // Iterate through SNPs
  size_t window_count;
  size_t left, right;
  double p_value;
  size_t filtered_d = 0;

  for (int i = res.find_first(); i != boost::dynamic_bitset<>::npos; i = res.find_next(i)) {
    // std::cout << i << std::endl;

    left = std::max(0, i - window_size_half);
    right = std::min(aln_length, i + window_size_half + 1);

     
    std::pair<size_t, size_t> window_count;
    window_count = range_count(res, left, right);

    // std::cout << "window_count: " << window_count.second << std::endl;


    if (window_count.second > 1){
      p_value = 1.0 - cached_binomial_cdf(window_count.first, p, window_count.second);
      // p_value = 1.0 - cached_binomial_cdf(right-left, p, window_count.second);

      if (p_value >= p_value_threshold) {
          filtered_d++;
          // std::cout << "p_value: " << p_value << std::endl;
          // std::cout << "kept!" << std::endl;
      }
    } else {
      filtered_d++;
    }

  }

  // std::cout << "filtered_d: " << filtered_d << std::endl;
  return filtered_d;
}

inline std::tuple<std::vector<uint64_t>, std::vector<uint64_t>,
                  std::vector<uint64_t>, std::vector<std::string>, 
                  std::vector<uint64_t>>
pairsnp(const std::vector<std::string> fastas, int n_threads, int dist, bool filter)
{

  // open filename and initialise kseq
  size_t n_seqs = 0;
  size_t seq_length;
  size_t i_end, j_start;

  // initialise bitmaps
  std::vector<std::string> seq_names;
  std::vector<boost::dynamic_bitset<>> A_snps;
  std::vector<boost::dynamic_bitset<>> C_snps;
  std::vector<boost::dynamic_bitset<>> G_snps;
  std::vector<boost::dynamic_bitset<>> T_snps;

  if ((fastas.size() < 0) || (fastas.size() > 2))
  {
    throw std::runtime_error("Invalid number of fasta files!");
  }

  // Load first sequence file
  std::pair<size_t, size_t> ls;
  ls = load_seqs(fastas[0], A_snps, C_snps, G_snps, T_snps, seq_names);
  n_seqs = i_end = ls.first;
  seq_length = ls.second;

  // deal with search range and load second file if provided
  if (fastas.size() == 1)
  {
    j_start = 0;
  }
  else
  {
    j_start = n_seqs;
    n_seqs += load_seqs(fastas[1], A_snps, C_snps, G_snps, T_snps, seq_names).first;
  }

  // Set up progress meter
  // static const uint64_t n_progress_ticks = 1000;
  // uint64_t update_every = 1;
  // if (n_seqs > n_progress_ticks) {
  //   update_every = n_seqs / n_progress_ticks;
  // }
  // ProgressMeter dist_progress(n_progress_ticks, true);
  // int progress = 0;

  // Shared variables for openmp loop
  std::vector<std::vector<uint64_t>> rows(n_seqs);
  std::vector<std::vector<uint64_t>> cols(n_seqs);
  std::vector<std::vector<uint64_t>> distances(n_seqs);
  std::vector<std::vector<uint64_t>> filt_distances(n_seqs);
  uint64_t len = 0;
  bool interrupt = false;

#pragma omp parallel for schedule(static) reduction(+ \
                                                    : len) num_threads(n_threads)
  for (uint64_t i = 0; i < i_end; i++)
  {
    // Cannot throw in an openmp block, short circuit instead
    if (interrupt || PyErr_CheckSignals() != 0)
    {
      interrupt = true;
    }
    else
    {
      int d;
      boost::dynamic_bitset<> res(seq_length);

      for (uint64_t j = std::max(j_start, i + 1); j < n_seqs; j++)
      {

        res = A_snps[i] & A_snps[j];
        res |= C_snps[i] & C_snps[j];
        res |= G_snps[i] & G_snps[j];
        res |= T_snps[i] & T_snps[j];

        d = seq_length - res.count();

        if (d <= dist)
        {
          rows[i].push_back(i);
          cols[i].push_back(j);
          distances[i].push_back(d);

          if (filter){
            uint64_t fd = filter_recomb(res);
            filt_distances[i].push_back(fd);
          }
        }
      }

      len += distances[i].size();

      //       if (i % update_every == 0) {
      // #pragma omp critical
      //         dist_progress.tick_count(++progress);
      //       }
    }
  }

  // // Finalise
  // if (interrupt) {
  //   check_interrupts();
  // } else {
  //   dist_progress.finalise();
  // }

  // Combine the lists from each thread
  std::vector<uint64_t> distances_all = combine_vectors(distances, len);
  std::vector<uint64_t> filt_distances_all = combine_vectors(filt_distances, len);
  std::vector<uint64_t> rows_all = combine_vectors(rows, len);
  std::vector<uint64_t> cols_all = combine_vectors(cols, len);

  return std::make_tuple(rows_all, cols_all, distances_all, seq_names, filt_distances_all);
}
