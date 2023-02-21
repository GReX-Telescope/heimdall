// -*- c++ -*-
/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include <iostream>
#include <memory>
#include <vector>

#include <algorithm>
#include <fstream>
#include <functional>
#include <iomanip>
#include <sstream>
#include <string>
#include <time.h>
#include <utility> // For std::pair

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include <thrust/copy.h>
#include <thrust/extrema.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/gather.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/version.h>

#include "hd/clean_filterbank_rfi.hpp"
#include "hd/maths.hpp"
#include "hd/pipeline.hpp"
#include "hd/strided_range.hpp"

#include "hd/find_giants.hpp"
#include "hd/get_rms.hpp"
#include "hd/label_candidate_clusters.hpp"
#include "hd/matched_filter.hpp"
#include "hd/merge_candidates.hpp"
#include "hd/remove_baseline.hpp"

#include "hd/DataSource.hpp"
#include "hd/stopwatch.hpp" // For benchmarking

#include <dedisp.h>

// Socket stuff
#include <arpa/inet.h>  // htons, inet_addr
#include <netinet/in.h> // sockaddr_in
#include <sys/socket.h> // socket, sendto
#include <sys/types.h>  // uint16_t
#include <unistd.h>     // close

struct hd_pipeline_t {
  hd_params params;
  dedisp_plan dedispersion_plan;
  // Memory buffers used during pipeline execution
  std::vector<hd_byte> h_clean_filterbank;
  thrust::host_vector<hd_byte> h_dm_series;
  thrust::device_vector<hd_float> d_time_series;
  thrust::device_vector<hd_float> d_filtered_series;
  // Exfil options
  int socket;
  sockaddr_in dest_addr;
  std::ofstream file;
};

void write_candidates(hd_pipeline pl, hd_size nsamps, hd_size first_idx,
                      hd_size nsamps_computed, const float *dm_list,
                      thrust::host_vector<hd_float> giant_peaks,
                      thrust::host_vector<hd_size> giant_inds,
                      thrust::host_vector<hd_size> giant_filter_inds,
                      thrust::host_vector<hd_size> giant_dm_inds) {

  hd_size filterbank_ind;
  hd_size overlap =
      pl->params.boxcar_max + dedisp_get_max_delay(pl->dedispersion_plan);
  hd_size block_size = nsamps - overlap;

  std::stringstream ss("");
  // Set some formatting parameters
  ss << std::setw(2) << std::setfill('0');

  if (first_idx > 0) {
    for (hd_size i = 0; i < giant_peaks.size(); ++i) {
      if (giant_peaks[i] > pl->params.detect_thresh) {
        hd_size giant_index = giant_inds[i] % nsamps;
        hd_size beam_no = giant_inds[i] / nsamps + pl->params.beam;
        hd_size samp_idx = first_idx + giant_index;
        hd_size block_no = (giant_index + first_idx) /
                           (nsamps - pl->params.boxcar_max -
                            dedisp_get_max_delay(pl->dedispersion_plan));
        if (giant_index < overlap)
          filterbank_ind = block_no * block_size * pl->params.nbeams +
                           (beam_no + 1) * block_size + giant_index - overlap;
        else
          filterbank_ind = block_no * block_size * pl->params.nbeams +
                           (beam_no - 1) * block_size + giant_index + nsamps -
                           2 * overlap;

        if (giant_index < nsamps_computed + pl->params.boxcar_max / 2) {
          // Serialize to the stringstream
          ss << giant_peaks[i] << "\t" << filterbank_ind << "\t" << samp_idx
             << "\t" << samp_idx * pl->params.dt << "\t" << giant_filter_inds[i]
             << "\t" << giant_dm_inds[i] << "\t" << dm_list[giant_dm_inds[i]]
             << "\t" << beam_no << std::endl;
        }
      }
    }
    // If we have a coincidencer, write output to that,
    // instead of to a file
    if (pl->params.coincidencer_host != NULL &&
        pl->params.coincidencer_port != -1) {
      int n_bytes = ::sendto(pl->socket, ss.str().c_str(), ss.str().length(), 0,
                             reinterpret_cast<sockaddr *>(&pl->dest_addr),
                             sizeof(pl->dest_addr));
      std::cout << "Wrote " << n_bytes << " bytes to the socket" << std::endl;
    } else {
      pl->file << ss.str();
    }
  }
}

void tfunc(std::vector<hd_byte> &vec) {
  hd_byte *ddata = (hd_byte *)malloc(sizeof(hd_byte) * 200);
  std::copy(vec.begin(), vec.end(), ddata);
  std::copy(vec.begin(), vec.end(), ddata + 100);

  for (int i = 0; i < 200; i++)
    std::cout << +ddata[i] << " ";
  std::cout << " " << std::endl;

  free(ddata);
}

#define HD_BENCHMARK

#ifdef HD_BENCHMARK
void start_timer(Stopwatch &timer) { timer.start(); }
void stop_timer(Stopwatch &timer) {
  cudaDeviceSynchronize();
  timer.stop();
}
#else
void start_timer(Stopwatch &timer) {}
void stop_timer(Stopwatch &timer) {}
#endif // HD_BENCHMARK

template <typename T, typename U> std::pair<T &, U &> tie(T &a, U &b) {
  return std::pair<T &, U &>(a, b);
}

hd_error allocate_gpu(const hd_pipeline pl) {
  // This is just a simple proc-->GPU heuristic to get us started
  int gpu_count;
  cudaGetDeviceCount(&gpu_count);
  int proc_idx = 0;
  int gpu_idx = pl->params.gpu_id;

  cudaError_t cerror = cudaSetDevice(gpu_idx);
  if (cerror != cudaSuccess) {
    std::cerr << "Could not setCudaDevice to " << gpu_idx << ": "
              << cudaGetErrorString(cerror) << std::endl;
    return throw_cuda_error(cerror);
  }

  if (pl->params.verbosity >= 1) {
    std::cout << "Process " << proc_idx << " using GPU " << gpu_idx
              << std::endl;
  }

  if (!pl->params.yield_cpu) {
    if (pl->params.verbosity >= 2) {
      std::cout << "\tProcess " << proc_idx << " setting CPU to spin"
                << std::endl;
    }
    cerror = cudaSetDeviceFlags(cudaDeviceScheduleSpin);
    if (cerror != cudaSuccess) {
      return throw_cuda_error(cerror);
    }
  } else {
    if (pl->params.verbosity >= 2) {
      std::cout << "\tProcess " << proc_idx << " setting CPU to yield"
                << std::endl;
    }
    // Note: This Yield flag doesn't seem to work properly.
    //   The BlockingSync flag does the job, although it may interfere
    //     with GPU/CPU overlapping (not currently used).
    // cerror = cudaSetDeviceFlags(cudaDeviceScheduleYield);
    cerror = cudaSetDeviceFlags(cudaDeviceBlockingSync);
    if (cerror != cudaSuccess) {
      return throw_cuda_error(cerror);
    }
  }

  return HD_NO_ERROR;
}

unsigned int get_filter_index(unsigned int filter_width) {
  // Uses the compiler builtin to count leading zeros to get the log2 of the
  // 32bit power of two filter_width. This is twice as fast as the LUT solution.
  return sizeof(unsigned int) * 8 - __builtin_clz(filter_width) - 1;
}

hd_error hd_create_pipeline(hd_pipeline *pipeline_, hd_params params) {
  *pipeline_ = 0;

  // Note: We use a smart pointer here to automatically clean up after errors
  typedef std::unique_ptr<hd_pipeline_t> smart_pipeline_ptr;
  smart_pipeline_ptr pipeline = smart_pipeline_ptr(new hd_pipeline_t());
  if (!pipeline.get()) {
    return throw_error(HD_MEM_ALLOC_FAILED);
  }

  // Setup the writing context
  if (params.coincidencer_host != NULL && params.coincidencer_port != -1) {
    // Create the socket
    int sock = ::socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    sockaddr_in dest_addr;
    dest_addr.sin_family = AF_INET;
    dest_addr.sin_port = htons(params.coincidencer_port);
    dest_addr.sin_addr.s_addr = inet_addr(params.coincidencer_host);
    pipeline->socket = sock;
    pipeline->dest_addr = dest_addr;
  } else {
    // Create the file
    std::string filename = std::string(params.output_dir) + "/giants.cand";
    pipeline->file = std::ofstream(filename.c_str(), std::ios::app);
  }

  pipeline->params = params;

  if (params.verbosity >= 2) {
    std::cout << "\tAllocating GPU..." << std::endl;
  }

  hd_error error = allocate_gpu(pipeline.get());
  if (error != HD_NO_ERROR) {
    return throw_error(error);
  }

  if (params.verbosity >= 1) {
    std::cout << "nchans = " << params.nchans << std::endl;
    std::cout << "dt     = " << params.dt << std::endl;
    std::cout << "f0     = " << params.f0 << std::endl;
    std::cout << "df     = " << params.df << std::endl;
    std::cout << "nsnap     = " << params.nsnap << std::endl;
  }

  if (params.verbosity >= 2) {
    std::cout << "\tCreating dedispersion plan..." << std::endl;
  }

  dedisp_error derror;
  derror = dedisp_create_plan(&pipeline->dedispersion_plan, params.nchans,
                              params.dt, params.f0, params.df);
  if (derror != DEDISP_NO_ERROR) {
    return throw_dedisp_error(derror);
  }
  // TODO: Consider loading a pre-generated DM list instead for flexibility
  derror = dedisp_generate_dm_list(
      pipeline->dedispersion_plan, pipeline->params.dm_min,
      pipeline->params.dm_max, pipeline->params.dm_pulse_width,
      pipeline->params.dm_tol);
  if (derror != DEDISP_NO_ERROR) {
    return throw_dedisp_error(derror);
  }

  if (pipeline->params.use_scrunching) {
    derror = dedisp_enable_adaptive_dt(pipeline->dedispersion_plan,
                                       pipeline->params.dm_pulse_width,
                                       pipeline->params.scrunch_tol);
    if (derror != DEDISP_NO_ERROR) {
      return throw_dedisp_error(derror);
    }
  }

  *pipeline_ = pipeline.release();

  if (params.verbosity >= 2) {
    std::cout << "\tInitialisation complete." << std::endl;
  }

  if (params.verbosity >= 1) {
    std::cout << "Using Thrust v" << THRUST_MAJOR_VERSION << "."
              << THRUST_MINOR_VERSION << "." << THRUST_SUBMINOR_VERSION
              << std::endl;
  }

  return HD_NO_ERROR;
}

hd_error hd_execute(hd_pipeline pl, const hd_byte *h_filterbank, hd_size nsamps,
                    hd_size nbits, hd_size first_idx, hd_size iidx,
                    hd_size *nsamps_processed, hd_size gulp_idx) {

  hd_error error = HD_NO_ERROR;

  Stopwatch total_timer;
  Stopwatch memory_timer;
  Stopwatch clean_timer;
  Stopwatch dedisp_timer;
  Stopwatch communicate_timer;
  Stopwatch copy_timer;
  Stopwatch baseline_timer;
  Stopwatch normalise_timer;
  Stopwatch filter_timer;
  Stopwatch coinc_timer;
  Stopwatch giants_timer;
  Stopwatch candidates_timer;

  start_timer(total_timer);

  start_timer(clean_timer);
  // Note: Filterbank cleaning must be done out-of-place
  hd_size nbytes = (nsamps)*pl->params.nchans * nbits / 8 * pl->params.nbeams;
  start_timer(memory_timer);
  pl->h_clean_filterbank.resize(nbytes);
  // std::vector<int>          h_killmask(pl->params.nchans, 1);
  stop_timer(memory_timer);

  // copy to clean filterbank
  std::copy(h_filterbank, h_filterbank + nbytes,
            pl->h_clean_filterbank.begin());

  // apply manual killmasks
  /*  error = apply_manual_killmasks (pl->dedispersion_plan,
                                  &h_killmask[0],
                                  pl->params.num_channel_zaps,
                                  pl->params.channel_zaps);
  if( error != HD_NO_ERROR ) {
    return throw_error(error);
  }

  hd_size good_chan_count = thrust::reduce(h_killmask.begin(),
                                           h_killmask.end());
  hd_size bad_chan_count = pl->params.nchans - good_chan_count;
  if( pl->params.verbosity >= 2 ) {
   std::cout << "Bad channel count = " << bad_chan_count << std::endl;
  }
  */
  stop_timer(clean_timer);

  if (pl->params.verbosity >= 2) {
    std::cout << "\tGenerating DM list..." << std::endl;
  }

  if (pl->params.verbosity >= 3) {
    std::cout << "dm_min = " << pl->params.dm_min << std::endl;
    std::cout << "dm_max = " << pl->params.dm_max << std::endl;
    std::cout << "dm_tol = " << pl->params.dm_tol << std::endl;
    std::cout << "dm_pulse_width = " << pl->params.dm_pulse_width << std::endl;
    std::cout << "nchans = " << pl->params.nchans << std::endl;
    std::cout << "dt = " << pl->params.dt << std::endl;

    std::cout << "dedisp nchans = "
              << dedisp_get_channel_count(pl->dedispersion_plan) << std::endl;
    std::cout << "dedisp dt = " << dedisp_get_dt(pl->dedispersion_plan)
              << std::endl;
    std::cout << "dedisp f0 = " << dedisp_get_f0(pl->dedispersion_plan)
              << std::endl;
    std::cout << "dedisp df = " << dedisp_get_df(pl->dedispersion_plan)
              << std::endl;
  }

  hd_size dm_count = dedisp_get_dm_count(pl->dedispersion_plan);
  const float *dm_list = dedisp_get_dm_list(pl->dedispersion_plan);

  const dedisp_size *scrunch_factors =
      dedisp_get_dt_factors(pl->dedispersion_plan);
  if (pl->params.verbosity >= 3) {
    std::cout << "DM List for " << pl->params.dm_min << " to "
              << pl->params.dm_max << std::endl;
    for (hd_size i = 0; i < dm_count; ++i) {
      std::cout << dm_list[i] << std::endl;
    }
  }

  if (pl->params.verbosity >= 2) {
    std::cout << "Scrunch factors:" << std::endl;
    for (hd_size i = 0; i < dm_count; ++i) {
      std::cout << scrunch_factors[i] << " ";
    }
    std::cout << std::endl;
  }

  // Set channel killmask for dedispersion
  // dedisp_set_killmask(pl->dedispersion_plan, &h_killmask[0]);
  hd_size nsamps_computed =
      nsamps - dedisp_get_max_delay(pl->dedispersion_plan);
  hd_size series_stride = nsamps_computed;

  // Report the number of samples that will be properly processed
  *nsamps_processed = nsamps_computed - pl->params.boxcar_max;

  if (pl->params.verbosity >= 3) {
    std::cout << "dm_count = " << dm_count << std::endl;
    std::cout << "max delay = " << dedisp_get_max_delay(pl->dedispersion_plan)
              << std::endl;
    std::cout << "nsamps_computed = " << nsamps_computed << std::endl;
  }

  start_timer(memory_timer);

  pl->h_dm_series.resize((nsamps * (pl->params.nbeams - 1) + series_stride) *
                         pl->params.dm_nbits / 8 * dm_count);
  pl->d_time_series.resize(series_stride + (pl->params.nbeams - 1) * nsamps);
  pl->d_filtered_series.resize(series_stride + (pl->params.nbeams - 1) * nsamps,
                               0);

  stop_timer(memory_timer);

  GetRMSPlan rms_getter;
  RemoveBaselinePlan baseline_remover;
  MatchedFilterPlan<hd_float> matched_filter_plan;
  GiantFinder giant_finder;

  thrust::device_vector<hd_float> d_giant_peaks;
  thrust::device_vector<hd_size> d_giant_inds;
  thrust::device_vector<hd_size> d_giant_begins;
  thrust::device_vector<hd_size> d_giant_ends;
  thrust::device_vector<hd_size> d_giant_filter_inds;
  thrust::device_vector<hd_size> d_giant_dm_inds;
  thrust::device_vector<hd_size> d_giant_members;

  typedef thrust::device_ptr<hd_float> dev_float_ptr;
  typedef thrust::device_ptr<hd_size> dev_size_ptr;

  if (pl->params.verbosity >= 2) {
    std::cout << "\tDedispersing for DMs " << dm_list[0] << " to "
              << dm_list[dm_count - 1] << "..." << std::endl;
  }

  // Dedisperse
  dedisp_error derror;
  const dedisp_byte *in = &pl->h_clean_filterbank[0];
  dedisp_byte *out = &pl->h_dm_series[0];
  dedisp_size in_nbits = nbits;
  dedisp_size in_stride = pl->params.nchans * in_nbits / 8;
  dedisp_size out_nbits = pl->params.dm_nbits;
  dedisp_size out_stride = series_stride * out_nbits / 8 +
                           (pl->params.nbeams - 1) * nsamps * out_nbits / 8;
  unsigned flags = 0;
  start_timer(dedisp_timer);
  derror = dedisp_execute_adv(pl->dedispersion_plan, nsamps * pl->params.nbeams,
                              in, in_nbits, in_stride, out, out_nbits,
                              out_stride, flags);

  stop_timer(dedisp_timer);
  if (derror != DEDISP_NO_ERROR) {
    return throw_dedisp_error(derror);
  }

  if (pl->params.verbosity >= 2) {
    std::cout << "\tBeginning inner pipeline..." << std::endl;
  }

  bool too_many_giants = false;
  hd_size dm_idx_output = 0;

  // try to ease into pipeline
  if (gulp_idx > 1) {

    // For each DM
    for (hd_size dm_idx = 0; dm_idx < dm_count; ++dm_idx) {
      dm_idx_output = dm_idx;

      hd_size cur_dm_scrunch = scrunch_factors[dm_idx];
      hd_size cur_nsamps =
          (nsamps_computed + nsamps * (pl->params.nbeams - 1)) / cur_dm_scrunch;
      hd_float cur_dt = pl->params.dt * cur_dm_scrunch;

      // Bail if the candidate rate is too high
      if (too_many_giants) {
        break;
      }

      if (pl->params.verbosity >= 4) {
        std::cout << "dm_idx     = " << dm_idx << std::endl;
        std::cout << "scrunch    = " << scrunch_factors[dm_idx] << std::endl;
        std::cout << "cur_nsamps = " << cur_nsamps << std::endl;
        std::cout << "dt0        = " << pl->params.dt << std::endl;
        std::cout << "cur_dt     = " << cur_dt << std::endl;

        std::cout << "\tBaselining and normalising each beam..." << std::endl;
      }

      hd_float *time_series = thrust::raw_pointer_cast(&pl->d_time_series[0]);

      // Copy the time series to the device and convert to floats
      hd_size offset = dm_idx *
                       (series_stride + (pl->params.nbeams - 1) * nsamps) *
                       pl->params.dm_nbits / 8;
      start_timer(copy_timer);
      switch (pl->params.dm_nbits) {
      case 8:
        thrust::copy((unsigned char *)&pl->h_dm_series[offset],
                     (unsigned char *)&pl->h_dm_series[offset] + cur_nsamps,
                     pl->d_time_series.begin());
        break;
      case 16:
        thrust::copy((unsigned short *)&pl->h_dm_series[offset],
                     (unsigned short *)&pl->h_dm_series[offset] + cur_nsamps,
                     pl->d_time_series.begin());
        break;
      case 32:
        // Note: 32-bit implies float, not unsigned int
        thrust::copy((float *)&pl->h_dm_series[offset],
                     (float *)&pl->h_dm_series[offset] + cur_nsamps,
                     pl->d_time_series.begin());
        break;
      default:
        return HD_INVALID_NBITS;
      }

      stop_timer(copy_timer);

      // Remove the baseline
      // -------------------
      // Note: Divided by 2 to form a smoothing radius
      hd_size nsamps_smooth =
          hd_size(pl->params.baseline_length / (2 * cur_dt));
      // Crop the smoothing length in case not enough samples
      start_timer(baseline_timer);

      // TESTING
      error = baseline_remover.exec(time_series, cur_nsamps, nsamps_smooth);
      stop_timer(baseline_timer);
      if (error != HD_NO_ERROR) {
        return throw_error(error);
      }

      // Normalise
      // ---------
      start_timer(normalise_timer);
      hd_float rms = rms_getter.exec(time_series, cur_nsamps);
      thrust::transform(pl->d_time_series.begin(), pl->d_time_series.end(),
                        thrust::make_constant_iterator(hd_float(1.0) / rms),
                        pl->d_time_series.begin(),
                        thrust::multiplies<hd_float>());
      stop_timer(normalise_timer);

      // Prepare the boxcar filters
      // --------------------------
      // We can't process the first and last max-filter-width/2 samples
      hd_size rel_boxcar_max = pl->params.boxcar_max / cur_dm_scrunch;

      hd_size max_nsamps_filtered = cur_nsamps + 1 - rel_boxcar_max;
      // This is the relative offset into the time series of the filtered data
      hd_size cur_filtered_offset = rel_boxcar_max / 2;
      // minimum filter width
      hd_size min_filter_width = std::max(cur_dm_scrunch, hd_size(1));

      // Create and prepare matched filtering operations
      start_timer(filter_timer);
      // Note: Filter width is relative to the current time resolution
      matched_filter_plan.prep(time_series, cur_nsamps, rel_boxcar_max);
      stop_timer(filter_timer);
      // --------------------------

      hd_float *filtered_series =
          thrust::raw_pointer_cast(&pl->d_filtered_series[0]);

      // Note: Filtering is done using a combination of tscrunching and
      //         'proper' boxcar convolution. The parameter min_tscrunch_width
      //         indicates how much of each to do. Raising min_tscrunch_width
      //         increases sensitivity but decreases performance and vice
      //         versa.

      // For each boxcar filter
      // Note: We cannot detect pulse widths < current time resolution

      // Will make it a command line option to double or linearly increase the
      // filter width? boxcar filter loop starts

      // This variable is unused
      // int boxcar_inc = pl->params.boxcar_max / pl->params.n_boxcar_inc;

      for (hd_size filter_width = min_filter_width;
           filter_width <= pl->params.boxcar_max; filter_width *= 2) {
        hd_size rel_filter_width = filter_width / cur_dm_scrunch;
        hd_size filter_idx = filter_width;

        if (pl->params.verbosity >= 4) {
          std::cout << "Filtering each beam at width of " << filter_width
                    << std::endl;
        }

        // Note: Filter width is relative to the current time resolution
        hd_size rel_min_tscrunch_width = std::max(
            pl->params.min_tscrunch_width / cur_dm_scrunch, hd_size(1));
        hd_size rel_tscrunch_width =
            std::max(2 * rel_filter_width / rel_min_tscrunch_width, hd_size(1));
        // Filter width relative to cur_dm_scrunch AND tscrunch
        hd_size rel_rel_filter_width = rel_filter_width / rel_tscrunch_width;

        start_timer(filter_timer);

        error = matched_filter_plan.exec(filtered_series, rel_filter_width,
                                         rel_tscrunch_width);

        if (error != HD_NO_ERROR) {
          return throw_error(error);
        }
        // Divide and round up
        hd_size cur_nsamps_filtered =
            ((max_nsamps_filtered - 1) / rel_tscrunch_width + 1);
        hd_size cur_scrunch = cur_dm_scrunch * rel_tscrunch_width;

        // TESTING Proper normalisation

        hd_float rms = rms_getter.exec(filtered_series, cur_nsamps_filtered);
        thrust::transform(thrust::device_ptr<hd_float>(filtered_series),
                          thrust::device_ptr<hd_float>(filtered_series) +
                              cur_nsamps_filtered,
                          thrust::make_constant_iterator(hd_float(1.0) / rms),
                          thrust::device_ptr<hd_float>(filtered_series),
                          thrust::multiplies<hd_float>());

        hd_size prev_giant_count = d_giant_peaks.size();

        start_timer(giants_timer);

        error = giant_finder.exec(
            filtered_series, cur_nsamps_filtered, pl->params.detect_thresh,
            pl->params.cand_sep_time * rel_rel_filter_width, d_giant_peaks,
            d_giant_inds, d_giant_begins, d_giant_ends);

        if (error != HD_NO_ERROR) {
          return throw_error(error);
        }

        hd_size rel_cur_filtered_offset =
            (cur_filtered_offset / rel_tscrunch_width);

        auto transform_fn = [rel_cur_filtered_offset,
                             cur_scrunch] __device__(hd_size a) {
          return (a + rel_cur_filtered_offset) * cur_scrunch;
        };

        thrust::transform(
            d_giant_inds.begin() + prev_giant_count, d_giant_inds.end(),
            d_giant_inds.begin() + prev_giant_count, transform_fn);
        thrust::transform(
            d_giant_begins.begin() + prev_giant_count, d_giant_begins.end(),
            d_giant_begins.begin() + prev_giant_count, transform_fn);
        thrust::transform(
            d_giant_ends.begin() + prev_giant_count, d_giant_ends.end(),
            d_giant_ends.begin() + prev_giant_count, transform_fn);

        d_giant_filter_inds.resize(d_giant_peaks.size(), filter_idx);
        d_giant_dm_inds.resize(d_giant_peaks.size(), dm_idx);
        // Note: This could be used to track total member samples if desired
        d_giant_members.resize(d_giant_peaks.size(), 1);

        stop_timer(giants_timer);

        // Bail if the candidate rate is too high
        hd_size total_giant_count = d_giant_peaks.size();
        hd_float data_length_mins = nsamps * pl->params.dt / 60.0;

        if (total_timer.getTime() > 3500 || total_giant_count > 1000000) {
          too_many_giants = true;
          float searched = ((float)dm_idx * 100) / (float)dm_count;
          std::cout
              << "WARNING: exceeded processing time / giant count: 3.5s, 10k "
                 "DM ["
              << dm_list[dm_idx] << "] space searched " << searched << "%"
              << std::endl;
          break;
        }

      } // close filter width loop

    } // close DM loop

  } // close gulp_idx condition

  hd_size giant_count = d_giant_peaks.size();
  std::cout << "Giant count = " << giant_count << std::endl;
  std::cout << "final_space_searched " << dm_list[dm_idx_output] << std::endl;

  // Fill in the data
  write_candidates(pl, nsamps, first_idx, nsamps_computed, dm_list,
                   d_giant_peaks, d_giant_inds, d_giant_filter_inds,
                   d_giant_dm_inds);

  stop_timer(candidates_timer);

  stop_timer(total_timer);

  std::cout << "Mem alloc time:          " << memory_timer.getTime()
            << std::endl;
  std::cout << "0-DM cleaning time:      " << clean_timer.getTime()
            << std::endl;
  std::cout << "Dedispersion time:       " << dedisp_timer.getTime()
            << std::endl;
  std::cout << "Copy time:               " << copy_timer.getTime() << std::endl;
  std::cout << "Baselining time:         " << baseline_timer.getTime()
            << std::endl;
  std::cout << "Normalisation time:      " << normalise_timer.getTime()
            << std::endl;
  std::cout << "Filtering time:          " << filter_timer.getTime()
            << std::endl;
  std::cout << "Find giants time:        " << giants_timer.getTime()
            << std::endl;
  std::cout << "Process candidates time: " << candidates_timer.getTime()
            << std::endl;
  std::cout << "Total time:              " << total_timer.getTime()
            << std::endl;

  if (too_many_giants) {
    return HD_TOO_MANY_EVENTS;
  } else {
    return HD_NO_ERROR;
  }
}

void hd_destroy_pipeline(hd_pipeline pipeline) {
  if (pipeline->params.verbosity >= 2) {
    std::cout << "\tDeleting pipeline object..." << std::endl;
  }

  dedisp_destroy_plan(pipeline->dedispersion_plan);

  // Close the writing contexts
  if (pipeline->params.coincidencer_host != NULL &&
      pipeline->params.coincidencer_port != -1) {
    ::close(pipeline->socket);
  } else {
    pipeline->file.close();
  }

  // Note: This assumes memory owned by pipeline cleans itself up
  if (pipeline) {
    delete pipeline;
  }
}
