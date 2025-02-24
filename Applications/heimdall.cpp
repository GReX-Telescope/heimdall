/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <cstdlib>
#include <cuda_profiler_api.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string.h>
#include <unistd.h>
#include <vector>

#include "hd/cached_allocator.cuh"
#include "hd/default_params.hpp"
#include "hd/error.hpp"
#include "hd/find_giants.hpp"
#include "hd/parse_command_line.hpp"
#include "hd/pipeline.hpp"
#include "hd/types.hpp"
#include <dedisp.h>

// input formats supported
#include "hd/DataSource.hpp"
#include "hd/SigprocFile.hpp"
#ifdef HAVE_PSRDADA
#include "hd/PSRDadaRingBuffer.hpp"
#endif

#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>

#include "hd/stopwatch.hpp"

namespace logging = boost::log;

void log_init(int verbosity) {
  if (verbosity >= 4) {
    logging::core::get()->set_filter(logging::trivial::severity >=
                                     logging::trivial::trace);
  } else if (verbosity >= 1) {
    logging::core::get()->set_filter(logging::trivial::severity >=
                                     logging::trivial::debug);
  } else {
    logging::core::get()->set_filter(logging::trivial::severity >=
                                     logging::trivial::info);
  }
}

int main(int argc, char *argv[]) {
  hd_params params;
  hd_set_default_params(&params);
  int ok = hd_parse_command_line(argc, argv, &params);

  log_init(params.verbosity);

  size_t nsamps_gulp = params.nsamps_gulp;

  if (ok < 0)
    return 1;
  DataSource *data_source = 0;
#ifdef HAVE_PSRDADA
  if (params.dada_id != 0) {

    if (params.verbosity)
      cerr << "Creating PSRDADA client" << endl;

    PSRDadaRingBuffer *d = new PSRDadaRingBuffer(params.dada_id);

    // Read from psrdada ring buffer
    if (!d || d->get_error()) {
      cerr << "ERROR: Failed to initialise connection to psrdada" << endl;
      return -1;
    }

    if (params.verbosity)
      cerr << "Connecting to ring buffer" << endl;

    // connect to PSRDADA ring buffer
    if (!d->connect()) {
      cerr << "ERROR: Failed to connection to psrdada ring buffer" << endl;
      return -1;
    }

    if (params.verbosity)
      cerr << "Waiting for next header / data" << endl;

    // wait for and then read next PSRDADA header/observation
    if (!d->read_header()) {
      cerr << "ERROR: Failed to connection to psrdada ring buffer" << endl;
      return -1;
    }

    data_source = (DataSource *)d;
    if (!params.override_beam)
      params.beam = d->get_beam() - 1;
  } else {
#endif
    // Read from filterbank file
    data_source = new SigprocFile(params.sigproc_file, params.fswap);
    if (!data_source || data_source->get_error()) {
      cerr << "ERROR: Failed to open data file" << endl;
      return -1;
    }
#ifdef HAVE_PSRDADA
  }
#endif

  if (!params.override_beam)
    if (data_source->get_beam() > 0)
      params.beam = data_source->get_beam() - 1;
    else
      params.beam = 0;

  params.f0 = data_source->get_f0();
  params.df = data_source->get_df();
  params.dt = data_source->get_tsamp();

  if (params.verbosity > 0)
    cout << "processing beam " << (params.beam + 1) << endl;

  float tsamp = data_source->get_tsamp() / 1000000;
  size_t stride = data_source->get_stride();
  size_t nbits = data_source->get_nbit();

  params.nchans = data_source->get_nchan();
  params.utc_start = data_source->get_utc_start();
  params.spectra_per_second = data_source->get_spectra_rate();

  bool stop_requested = false;

  // Create the pipeline object
  // --------------------------
  hd_pipeline pipeline;
  hd_error error;
  error = hd_create_pipeline(&pipeline, params);
  if (error != HD_NO_ERROR) {
    cerr << "ERROR: Pipeline creation failed" << endl;
    cerr << "       " << hd_get_error_string(error) << endl;
    return -1;
  }
  // --------------------------

  size_t derror;
  dedisp_plan dedispersion_plan;
  derror = dedisp_create_plan(&dedispersion_plan, params.nchans, params.dt,
                              params.f0, params.df);

  derror =
      dedisp_generate_dm_list(dedispersion_plan, params.dm_min, params.dm_max,
                              params.dm_pulse_width, params.dm_tol);

  size_t max_delay = dedisp_get_max_delay(dedispersion_plan);
  size_t boxcar_max = params.boxcar_max;

  if (params.verbosity >= 2)
    cout << "allocating filterbank data vector for "
         << (nsamps_gulp + max_delay + boxcar_max) << " samples with size "
         << ((nsamps_gulp + max_delay + boxcar_max) * stride * params.nbeams)
         << " bytes"
         << " with " << params.nbeams << " beams." << endl;
  std::vector<hd_byte> filterbank((nsamps_gulp + max_delay + boxcar_max) *
                                  stride * params.nbeams);

  if (params.verbosity >= 1) {
    cout << "Beginning data processing, requesting " << nsamps_gulp
         << " samples" << endl;
  }

  int fseq = 0;
  size_t cur_nsamps = 0;
  size_t total_nsamps = 0;
  size_t nsamps_read = 0;
  size_t overlap = 0;
  size_t gulp_idx = 0;
  hd_size nsamps_processed;

  // for dealing with memory allocation a priori
  error =
      hd_execute(pipeline, &filterbank[0], nsamps_gulp + max_delay + boxcar_max,
                 nbits, total_nsamps, cur_nsamps, &nsamps_processed, 0);

  for (int i = 0; i < params.nbeams; i++) {
    for (int j = i * (nsamps_gulp + max_delay + boxcar_max) * stride;
         j <
         (i * (nsamps_gulp + max_delay + boxcar_max) + max_delay + boxcar_max) *
             stride;
         j++)
      filterbank[j] = 128;
    nsamps_read += data_source->get_data(
        nsamps_gulp,
        (char *)&filterbank[(i * (nsamps_gulp + max_delay + boxcar_max) +
                             max_delay + boxcar_max) *
                            stride]);
  }
  nsamps_read = nsamps_read / params.nbeams;

  while (nsamps_read && !stop_requested) {

    if (params.verbosity >= 1) {
      cout << "Executing pipeline on new gulp of "
           << nsamps_gulp + max_delay + boxcar_max << " samples..." << endl;
      cout << "total_nsamps =" << total_nsamps << endl;
    }

    cudaProfilerStart();
    error = hd_execute(pipeline, &filterbank[0],
                       nsamps_gulp + max_delay + boxcar_max, nbits,
                       total_nsamps, cur_nsamps, &nsamps_processed, gulp_idx);
    cudaProfilerStop();
    gulp_idx++;
    if (error == HD_NO_ERROR) {
      if (params.verbosity >= 1)
        cout << "Processed " << nsamps_processed << " samples." << endl;
    } else if (error == HD_TOO_MANY_EVENTS) {
      if (params.verbosity >= 1)
        cerr
            << "WARNING: hd_execute produces too many events, some data skipped"
            << endl;
    } else {
      cerr << "ERROR: Pipeline execution failed" << endl;
      cerr << "       " << hd_get_error_string(error) << endl;
      hd_destroy_pipeline(pipeline);
      dedisp_destroy_plan(dedispersion_plan);
      delete data_source;
      return -1;
    }

    if (params.verbosity >= 1)
      cout << "Main: nsamps_processed=" << nsamps_processed << endl;

    if (total_nsamps == 0)
      total_nsamps += nsamps_gulp - max_delay - boxcar_max;
    else
      total_nsamps += nsamps_processed;

    for (int i = 0; i < params.nbeams; i++) {
      std::copy(
          &filterbank[((i * (nsamps_gulp + max_delay + boxcar_max) +
                        nsamps_gulp)) *
                      stride],
          &filterbank[((i + 1) * (nsamps_gulp + max_delay + boxcar_max)) *
                      stride],
          &filterbank[i * (nsamps_gulp + max_delay + boxcar_max) * stride]);
      nsamps_read = data_source->get_data(
          (nsamps_gulp),
          (char *)&filterbank[(max_delay + boxcar_max +
                               i * (nsamps_gulp + max_delay + boxcar_max)) *
                              stride]);
    }

    // at the end of data, never execute the pipeline
    if (nsamps_read < (nsamps_gulp - overlap))
      stop_requested = 1;
  }

  if (params.verbosity >= 1) {
    cout << "Successfully processed a total of " << total_nsamps << " samples."
         << endl;
    cout << "Shutting down..." << endl;
  }

  hd_destroy_pipeline(pipeline);
  dedisp_destroy_plan(dedispersion_plan);
  delete data_source;
  free((void *)(params.sigproc_file));
  free((void *)(params.output_dir));
  free((void *)(params.coincidencer_host));

  if (params.verbosity >= 1) {
    cout << "All done." << endl;
  }
  g_allocator.free_all();
}
