/***************************************************************************
 *
 *   Copyright (C) 2012 by Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#ifndef PSR_DADA_RING_BUFFER_HPP
#define PSR_DADA_RING_BUFFER_HPP

#include <dada_hdu.h>

#include "DataSource.hpp"

class PSRDadaRingBuffer : public DataSource {

public:
  PSRDadaRingBuffer(key_t dada_id);
  ~PSRDadaRingBuffer();

  bool connect();
  bool disconnect();
  bool read_header();

  bool get_error() const { return dada_error; }
  size_t get_data(size_t nsamps, char *data);

private:
  key_t dada_key;

  dada_hdu_t *hdu;

  bool dada_error;

  bool printed_first_line;
};

#endif