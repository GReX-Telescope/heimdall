/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#ifndef MERGE_CANDIDATES_HPP
#define MEREG_CANDIDATES_HPP

#include "error.hpp"
#include "types.hpp"

hd_error merge_candidates(hd_size count, hd_size *d_labels,
                          ConstRawCandidates d_cands, RawCandidates d_groups);

#endif