/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#ifndef REMOVE_BASELINE_HPP
#define REMOVE_BASELINE_HPP

#include "error.hpp"
#include "types.hpp"

#include <boost/shared_ptr.hpp>

struct RemoveBaselinePlan_impl;

struct RemoveBaselinePlan {
  RemoveBaselinePlan();
  hd_error exec(hd_float *d_data, hd_size count, hd_size smooth_radius);

private:
  boost::shared_ptr<RemoveBaselinePlan_impl> m_impl;
};

#endif