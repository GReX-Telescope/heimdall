/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#ifndef PARSE_COMMAND_LINE_HPP
#define PARSE_COMMAND_LINE_HPP

#ifdef __cplusplus
extern "C" {
#endif

#include "params.hpp"

int hd_parse_command_line(int argc, char *argv[], hd_params *params);
void hd_print_usage();

#ifdef __cplusplus
} // closing brace for extern "C"
#endif

#endif