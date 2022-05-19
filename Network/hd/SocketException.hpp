/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#ifndef SOCKET_EXCEPTION_HPP
#define SOCKET_EXCEPTION_HPP

#include <string>

class SocketException {
public:
  SocketException(std::string s) : m_s(s){};
  ~SocketException(){};

  std::string description() { return m_s; }

private:
  std::string m_s;
};

#endif
