/* -*- C++ -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library
 *
 * Copyright (C) 2003-2008
 * Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
 * (Egervary Research Group on Combinatorial Optimization, EGRES).
 *
 * Permission to use, modify and distribute this software is granted
 * provided that this copyright notice appears in all copies. For
 * precise terms see the accompanying LICENSE file.
 *
 * This software is provided "AS IS" with no warranty of any kind,
 * express or implied, and with no claim as to its suitability for any
 * purpose.
 *
 */

///\file
///\brief The implementation of the LP solver interface.

#include <lemon/lp_base.h>
namespace lemon {
  
  const LpSolverBase::Value
  LpSolverBase::INF = std::numeric_limits<Value>::infinity();
  const LpSolverBase::Value
  LpSolverBase::NaN = std::numeric_limits<Value>::quiet_NaN();

//   const LpSolverBase::Constr::Value
//   LpSolverBase::Constr::INF = std::numeric_limits<Value>::infinity();
//   const LpSolverBase::Constr::Value
//   LpSolverBase::Constr::NaN = std::numeric_limits<Value>::quiet_NaN();
  
} //namespace lemon
