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

#include <lemon/lp_utils.h>

namespace lemon {

  LpResultMap lpResultMap(const LpSolverBase& lp) {
    return LpResultMap(lp);
  }

  LpColNameMap lpColNameMap(const LpSolverBase& lp) {
    return LpColNameMap(lp);
  }

  LpColNameWriteMap lpColNameMap(LpSolverBase& lp) {
    return LpColNameWriteMap(lp);
  }

}
