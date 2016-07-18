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

#ifndef LEMON_MIN_COST_FLOW_H
#define LEMON_MIN_COST_FLOW_H

/// \ingroup min_cost_flow
///
/// \file
/// \brief An efficient algorithm for finding a minimum cost flow.

#include <lemon/network_simplex.h>

namespace lemon {

  /// \addtogroup min_cost_flow
  /// @{

  /// \brief An efficient algorithm for finding a minimum cost flow.
  ///
  /// \ref MinCostFlow provides an efficient algorithm for finding
  /// a minimum cost flow.
  ///
  /// This class is just an alias for \ref NetworkSimplex,
  /// which is the most efficient algorithm for the minimum cost
  /// flow problem in LEMON according to our benchmark tests.
  /// For the detailed documentation of this class see
  /// \ref NetworkSimplex.
  ///
  /// There are four implementations for the minimum cost flow problem,
  /// which can be used exactly the same way.
  /// - \ref CapacityScaling
  /// - \ref CostScaling
  /// - \ref CycleCanceling
  /// - \ref NetworkSimplex
  ///
  /// \tparam Graph The directed graph type the algorithm runs on.
  /// \tparam LowerMap The type of the lower bound map.
  /// \tparam CapacityMap The type of the capacity (upper bound) map.
  /// \tparam CostMap The type of the cost (length) map.
  /// \tparam SupplyMap The type of the supply map.
  ///
  /// \warning
  /// - Edge capacities and costs should be \e non-negative \e integers.
  /// - Supply values should be \e signed \e integers.
  /// - The value types of the maps should be convertible to each other.
  /// - \c CostMap::Value must be signed type.
  ///
  /// \author Peter Kovacs

  template < typename Graph,
             typename LowerMap = typename Graph::template EdgeMap<int>,
             typename CapacityMap = typename Graph::template EdgeMap<int>,
             typename CostMap = typename Graph::template EdgeMap<int>,
             typename SupplyMap = typename Graph::template NodeMap<int> >
  class MinCostFlow :
    public NetworkSimplex< Graph, LowerMap, CapacityMap,
                           CostMap, SupplyMap >
  {
    typedef NetworkSimplex< Graph, LowerMap, CapacityMap,
                            CostMap, SupplyMap > MinCostFlowImpl;
    typedef typename Graph::Node Node;
    typedef typename SupplyMap::Value Supply;

  public:

    /// General constructor (with lower bounds).
    MinCostFlow( const Graph &graph,
                 const LowerMap &lower,
                 const CapacityMap &capacity,
                 const CostMap &cost,
                 const SupplyMap &supply ) :
      MinCostFlowImpl(graph, lower, capacity, cost, supply) {}

    /// General constructor of the class (without lower bounds).
    MinCostFlow( const Graph &graph,
                 const CapacityMap &capacity,
                 const CostMap &cost,
                 const SupplyMap &supply ) :
      MinCostFlowImpl(graph, capacity, cost, supply) {}

    /// Simple constructor (with lower bounds).
    MinCostFlow( const Graph &graph,
                 const LowerMap &lower,
                 const CapacityMap &capacity,
                 const CostMap &cost,
                 Node s, Node t,
                 Supply flow_value ) :
      MinCostFlowImpl( graph, lower, capacity, cost,
                       s, t, flow_value ) {}

    /// Simple constructor (without lower bounds).
    MinCostFlow( const Graph &graph,
                 const CapacityMap &capacity,
                 const CostMap &cost,
                 Node s, Node t,
                 Supply flow_value ) :
      MinCostFlowImpl( graph, capacity, cost,
                       s, t, flow_value ) {}

  }; //class MinCostFlow

  ///@}

} //namespace lemon

#endif //LEMON_MIN_COST_FLOW_H
