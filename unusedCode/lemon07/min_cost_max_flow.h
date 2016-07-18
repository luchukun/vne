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

#ifndef LEMON_MIN_COST_MAX_FLOW_H
#define LEMON_MIN_COST_MAX_FLOW_H

/// \ingroup min_cost_flow
///
/// \file
/// \brief An efficient algorithm for finding a minimum cost maximum flow.

#include <lemon/preflow.h>
#include <lemon/network_simplex.h>
#include <lemon/maps.h>

namespace lemon {

  /// \addtogroup min_cost_flow
  /// @{

  /// \brief An efficient algorithm for finding a minimum cost
  /// maximum flow.
  ///
  /// \ref MinCostMaxFlow implements an efficient algorithm for
  /// finding a maximum flow having minimal total cost from a given
  /// source node to a given target node in a directed graph.
  ///
  /// \ref MinCostMaxFlow uses \ref Preflow for finding the maximum
  /// flow value and \ref NetworkSimplex for finding a minimum cost
  /// flow of that value.
  /// According to our benchmark tests \ref Preflow is generally the
  /// most efficient algorithm for the maximum flow problem and
  /// \ref NetworkSimplex is the most efficient for the minimum cost
  /// flow problem in LEMON.
  ///
  /// \tparam Graph The directed graph type the algorithm runs on.
  /// \tparam CapacityMap The type of the capacity (upper bound) map.
  /// \tparam CostMap The type of the cost (length) map.
  ///
  /// \warning
  /// - Edge capacities and costs should be \e non-negative \e integers.
  /// - \c CapacityMap::Value must be convertible to \c CostMap::Value.
  /// - \c CostMap::Value must be signed type.
  ///
  /// \author Peter Kovacs

  template < typename Graph,
             typename CapacityMap = typename Graph::template EdgeMap<int>,
             typename CostMap = typename Graph::template EdgeMap<int> >
  class MinCostMaxFlow
  {
    GRAPH_TYPEDEFS(typename Graph);

    typedef typename CapacityMap::Value Capacity;
    typedef typename CostMap::Value Cost;
    typedef typename Graph::template NodeMap<Cost> SupplyMap;

    typedef Preflow<Graph, CapacityMap> MaxFlowImpl;
    typedef NetworkSimplex< Graph, CapacityMap, CapacityMap,
                            CostMap, SupplyMap > MinCostFlowImpl;

  public:

    /// The type of the flow map.
    typedef typename Graph::template EdgeMap<Capacity> FlowMap;
    /// The type of the potential map.
    typedef typename Graph::template NodeMap<Cost> PotentialMap;

  private:

    // The directed graph the algorithm runs on
    const Graph &_graph;
    // The capacity map
    const CapacityMap &_capacity;
    // The cost map
    const CostMap &_cost;

    // Edge map of the found flow
    FlowMap *_flow;
    bool _local_flow;
    // Node map of the current potentials
    PotentialMap *_potential;
    bool _local_potential;

    // The source node
    Node _source;
    // The target node
    Node _target;

  public:

    /// \brief Constructor.
    ///
    /// Constructor.
    ///
    /// \param graph The directed graph the algorithm runs on.
    /// \param capacity The capacities (upper bounds) of the edges.
    /// \param cost The cost (length) values of the edges.
    /// \param s The source node.
    /// \param t The target node.
    MinCostMaxFlow( const Graph &graph,
                    const CapacityMap &capacity,
                    const CostMap &cost,
                    Node s, Node t ) :
      _graph(graph), _capacity(capacity), _cost(cost), _flow(0),
      _local_flow(false), _potential(0), _local_potential(false),
      _source(s), _target(t) {}

    /// Destructor.
    ~MinCostMaxFlow() {
      if (_local_flow) delete _flow;
      if (_local_potential) delete _potential;
    }

    /// \brief Sets the flow map.
    ///
    /// Sets the flow map.
    ///
    /// \return \c (*this)
    MinCostMaxFlow& flowMap(FlowMap &map) {
      if (_local_flow) {
        delete _flow;
        _local_flow = false;
      }
      _flow = &map;
      return *this;
    }

    /// \brief Sets the potential map.
    ///
    /// Sets the potential map.
    ///
    /// \return \c (*this)
    MinCostMaxFlow& potentialMap(PotentialMap &map) {
      if (_local_potential) {
        delete _potential;
        _local_potential = false;
      }
      _potential = &map;
      return *this;
    }

    /// \name Execution control
    /// The only way to execute the algorithm is to call the run()
    /// function.

    /// @{

    /// \brief Runs the algorithm.
    ///
    /// Runs the algorithm.
    void run() {
      // Initializing maps
      if (!_flow) {
        _flow = new FlowMap(_graph);
        _local_flow = true;
      }
      if (!_potential) {
        _potential = new PotentialMap(_graph);
        _local_potential = true;
      }
      // Running Preflow
      MaxFlowImpl preflow(_graph, _capacity, _source, _target);
      preflow.flowMap(*_flow).runMinCut();
      // Running NetworkSimplex
      MinCostFlowImpl mcf( _graph, _capacity, _cost,
                           _source, _target, preflow.flowValue() );
      mcf.flowMap(*_flow).potentialMap(*_potential).run();
    }

    /// @}

    /// \name Query Functions
    /// The result of the algorithm can be obtained using these
    /// functions.
    /// \n run() must be called before using them.

    /// @{

    /// \brief Returns a const reference to the edge map storing the
    /// found flow.
    ///
    /// Returns a const reference to the edge map storing the found flow.
    ///
    /// \pre \ref run() must be called before using this function.
    const FlowMap& flowMap() const {
      return *_flow;
    }

    /// \brief Returns a const reference to the node map storing the
    /// found potentials (the dual solution).
    ///
    /// Returns a const reference to the node map storing the found
    /// potentials (the dual solution).
    ///
    /// \pre \ref run() must be called before using this function.
    const PotentialMap& potentialMap() const {
      return *_potential;
    }

    /// \brief Returns the flow on the given edge.
    ///
    /// Returns the flow on the given edge.
    ///
    /// \pre \ref run() must be called before using this function.
    Capacity flow(const Edge& edge) const {
      return (*_flow)[edge];
    }

    /// \brief Returns the potential of the given node.
    ///
    /// Returns the potential of the given node.
    ///
    /// \pre \ref run() must be called before using this function.
    Cost potential(const Node& node) const {
      return (*_potential)[node];
    }

    /// \brief Returns the total cost of the found flow.
    ///
    /// Returns the total cost of the found flow. The complexity of the
    /// function is \f$ O(e) \f$.
    ///
    /// \pre \ref run() must be called before using this function.
    Cost totalCost() const {
      Cost c = 0;
      for (EdgeIt e(_graph); e != INVALID; ++e)
        c += (*_flow)[e] * _cost[e];
      return c;
    }

    /// @}

  }; //class MinCostMaxFlow

  ///@}

} //namespace lemon

#endif //LEMON_MIN_COST_MAX_FLOW_H
