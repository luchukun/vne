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

#ifndef LEMON_CYCLE_CANCELING_H
#define LEMON_CYCLE_CANCELING_H

/// \ingroup min_cost_flow
///
/// \file
/// \brief Cycle-canceling algorithm for finding a minimum cost flow.

#include <vector>
#include <lemon/graph_adaptor.h>
#include <lemon/path.h>

#include <lemon/circulation.h>
#include <lemon/bellman_ford.h>
#include <lemon/min_mean_cycle.h>

namespace lemon {

  /// \addtogroup min_cost_flow
  /// @{

  /// \brief Implementation of a cycle-canceling algorithm for
  /// finding a minimum cost flow.
  ///
  /// \ref CycleCanceling implements a cycle-canceling algorithm for
  /// finding a minimum cost flow.
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
  /// \note By default the \ref BellmanFord "Bellman-Ford" algorithm is
  /// used for negative cycle detection with limited iteration number.
  /// However \ref CycleCanceling also provides the "Minimum Mean
  /// Cycle-Canceling" algorithm, which is \e strongly \e polynomial,
  /// but rather slower in practice.
  /// To use this version of the algorithm, call \ref run() with \c true
  /// parameter.
  ///
  /// \author Peter Kovacs

  template < typename Graph,
             typename LowerMap = typename Graph::template EdgeMap<int>,
             typename CapacityMap = typename Graph::template EdgeMap<int>,
             typename CostMap = typename Graph::template EdgeMap<int>,
             typename SupplyMap = typename Graph::template NodeMap<int> >
  class CycleCanceling
  {
    GRAPH_TYPEDEFS(typename Graph);

    typedef typename CapacityMap::Value Capacity;
    typedef typename CostMap::Value Cost;
    typedef typename SupplyMap::Value Supply;
    typedef typename Graph::template EdgeMap<Capacity> CapacityEdgeMap;
    typedef typename Graph::template NodeMap<Supply> SupplyNodeMap;

    typedef ResGraphAdaptor< const Graph, Capacity,
                             CapacityEdgeMap, CapacityEdgeMap > ResGraph;
    typedef typename ResGraph::Node ResNode;
    typedef typename ResGraph::NodeIt ResNodeIt;
    typedef typename ResGraph::Edge ResEdge;
    typedef typename ResGraph::EdgeIt ResEdgeIt;

  public:

    /// The type of the flow map.
    typedef typename Graph::template EdgeMap<Capacity> FlowMap;
    /// The type of the potential map.
    typedef typename Graph::template NodeMap<Cost> PotentialMap;

  private:

    /// \brief Map adaptor class for handling residual edge costs.
    ///
    /// \ref ResidualCostMap is a map adaptor class for handling
    /// residual edge costs.
    class ResidualCostMap : public MapBase<ResEdge, Cost>
    {
    private:

      const CostMap &_cost_map;

    public:

      ///\e
      ResidualCostMap(const CostMap &cost_map) : _cost_map(cost_map) {}

      ///\e
      Cost operator[](const ResEdge &e) const {
        return ResGraph::forward(e) ? _cost_map[e] : -_cost_map[e];
      }

    }; //class ResidualCostMap

  private:

    // The maximum number of iterations for the first execution of the
    // Bellman-Ford algorithm. It should be at least 2.
    static const int BF_FIRST_LIMIT  = 2;
    // The iteration limit for the Bellman-Ford algorithm is multiplied
    // by BF_LIMIT_FACTOR/100 in every round.
    static const int BF_LIMIT_FACTOR = 150;

  private:

    // The directed graph the algorithm runs on
    const Graph &_graph;
    // The original lower bound map
    const LowerMap *_lower;
    // The modified capacity map
    CapacityEdgeMap _capacity;
    // The original cost map
    const CostMap &_cost;
    // The modified supply map
    SupplyNodeMap _supply;
    bool _valid_supply;

    // Edge map of the current flow
    FlowMap *_flow;
    bool _local_flow;
    // Node map of the current potentials
    PotentialMap *_potential;
    bool _local_potential;

    // The residual graph
    ResGraph *_res_graph;
    // The residual cost map
    ResidualCostMap _res_cost;

  public:

    /// \brief General constructor (with lower bounds).
    ///
    /// General constructor (with lower bounds).
    ///
    /// \param graph The directed graph the algorithm runs on.
    /// \param lower The lower bounds of the edges.
    /// \param capacity The capacities (upper bounds) of the edges.
    /// \param cost The cost (length) values of the edges.
    /// \param supply The supply values of the nodes (signed).
    CycleCanceling( const Graph &graph,
                    const LowerMap &lower,
                    const CapacityMap &capacity,
                    const CostMap &cost,
                    const SupplyMap &supply ) :
      _graph(graph), _lower(&lower), _capacity(graph), _cost(cost),
      _supply(graph), _flow(0), _local_flow(false),
      _potential(0), _local_potential(false), _res_cost(_cost)
    {
      // Removing non-zero lower bounds
      _capacity = subMap(capacity, lower);
      Supply sum = 0;
      for (NodeIt n(_graph); n != INVALID; ++n) {
        Supply s = supply[n];
        for (InEdgeIt e(_graph, n); e != INVALID; ++e)
          s += lower[e];
        for (OutEdgeIt e(_graph, n); e != INVALID; ++e)
          s -= lower[e];
        sum += (_supply[n] = s);
      }
      _valid_supply = sum == 0;
    }

    /// \brief General constructor (without lower bounds).
    ///
    /// General constructor (without lower bounds).
    ///
    /// \param graph The directed graph the algorithm runs on.
    /// \param capacity The capacities (upper bounds) of the edges.
    /// \param cost The cost (length) values of the edges.
    /// \param supply The supply values of the nodes (signed).
    CycleCanceling( const Graph &graph,
                    const CapacityMap &capacity,
                    const CostMap &cost,
                    const SupplyMap &supply ) :
      _graph(graph), _lower(NULL), _capacity(capacity), _cost(cost),
      _supply(supply), _flow(0), _local_flow(false),
      _potential(0), _local_potential(false), _res_cost(_cost)
    {
      // Checking the sum of supply values
      Supply sum = 0;
      for (NodeIt n(_graph); n != INVALID; ++n) sum += _supply[n];
      _valid_supply = sum == 0;
    }

    /// \brief Simple constructor (with lower bounds).
    ///
    /// Simple constructor (with lower bounds).
    ///
    /// \param graph The directed graph the algorithm runs on.
    /// \param lower The lower bounds of the edges.
    /// \param capacity The capacities (upper bounds) of the edges.
    /// \param cost The cost (length) values of the edges.
    /// \param s The source node.
    /// \param t The target node.
    /// \param flow_value The required amount of flow from node \c s
    /// to node \c t (i.e. the supply of \c s and the demand of \c t).
    CycleCanceling( const Graph &graph,
                    const LowerMap &lower,
                    const CapacityMap &capacity,
                    const CostMap &cost,
                    Node s, Node t,
                    Supply flow_value ) :
      _graph(graph), _lower(&lower), _capacity(graph), _cost(cost),
      _supply(graph), _flow(0), _local_flow(false),
      _potential(0), _local_potential(false), _res_cost(_cost)
    {
      // Removing non-zero lower bounds
      _capacity = subMap(capacity, lower);
      for (NodeIt n(_graph); n != INVALID; ++n) {
        Supply sum = 0;
        if (n == s) sum =  flow_value;
        if (n == t) sum = -flow_value;
        for (InEdgeIt e(_graph, n); e != INVALID; ++e)
          sum += lower[e];
        for (OutEdgeIt e(_graph, n); e != INVALID; ++e)
          sum -= lower[e];
        _supply[n] = sum;
      }
      _valid_supply = true;
    }

    /// \brief Simple constructor (without lower bounds).
    ///
    /// Simple constructor (without lower bounds).
    ///
    /// \param graph The directed graph the algorithm runs on.
    /// \param capacity The capacities (upper bounds) of the edges.
    /// \param cost The cost (length) values of the edges.
    /// \param s The source node.
    /// \param t The target node.
    /// \param flow_value The required amount of flow from node \c s
    /// to node \c t (i.e. the supply of \c s and the demand of \c t).
    CycleCanceling( const Graph &graph,
                    const CapacityMap &capacity,
                    const CostMap &cost,
                    Node s, Node t,
                    Supply flow_value ) :
      _graph(graph), _lower(NULL), _capacity(capacity), _cost(cost),
      _supply(graph, 0), _flow(0), _local_flow(false),
      _potential(0), _local_potential(false), _res_cost(_cost)
    {
      _supply[s] =  flow_value;
      _supply[t] = -flow_value;
      _valid_supply = true;
    }

    /// Destructor.
    ~CycleCanceling() {
      if (_local_flow) delete _flow;
      if (_local_potential) delete _potential;
      delete _res_graph;
    }

    /// \brief Sets the flow map.
    ///
    /// Sets the flow map.
    ///
    /// \return \c (*this)
    CycleCanceling& flowMap(FlowMap &map) {
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
    CycleCanceling& potentialMap(PotentialMap &map) {
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
    ///
    /// \param min_mean_cc Set this parameter to \c true to run the
    /// "Minimum Mean Cycle-Canceling" algorithm, which is strongly
    /// polynomial, but rather slower in practice.
    ///
    /// \return \c true if a feasible flow can be found.
    bool run(bool min_mean_cc = false) {
      return init() && start(min_mean_cc);
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

  private:

    /// Initializes the algorithm.
    bool init() {
      if (!_valid_supply) return false;

      // Initializing flow and potential maps
      if (!_flow) {
        _flow = new FlowMap(_graph);
        _local_flow = true;
      }
      if (!_potential) {
        _potential = new PotentialMap(_graph);
        _local_potential = true;
      }

      _res_graph = new ResGraph(_graph, _capacity, *_flow);

      // Finding a feasible flow using Circulation
      Circulation< Graph, ConstMap<Edge, Capacity>, CapacityEdgeMap,
                   SupplyMap >
        circulation( _graph, constMap<Edge>(Capacity(0)), _capacity,
                     _supply );
      return circulation.flowMap(*_flow).run();
    }

    bool start(bool min_mean_cc) {
      if (min_mean_cc)
        startMinMean();
      else
        start();

      // Handling non-zero lower bounds
      if (_lower) {
        for (EdgeIt e(_graph); e != INVALID; ++e)
          (*_flow)[e] += (*_lower)[e];
      }
      return true;
    }

    /// \brief Executes the algorithm using \ref BellmanFord.
    ///
    /// Executes the algorithm using the \ref BellmanFord
    /// "Bellman-Ford" algorithm for negative cycle detection with
    /// successively larger limit for the number of iterations.
    void start() {
      typename BellmanFord<ResGraph, ResidualCostMap>::PredMap pred(*_res_graph);
      typename ResGraph::template NodeMap<int> visited(*_res_graph);
      std::vector<ResEdge> cycle;
      int node_num = countNodes(_graph);

      int length_bound = BF_FIRST_LIMIT;
      bool optimal = false;
      while (!optimal) {
        BellmanFord<ResGraph, ResidualCostMap> bf(*_res_graph, _res_cost);
        bf.predMap(pred);
        bf.init(0);
        int iter_num = 0;
        bool cycle_found = false;
        while (!cycle_found) {
          int curr_iter_num = iter_num + length_bound <= node_num ?
                              length_bound : node_num - iter_num;
          iter_num += curr_iter_num;
          int real_iter_num = curr_iter_num;
          for (int i = 0; i < curr_iter_num; ++i) {
            if (bf.processNextWeakRound()) {
              real_iter_num = i;
              break;
            }
          }
          if (real_iter_num < curr_iter_num) {
            // Optimal flow is found
            optimal = true;
            // Setting node potentials
            for (NodeIt n(_graph); n != INVALID; ++n)
              (*_potential)[n] = bf.dist(n);
            break;
          } else {
            // Searching for node disjoint negative cycles
            for (ResNodeIt n(*_res_graph); n != INVALID; ++n)
              visited[n] = 0;
            int id = 0;
            for (ResNodeIt n(*_res_graph); n != INVALID; ++n) {
              if (visited[n] > 0) continue;
              visited[n] = ++id;
              ResNode u = pred[n] == INVALID ?
                          INVALID : _res_graph->source(pred[n]);
              while (u != INVALID && visited[u] == 0) {
                visited[u] = id;
                u = pred[u] == INVALID ?
                    INVALID : _res_graph->source(pred[u]);
              }
              if (u != INVALID && visited[u] == id) {
                // Finding the negative cycle
                cycle_found = true;
                cycle.clear();
                ResEdge e = pred[u];
                cycle.push_back(e);
                Capacity d = _res_graph->rescap(e);
                while (_res_graph->source(e) != u) {
                  cycle.push_back(e = pred[_res_graph->source(e)]);
                  if (_res_graph->rescap(e) < d)
                    d = _res_graph->rescap(e);
                }

                // Augmenting along the cycle
                for (int i = 0; i < int(cycle.size()); ++i)
                  _res_graph->augment(cycle[i], d);
              }
            }
          }

          if (!cycle_found)
            length_bound = length_bound * BF_LIMIT_FACTOR / 100;
        }
      }
    }

    /// \brief Executes the algorithm using \ref MinMeanCycle.
    ///
    /// Executes the algorithm using \ref MinMeanCycle for negative
    /// cycle detection.
    void startMinMean() {
      typedef Path<ResGraph> ResPath;
      MinMeanCycle<ResGraph, ResidualCostMap> mmc(*_res_graph, _res_cost);
      ResPath cycle;

      mmc.cyclePath(cycle).init();
      if (mmc.findMinMean()) {
        while (mmc.cycleLength() < 0) {
          // Finding the cycle
          mmc.findCycle();

          // Finding the largest flow amount that can be augmented
          // along the cycle
          Capacity delta = 0;
          for (typename ResPath::EdgeIt e(cycle); e != INVALID; ++e) {
            if (delta == 0 || _res_graph->rescap(e) < delta)
              delta = _res_graph->rescap(e);
          }

          // Augmenting along the cycle
          for (typename ResPath::EdgeIt e(cycle); e != INVALID; ++e)
            _res_graph->augment(e, delta);

          // Finding the minimum cycle mean for the modified residual
          // graph
          mmc.reset();
          if (!mmc.findMinMean()) break;
        }
      }

      // Computing node potentials
      BellmanFord<ResGraph, ResidualCostMap> bf(*_res_graph, _res_cost);
      bf.init(0); bf.start();
      for (NodeIt n(_graph); n != INVALID; ++n)
        (*_potential)[n] = bf.dist(n);
    }

  }; //class CycleCanceling

  ///@}

} //namespace lemon

#endif //LEMON_CYCLE_CANCELING_H
