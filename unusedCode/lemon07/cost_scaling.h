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

#ifndef LEMON_COST_SCALING_H
#define LEMON_COST_SCALING_H

/// \ingroup min_cost_flow
///
/// \file
/// \brief Cost scaling algorithm for finding a minimum cost flow.

#include <deque>
#include <lemon/graph_adaptor.h>
#include <lemon/graph_utils.h>
#include <lemon/maps.h>
#include <lemon/math.h>

#include <lemon/circulation.h>
#include <lemon/bellman_ford.h>

namespace lemon {

  /// \addtogroup min_cost_flow
  /// @{

  /// \brief Implementation of the cost scaling algorithm for finding a
  /// minimum cost flow.
  ///
  /// \ref CostScaling implements the cost scaling algorithm performing
  /// generalized push-relabel operations for finding a minimum cost
  /// flow.
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
  /// \note Edge costs are multiplied with the number of nodes during
  /// the algorithm so overflow problems may arise more easily than with
  /// other minimum cost flow algorithms.
  /// If it is available, <tt>long long int</tt> type is used instead of
  /// <tt>long int</tt> in the inside computations.
  ///
  /// \author Peter Kovacs

  template < typename Graph,
             typename LowerMap = typename Graph::template EdgeMap<int>,
             typename CapacityMap = typename Graph::template EdgeMap<int>,
             typename CostMap = typename Graph::template EdgeMap<int>,
             typename SupplyMap = typename Graph::template NodeMap<int> >
  class CostScaling
  {
    GRAPH_TYPEDEFS(typename Graph);

    typedef typename CapacityMap::Value Capacity;
    typedef typename CostMap::Value Cost;
    typedef typename SupplyMap::Value Supply;
    typedef typename Graph::template EdgeMap<Capacity> CapacityEdgeMap;
    typedef typename Graph::template NodeMap<Supply> SupplyNodeMap;

    typedef ResGraphAdaptor< const Graph, Capacity,
                             CapacityEdgeMap, CapacityEdgeMap > ResGraph;
    typedef typename ResGraph::Edge ResEdge;

#if defined __GNUC__ && !defined __STRICT_ANSI__
    typedef long long int LCost;
#else
    typedef long int LCost;
#endif
    typedef typename Graph::template EdgeMap<LCost> LargeCostMap;

  public:

    /// The type of the flow map.
    typedef typename Graph::template EdgeMap<Capacity> FlowMap;
    /// The type of the potential map.
    typedef typename Graph::template NodeMap<LCost> PotentialMap;

  private:

    /// \brief Map adaptor class for handling residual edge costs.
    ///
    /// \ref ResidualCostMap is a map adaptor class for handling
    /// residual edge costs.
    template <typename Map>
    class ResidualCostMap : public MapBase<ResEdge, typename Map::Value>
    {
    private:

      const Map &_cost_map;

    public:

      ///\e
      ResidualCostMap(const Map &cost_map) :
        _cost_map(cost_map) {}

      ///\e
      typename Map::Value operator[](const ResEdge &e) const {
        return ResGraph::forward(e) ?  _cost_map[e] : -_cost_map[e];
      }

    }; //class ResidualCostMap

    /// \brief Map adaptor class for handling reduced edge costs.
    ///
    /// \ref ReducedCostMap is a map adaptor class for handling reduced
    /// edge costs.
    class ReducedCostMap : public MapBase<Edge, LCost>
    {
    private:

      const Graph &_gr;
      const LargeCostMap &_cost_map;
      const PotentialMap &_pot_map;

    public:

      ///\e
      ReducedCostMap( const Graph &gr,
                      const LargeCostMap &cost_map,
                      const PotentialMap &pot_map ) :
        _gr(gr), _cost_map(cost_map), _pot_map(pot_map) {}

      ///\e
      LCost operator[](const Edge &e) const {
        return _cost_map[e] + _pot_map[_gr.source(e)]
                            - _pot_map[_gr.target(e)];
      }

    }; //class ReducedCostMap

  private:

    // Scaling factor
    static const int ALPHA = 4;

    // Paramters for heuristics
    static const int BF_HEURISTIC_EPSILON_BOUND = 5000;
    static const int BF_HEURISTIC_BOUND_FACTOR  = 3;

  private:

    // The directed graph the algorithm runs on
    const Graph &_graph;
    // The original lower bound map
    const LowerMap *_lower;
    // The modified capacity map
    CapacityEdgeMap _capacity;
    // The original cost map
    const CostMap &_orig_cost;
    // The scaled cost map
    LargeCostMap _cost;
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
    ResidualCostMap<LargeCostMap> _res_cost;
    // The reduced cost map
    ReducedCostMap *_red_cost;
    // The excess map
    SupplyNodeMap _excess;
    // The epsilon parameter used for cost scaling
    LCost _epsilon;

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
    CostScaling( const Graph &graph,
                 const LowerMap &lower,
                 const CapacityMap &capacity,
                 const CostMap &cost,
                 const SupplyMap &supply ) :
      _graph(graph), _lower(&lower), _capacity(graph), _orig_cost(cost),
      _cost(graph), _supply(graph), _flow(0), _local_flow(false),
      _potential(0), _local_potential(false), _res_cost(_cost),
      _excess(graph, 0)
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
        _supply[n] = s;
        sum += s;
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
    CostScaling( const Graph &graph,
                 const CapacityMap &capacity,
                 const CostMap &cost,
                 const SupplyMap &supply ) :
      _graph(graph), _lower(NULL), _capacity(capacity), _orig_cost(cost),
      _cost(graph), _supply(supply), _flow(0), _local_flow(false),
      _potential(0), _local_potential(false), _res_cost(_cost),
      _excess(graph, 0)
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
    CostScaling( const Graph &graph,
                 const LowerMap &lower,
                 const CapacityMap &capacity,
                 const CostMap &cost,
                 Node s, Node t,
                 Supply flow_value ) :
      _graph(graph), _lower(&lower), _capacity(graph), _orig_cost(cost),
      _cost(graph), _supply(graph), _flow(0), _local_flow(false),
      _potential(0), _local_potential(false), _res_cost(_cost),
      _excess(graph, 0)
    {
      // Removing nonzero lower bounds
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
    CostScaling( const Graph &graph,
                 const CapacityMap &capacity,
                 const CostMap &cost,
                 Node s, Node t,
                 Supply flow_value ) :
      _graph(graph), _lower(NULL), _capacity(capacity), _orig_cost(cost),
      _cost(graph), _supply(graph, 0), _flow(0), _local_flow(false),
      _potential(0), _local_potential(false), _res_cost(_cost),
      _excess(graph, 0)
    {
      _supply[s] =  flow_value;
      _supply[t] = -flow_value;
      _valid_supply = true;
    }

    /// Destructor.
    ~CostScaling() {
      if (_local_flow) delete _flow;
      if (_local_potential) delete _potential;
      delete _res_graph;
      delete _red_cost;
    }

    /// \brief Sets the flow map.
    ///
    /// Sets the flow map.
    ///
    /// \return \c (*this)
    CostScaling& flowMap(FlowMap &map) {
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
    CostScaling& potentialMap(PotentialMap &map) {
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
    /// \return \c true if a feasible flow can be found.
    bool run() {
      return init() && start();
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
        c += (*_flow)[e] * _orig_cost[e];
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

      _red_cost = new ReducedCostMap(_graph, _cost, *_potential);
      _res_graph = new ResGraph(_graph, _capacity, *_flow);

      // Initializing the scaled cost map and the epsilon parameter
      Cost max_cost = 0;
      int node_num = countNodes(_graph);
      for (EdgeIt e(_graph); e != INVALID; ++e) {
        _cost[e] = LCost(_orig_cost[e]) * node_num * ALPHA;
        if (_orig_cost[e] > max_cost) max_cost = _orig_cost[e];
      }
      _epsilon = max_cost * node_num;

      // Finding a feasible flow using Circulation
      Circulation< Graph, ConstMap<Edge, Capacity>, CapacityEdgeMap,
                   SupplyMap >
        circulation( _graph, constMap<Edge>(Capacity(0)), _capacity,
                     _supply );
      return circulation.flowMap(*_flow).run();
    }


    /// Executes the algorithm.
    bool start() {
      std::deque<Node> active_nodes;
      typename Graph::template NodeMap<bool> hyper(_graph, false);

      int node_num = countNodes(_graph);
      for ( ; _epsilon >= 1; _epsilon = _epsilon < ALPHA && _epsilon > 1 ?
                                        1 : _epsilon / ALPHA )
      {
        // Performing price refinement heuristic using Bellman-Ford
        // algorithm
        if (_epsilon <= BF_HEURISTIC_EPSILON_BOUND) {
          typedef ShiftMap< ResidualCostMap<LargeCostMap> > ShiftCostMap;
          ShiftCostMap shift_cost(_res_cost, _epsilon);
          BellmanFord<ResGraph, ShiftCostMap> bf(*_res_graph, shift_cost);
          bf.init(0);
          bool done = false;
          int K = int(BF_HEURISTIC_BOUND_FACTOR * sqrt(node_num));
          for (int i = 0; i < K && !done; ++i)
            done = bf.processNextWeakRound();
          if (done) {
            for (NodeIt n(_graph); n != INVALID; ++n)
              (*_potential)[n] = bf.dist(n);
            continue;
          }
        }

        // Saturating edges not satisfying the optimality condition
        Capacity delta;
        for (EdgeIt e(_graph); e != INVALID; ++e) {
          if (_capacity[e] - (*_flow)[e] > 0 && (*_red_cost)[e] < 0) {
            delta = _capacity[e] - (*_flow)[e];
            _excess[_graph.source(e)] -= delta;
            _excess[_graph.target(e)] += delta;
            (*_flow)[e] = _capacity[e];
          }
          if ((*_flow)[e] > 0 && -(*_red_cost)[e] < 0) {
            _excess[_graph.target(e)] -= (*_flow)[e];
            _excess[_graph.source(e)] += (*_flow)[e];
            (*_flow)[e] = 0;
          }
        }

        // Finding active nodes (i.e. nodes with positive excess)
        for (NodeIt n(_graph); n != INVALID; ++n)
          if (_excess[n] > 0) active_nodes.push_back(n);

        // Performing push and relabel operations
        while (active_nodes.size() > 0) {
          Node n = active_nodes[0], t;
          bool relabel_enabled = true;

          // Performing push operations if there are admissible edges
          if (_excess[n] > 0) {
            for (OutEdgeIt e(_graph, n); e != INVALID; ++e) {
              if (_capacity[e] - (*_flow)[e] > 0 && (*_red_cost)[e] < 0) {
                delta = _capacity[e] - (*_flow)[e] <= _excess[n] ?
                        _capacity[e] - (*_flow)[e] : _excess[n];
                t = _graph.target(e);

                // Push-look-ahead heuristic
                Capacity ahead = -_excess[t];
                for (OutEdgeIt oe(_graph, t); oe != INVALID; ++oe) {
                  if (_capacity[oe] - (*_flow)[oe] > 0 && (*_red_cost)[oe] < 0)
                    ahead += _capacity[oe] - (*_flow)[oe];
                }
                for (InEdgeIt ie(_graph, t); ie != INVALID; ++ie) {
                  if ((*_flow)[ie] > 0 && -(*_red_cost)[ie] < 0)
                    ahead += (*_flow)[ie];
                }
                if (ahead < 0) ahead = 0;

                // Pushing flow along the edge
                if (ahead < delta) {
                  (*_flow)[e] += ahead;
                  _excess[n] -= ahead;
                  _excess[t] += ahead;
                  active_nodes.push_front(t);
                  hyper[t] = true;
                  relabel_enabled = false;
                  break;
                } else {
                  (*_flow)[e] += delta;
                  _excess[n] -= delta;
                  _excess[t] += delta;
                  if (_excess[t] > 0 && _excess[t] <= delta)
                    active_nodes.push_back(t);
                }

                if (_excess[n] == 0) break;
              }
            }
          }

          if (_excess[n] > 0) {
            for (InEdgeIt e(_graph, n); e != INVALID; ++e) {
              if ((*_flow)[e] > 0 && -(*_red_cost)[e] < 0) {
                delta = (*_flow)[e] <= _excess[n] ? (*_flow)[e] : _excess[n];
                t = _graph.source(e);

                // Push-look-ahead heuristic
                Capacity ahead = -_excess[t];
                for (OutEdgeIt oe(_graph, t); oe != INVALID; ++oe) {
                  if (_capacity[oe] - (*_flow)[oe] > 0 && (*_red_cost)[oe] < 0)
                    ahead += _capacity[oe] - (*_flow)[oe];
                }
                for (InEdgeIt ie(_graph, t); ie != INVALID; ++ie) {
                  if ((*_flow)[ie] > 0 && -(*_red_cost)[ie] < 0)
                    ahead += (*_flow)[ie];
                }
                if (ahead < 0) ahead = 0;

                // Pushing flow along the edge
                if (ahead < delta) {
                  (*_flow)[e] -= ahead;
                  _excess[n] -= ahead;
                  _excess[t] += ahead;
                  active_nodes.push_front(t);
                  hyper[t] = true;
                  relabel_enabled = false;
                  break;
                } else {
                  (*_flow)[e] -= delta;
                  _excess[n] -= delta;
                  _excess[t] += delta;
                  if (_excess[t] > 0 && _excess[t] <= delta)
                    active_nodes.push_back(t);
                }

                if (_excess[n] == 0) break;
              }
            }
          }

          if (relabel_enabled && (_excess[n] > 0 || hyper[n])) {
            // Performing relabel operation if the node is still active
            LCost min_red_cost = std::numeric_limits<LCost>::max();
            for (OutEdgeIt oe(_graph, n); oe != INVALID; ++oe) {
              if ( _capacity[oe] - (*_flow)[oe] > 0 &&
                   (*_red_cost)[oe] < min_red_cost )
                min_red_cost = (*_red_cost)[oe];
            }
            for (InEdgeIt ie(_graph, n); ie != INVALID; ++ie) {
              if ((*_flow)[ie] > 0 && -(*_red_cost)[ie] < min_red_cost)
                min_red_cost = -(*_red_cost)[ie];
            }
            (*_potential)[n] -= min_red_cost + _epsilon;
            hyper[n] = false;
          }

          // Removing active nodes with non-positive excess
          while ( active_nodes.size() > 0 &&
                  _excess[active_nodes[0]] <= 0 &&
                  !hyper[active_nodes[0]] ) {
            active_nodes.pop_front();
          }
        }
      }

      // Computing node potentials for the original costs
      ResidualCostMap<CostMap> res_cost(_orig_cost);
      BellmanFord< ResGraph, ResidualCostMap<CostMap> >
        bf(*_res_graph, res_cost);
      bf.init(0); bf.start();
      for (NodeIt n(_graph); n != INVALID; ++n)
        (*_potential)[n] = bf.dist(n);

      // Handling non-zero lower bounds
      if (_lower) {
        for (EdgeIt e(_graph); e != INVALID; ++e)
          (*_flow)[e] += (*_lower)[e];
      }
      return true;
    }

  }; //class CostScaling

  ///@}

} //namespace lemon

#endif //LEMON_COST_SCALING_H
