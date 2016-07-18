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

#ifndef LEMON_CAPACITY_SCALING_H
#define LEMON_CAPACITY_SCALING_H

/// \ingroup min_cost_flow
///
/// \file
/// \brief Capacity scaling algorithm for finding a minimum cost flow.

#include <vector>
#include <lemon/bin_heap.h>

namespace lemon {

  /// \addtogroup min_cost_flow
  /// @{

  /// \brief Implementation of the capacity scaling algorithm for
  /// finding a minimum cost flow.
  ///
  /// \ref CapacityScaling implements the capacity scaling version
  /// of the successive shortest path algorithm for finding a minimum
  /// cost flow.
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
  class CapacityScaling
  {
    GRAPH_TYPEDEFS(typename Graph);

    typedef typename CapacityMap::Value Capacity;
    typedef typename CostMap::Value Cost;
    typedef typename SupplyMap::Value Supply;
    typedef typename Graph::template EdgeMap<Capacity> CapacityEdgeMap;
    typedef typename Graph::template NodeMap<Supply> SupplyNodeMap;
    typedef typename Graph::template NodeMap<Edge> PredMap;

  public:

    /// The type of the flow map.
    typedef typename Graph::template EdgeMap<Capacity> FlowMap;
    /// The type of the potential map.
    typedef typename Graph::template NodeMap<Cost> PotentialMap;

  private:

    /// \brief Special implementation of the \ref Dijkstra algorithm
    /// for finding shortest paths in the residual network.
    ///
    /// \ref ResidualDijkstra is a special implementation of the
    /// \ref Dijkstra algorithm for finding shortest paths in the
    /// residual network of the graph with respect to the reduced edge
    /// costs and modifying the node potentials according to the
    /// distance of the nodes.
    class ResidualDijkstra
    {
      typedef typename Graph::template NodeMap<int> HeapCrossRef;
      typedef BinHeap<Cost, HeapCrossRef> Heap;

    private:

      // The directed graph the algorithm runs on
      const Graph &_graph;

      // The main maps
      const FlowMap &_flow;
      const CapacityEdgeMap &_res_cap;
      const CostMap &_cost;
      const SupplyNodeMap &_excess;
      PotentialMap &_potential;

      // The distance map
      PotentialMap _dist;
      // The pred edge map
      PredMap &_pred;
      // The processed (i.e. permanently labeled) nodes
      std::vector<Node> _proc_nodes;

    public:

      /// Constructor.
      ResidualDijkstra( const Graph &graph,
                        const FlowMap &flow,
                        const CapacityEdgeMap &res_cap,
                        const CostMap &cost,
                        const SupplyMap &excess,
                        PotentialMap &potential,
                        PredMap &pred ) :
        _graph(graph), _flow(flow), _res_cap(res_cap), _cost(cost),
        _excess(excess), _potential(potential), _dist(graph),
        _pred(pred)
      {}

      /// Runs the algorithm from the given source node.
      Node run(Node s, Capacity delta = 1) {
        HeapCrossRef heap_cross_ref(_graph, Heap::PRE_HEAP);
        Heap heap(heap_cross_ref);
        heap.push(s, 0);
        _pred[s] = INVALID;
        _proc_nodes.clear();

        // Processing nodes
        while (!heap.empty() && _excess[heap.top()] > -delta) {
          Node u = heap.top(), v;
          Cost d = heap.prio() + _potential[u], nd;
          _dist[u] = heap.prio();
          heap.pop();
          _proc_nodes.push_back(u);

          // Traversing outgoing edges
          for (OutEdgeIt e(_graph, u); e != INVALID; ++e) {
            if (_res_cap[e] >= delta) {
              v = _graph.target(e);
              switch(heap.state(v)) {
              case Heap::PRE_HEAP:
                heap.push(v, d + _cost[e] - _potential[v]);
                _pred[v] = e;
                break;
              case Heap::IN_HEAP:
                nd = d + _cost[e] - _potential[v];
                if (nd < heap[v]) {
                  heap.decrease(v, nd);
                  _pred[v] = e;
                }
                break;
              case Heap::POST_HEAP:
                break;
              }
            }
          }

          // Traversing incoming edges
          for (InEdgeIt e(_graph, u); e != INVALID; ++e) {
            if (_flow[e] >= delta) {
              v = _graph.source(e);
              switch(heap.state(v)) {
              case Heap::PRE_HEAP:
                heap.push(v, d - _cost[e] - _potential[v]);
                _pred[v] = e;
                break;
              case Heap::IN_HEAP:
                nd = d - _cost[e] - _potential[v];
                if (nd < heap[v]) {
                  heap.decrease(v, nd);
                  _pred[v] = e;
                }
                break;
              case Heap::POST_HEAP:
                break;
              }
            }
          }
        }
        if (heap.empty()) return INVALID;

        // Updating potentials of processed nodes
        Node t = heap.top();
        Cost t_dist = heap.prio();
        for (int i = 0; i < int(_proc_nodes.size()); ++i)
          _potential[_proc_nodes[i]] += _dist[_proc_nodes[i]] - t_dist;

        return t;
      }

    }; //class ResidualDijkstra

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

    // The residual capacity map
    CapacityEdgeMap _res_cap;
    // The excess map
    SupplyNodeMap _excess;
    // The excess nodes (i.e. nodes with positive excess)
    std::vector<Node> _excess_nodes;
    // The deficit nodes (i.e. nodes with negative excess)
    std::vector<Node> _deficit_nodes;

    // The delta parameter used for capacity scaling
    Capacity _delta;
    // The maximum number of phases
    int _phase_num;

    // The pred edge map
    PredMap _pred;
    // Implementation of the Dijkstra algorithm for finding augmenting
    // shortest paths in the residual network
    ResidualDijkstra *_dijkstra;

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
    CapacityScaling( const Graph &graph,
                     const LowerMap &lower,
                     const CapacityMap &capacity,
                     const CostMap &cost,
                     const SupplyMap &supply ) :
      _graph(graph), _lower(&lower), _capacity(graph), _cost(cost),
      _supply(graph), _flow(0), _local_flow(false),
      _potential(0), _local_potential(false),
      _res_cap(graph), _excess(graph), _pred(graph)
    {
      // Removing non-zero lower bounds
      _capacity = subMap(capacity, lower);
      _res_cap = _capacity;
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
    CapacityScaling( const Graph &graph,
                     const CapacityMap &capacity,
                     const CostMap &cost,
                     const SupplyMap &supply ) :
      _graph(graph), _lower(NULL), _capacity(capacity), _cost(cost),
      _supply(supply), _flow(0), _local_flow(false),
      _potential(0), _local_potential(false),
      _res_cap(capacity), _excess(graph), _pred(graph)
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
    CapacityScaling( const Graph &graph,
                     const LowerMap &lower,
                     const CapacityMap &capacity,
                     const CostMap &cost,
                     Node s, Node t,
                     Supply flow_value ) :
      _graph(graph), _lower(&lower), _capacity(graph), _cost(cost),
      _supply(graph), _flow(0), _local_flow(false),
      _potential(0), _local_potential(false),
      _res_cap(graph), _excess(graph), _pred(graph)
    {
      // Removing non-zero lower bounds
      _capacity = subMap(capacity, lower);
      _res_cap = _capacity;
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
    CapacityScaling( const Graph &graph,
                     const CapacityMap &capacity,
                     const CostMap &cost,
                     Node s, Node t,
                     Supply flow_value ) :
      _graph(graph), _lower(NULL), _capacity(capacity), _cost(cost),
      _supply(graph, 0), _flow(0), _local_flow(false),
      _potential(0), _local_potential(false),
      _res_cap(capacity), _excess(graph), _pred(graph)
    {
      _supply[s] =  flow_value;
      _supply[t] = -flow_value;
      _valid_supply = true;
    }

    /// Destructor.
    ~CapacityScaling() {
      if (_local_flow) delete _flow;
      if (_local_potential) delete _potential;
      delete _dijkstra;
    }

    /// \brief Sets the flow map.
    ///
    /// Sets the flow map.
    ///
    /// \return \c (*this)
    CapacityScaling& flowMap(FlowMap &map) {
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
    CapacityScaling& potentialMap(PotentialMap &map) {
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
    /// \param scaling Enable or disable capacity scaling.
    /// If the maximum edge capacity and/or the amount of total supply
    /// is rather small, the algorithm could be slightly faster without
    /// scaling.
    ///
    /// \return \c true if a feasible flow can be found.
    bool run(bool scaling = true) {
      return init(scaling) && start();
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
    bool init(bool scaling) {
      if (!_valid_supply) return false;

      // Initializing maps
      if (!_flow) {
        _flow = new FlowMap(_graph);
        _local_flow = true;
      }
      if (!_potential) {
        _potential = new PotentialMap(_graph);
        _local_potential = true;
      }
      for (EdgeIt e(_graph); e != INVALID; ++e) (*_flow)[e] = 0;
      for (NodeIt n(_graph); n != INVALID; ++n) (*_potential)[n] = 0;
      _excess = _supply;

      _dijkstra = new ResidualDijkstra( _graph, *_flow, _res_cap, _cost,
                                        _excess, *_potential, _pred );

      // Initializing delta value
      if (scaling) {
        // With scaling
        Supply max_sup = 0, max_dem = 0;
        for (NodeIt n(_graph); n != INVALID; ++n) {
          if ( _supply[n] > max_sup) max_sup =  _supply[n];
          if (-_supply[n] > max_dem) max_dem = -_supply[n];
        }
        Capacity max_cap = 0;
        for (EdgeIt e(_graph); e != INVALID; ++e) {
          if (_capacity[e] > max_cap) max_cap = _capacity[e];
        }
        max_sup = std::min(std::min(max_sup, max_dem), max_cap);
        _phase_num = 0;
        for (_delta = 1; 2 * _delta <= max_sup; _delta *= 2)
          ++_phase_num;
      } else {
        // Without scaling
        _delta = 1;
      }

      return true;
    }

    bool start() {
      if (_delta > 1)
        return startWithScaling();
      else
        return startWithoutScaling();
    }

    /// Executes the capacity scaling algorithm.
    bool startWithScaling() {
      // Processing capacity scaling phases
      Node s, t;
      int phase_cnt = 0;
      int factor = 4;
      while (true) {
        // Saturating all edges not satisfying the optimality condition
        for (EdgeIt e(_graph); e != INVALID; ++e) {
          Node u = _graph.source(e), v = _graph.target(e);
          Cost c = _cost[e] + (*_potential)[u] - (*_potential)[v];
          if (c < 0 && _res_cap[e] >= _delta) {
            _excess[u] -= _res_cap[e];
            _excess[v] += _res_cap[e];
            (*_flow)[e] = _capacity[e];
            _res_cap[e] = 0;
          }
          else if (c > 0 && (*_flow)[e] >= _delta) {
            _excess[u] += (*_flow)[e];
            _excess[v] -= (*_flow)[e];
            (*_flow)[e] = 0;
            _res_cap[e] = _capacity[e];
          }
        }

        // Finding excess nodes and deficit nodes
        _excess_nodes.clear();
        _deficit_nodes.clear();
        for (NodeIt n(_graph); n != INVALID; ++n) {
          if (_excess[n] >=  _delta) _excess_nodes.push_back(n);
          if (_excess[n] <= -_delta) _deficit_nodes.push_back(n);
        }
        int next_node = 0;

        // Finding augmenting shortest paths
        while (next_node < int(_excess_nodes.size())) {
          // Checking deficit nodes
          if (_delta > 1) {
            bool delta_deficit = false;
            for (int i = 0; i < int(_deficit_nodes.size()); ++i) {
              if (_excess[_deficit_nodes[i]] <= -_delta) {
                delta_deficit = true;
                break;
              }
            }
            if (!delta_deficit) break;
          }

          // Running Dijkstra
          s = _excess_nodes[next_node];
          if ((t = _dijkstra->run(s, _delta)) == INVALID) {
            if (_delta > 1) {
              ++next_node;
              continue;
            }
            return false;
          }

          // Augmenting along a shortest path from s to t.
          Capacity d = std::min(_excess[s], -_excess[t]);
          Node u = t;
          Edge e;
          if (d > _delta) {
            while ((e = _pred[u]) != INVALID) {
              Capacity rc;
              if (u == _graph.target(e)) {
                rc = _res_cap[e];
                u = _graph.source(e);
              } else {
                rc = (*_flow)[e];
                u = _graph.target(e);
              }
              if (rc < d) d = rc;
            }
          }
          u = t;
          while ((e = _pred[u]) != INVALID) {
            if (u == _graph.target(e)) {
              (*_flow)[e] += d;
              _res_cap[e] -= d;
              u = _graph.source(e);
            } else {
              (*_flow)[e] -= d;
              _res_cap[e] += d;
              u = _graph.target(e);
            }
          }
          _excess[s] -= d;
          _excess[t] += d;

          if (_excess[s] < _delta) ++next_node;
        }

        if (_delta == 1) break;
        if (++phase_cnt > _phase_num / 4) factor = 2;
        _delta = _delta <= factor ? 1 : _delta / factor;
      }

      // Handling non-zero lower bounds
      if (_lower) {
        for (EdgeIt e(_graph); e != INVALID; ++e)
          (*_flow)[e] += (*_lower)[e];
      }
      return true;
    }

    /// Executes the successive shortest path algorithm.
    bool startWithoutScaling() {
      // Finding excess nodes
      for (NodeIt n(_graph); n != INVALID; ++n)
        if (_excess[n] > 0) _excess_nodes.push_back(n);
      if (_excess_nodes.size() == 0) return true;
      int next_node = 0;

      // Finding shortest paths
      Node s, t;
      while ( _excess[_excess_nodes[next_node]] > 0 ||
              ++next_node < int(_excess_nodes.size()) )
      {
        // Running Dijkstra
        s = _excess_nodes[next_node];
        if ((t = _dijkstra->run(s)) == INVALID) return false;

        // Augmenting along a shortest path from s to t
        Capacity d = std::min(_excess[s], -_excess[t]);
        Node u = t;
        Edge e;
        if (d > 1) {
          while ((e = _pred[u]) != INVALID) {
            Capacity rc;
            if (u == _graph.target(e)) {
              rc = _res_cap[e];
              u = _graph.source(e);
            } else {
              rc = (*_flow)[e];
              u = _graph.target(e);
            }
            if (rc < d) d = rc;
          }
        }
        u = t;
        while ((e = _pred[u]) != INVALID) {
          if (u == _graph.target(e)) {
            (*_flow)[e] += d;
            _res_cap[e] -= d;
            u = _graph.source(e);
          } else {
            (*_flow)[e] -= d;
            _res_cap[e] += d;
            u = _graph.target(e);
          }
        }
        _excess[s] -= d;
        _excess[t] += d;
      }

      // Handling non-zero lower bounds
      if (_lower) {
        for (EdgeIt e(_graph); e != INVALID; ++e)
          (*_flow)[e] += (*_lower)[e];
      }
      return true;
    }

  }; //class CapacityScaling

  ///@}

} //namespace lemon

#endif //LEMON_CAPACITY_SCALING_H
