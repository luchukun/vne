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

#ifndef LEMON_NETWORK_SIMPLEX_H
#define LEMON_NETWORK_SIMPLEX_H

/// \ingroup min_cost_flow
///
/// \file
/// \brief Network simplex algorithm for finding a minimum cost flow.

#include <vector>
#include <limits>

#include <lemon/graph_adaptor.h>
#include <lemon/graph_utils.h>
#include <lemon/smart_graph.h>
#include <lemon/math.h>

namespace lemon {

  /// \addtogroup min_cost_flow
  /// @{

  /// \brief Implementation of the network simplex algorithm for
  /// finding a minimum cost flow.
  ///
  /// \ref NetworkSimplex implements the network simplex algorithm for
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
  /// \note \ref NetworkSimplex provides six different pivot rule
  /// implementations that significantly affect the efficiency of the
  /// algorithm.
  /// By default a combined pivot rule is used, which is the fastest
  /// implementation according to our benchmark tests.
  /// Another pivot rule can be selected using \ref run() function
  /// with the proper parameter.
  ///
  /// \author Peter Kovacs

  template < typename Graph,
             typename LowerMap = typename Graph::template EdgeMap<int>,
             typename CapacityMap = typename Graph::template EdgeMap<int>,
             typename CostMap = typename Graph::template EdgeMap<int>,
             typename SupplyMap = typename Graph::template NodeMap<int> >
  class NetworkSimplex
  {
    typedef typename CapacityMap::Value Capacity;
    typedef typename CostMap::Value Cost;
    typedef typename SupplyMap::Value Supply;

    typedef SmartGraph SGraph;
    GRAPH_TYPEDEFS(typename SGraph);

    typedef typename SGraph::template EdgeMap<Capacity> SCapacityMap;
    typedef typename SGraph::template EdgeMap<Cost> SCostMap;
    typedef typename SGraph::template NodeMap<Supply> SSupplyMap;
    typedef typename SGraph::template NodeMap<Cost> SPotentialMap;

    typedef typename SGraph::template NodeMap<int> IntNodeMap;
    typedef typename SGraph::template NodeMap<bool> BoolNodeMap;
    typedef typename SGraph::template NodeMap<Node> NodeNodeMap;
    typedef typename SGraph::template NodeMap<Edge> EdgeNodeMap;
    typedef typename SGraph::template EdgeMap<int> IntEdgeMap;

    typedef typename Graph::template NodeMap<Node> NodeRefMap;
    typedef typename Graph::template EdgeMap<Edge> EdgeRefMap;

  public:

    /// The type of the flow map.
    typedef typename Graph::template EdgeMap<Capacity> FlowMap;
    /// The type of the potential map.
    typedef typename Graph::template NodeMap<Cost> PotentialMap;

  public:

    /// Enum type to select the pivot rule used by \ref run().
    enum PivotRuleEnum {
      FIRST_ELIGIBLE_PIVOT,
      BEST_ELIGIBLE_PIVOT,
      BLOCK_SEARCH_PIVOT,
      LIMITED_SEARCH_PIVOT,
      CANDIDATE_LIST_PIVOT,
      COMBINED_PIVOT
    };

  private:

    /// \brief Map adaptor class for handling reduced edge costs.
    ///
    /// Map adaptor class for handling reduced edge costs.
    class ReducedCostMap : public MapBase<Edge, Cost>
    {
    private:

      const SGraph &_gr;
      const SCostMap &_cost_map;
      const SPotentialMap &_pot_map;

    public:

      ///\e
      ReducedCostMap( const SGraph &gr,
                      const SCostMap &cost_map,
                      const SPotentialMap &pot_map ) :
        _gr(gr), _cost_map(cost_map), _pot_map(pot_map) {}

      ///\e
      Cost operator[](const Edge &e) const {
        return _cost_map[e] + _pot_map[_gr.source(e)]
                           - _pot_map[_gr.target(e)];
      }

    }; //class ReducedCostMap

  private:

    /// \brief Implementation of the "First Eligible" pivot rule for the
    /// \ref NetworkSimplex "network simplex" algorithm.
    ///
    /// This class implements the "First Eligible" pivot rule
    /// for the \ref NetworkSimplex "network simplex" algorithm.
    class FirstEligiblePivotRule
    {
    private:

      NetworkSimplex &_ns;
      EdgeIt _next_edge;

    public:

      /// Constructor.
      FirstEligiblePivotRule(NetworkSimplex &ns) :
        _ns(ns), _next_edge(ns._graph) {}

      /// Finds the next entering edge.
      bool findEnteringEdge() {
        for (EdgeIt e = _next_edge; e != INVALID; ++e) {
          if (_ns._state[e] * _ns._red_cost[e] < 0) {
            _ns._in_edge = e;
            _next_edge = ++e;
            return true;
          }
        }
        for (EdgeIt e(_ns._graph); e != _next_edge; ++e) {
          if (_ns._state[e] * _ns._red_cost[e] < 0) {
            _ns._in_edge = e;
            _next_edge = ++e;
            return true;
          }
        }
        return false;
      }
    }; //class FirstEligiblePivotRule

    /// \brief Implementation of the "Best Eligible" pivot rule for the
    /// \ref NetworkSimplex "network simplex" algorithm.
    ///
    /// This class implements the "Best Eligible" pivot rule
    /// for the \ref NetworkSimplex "network simplex" algorithm.
    class BestEligiblePivotRule
    {
    private:

      NetworkSimplex &_ns;

    public:

      /// Constructor.
      BestEligiblePivotRule(NetworkSimplex &ns) : _ns(ns) {}

      /// Finds the next entering edge.
      bool findEnteringEdge() {
        Cost min = 0;
        for (EdgeIt e(_ns._graph); e != INVALID; ++e) {
          if (_ns._state[e] * _ns._red_cost[e] < min) {
            min = _ns._state[e] * _ns._red_cost[e];
            _ns._in_edge = e;
          }
        }
        return min < 0;
      }
    }; //class BestEligiblePivotRule

    /// \brief Implementation of the "Block Search" pivot rule for the
    /// \ref NetworkSimplex "network simplex" algorithm.
    ///
    /// This class implements the "Block Search" pivot rule
    /// for the \ref NetworkSimplex "network simplex" algorithm.
    class BlockSearchPivotRule
    {
    private:

      NetworkSimplex &_ns;
      EdgeIt _next_edge, _min_edge;
      int _block_size;

      static const int MIN_BLOCK_SIZE = 10;

    public:

      /// Constructor.
      BlockSearchPivotRule(NetworkSimplex &ns) :
        _ns(ns), _next_edge(ns._graph), _min_edge(ns._graph)
      {
        _block_size = 2 * int(sqrt(countEdges(_ns._graph)));
        if (_block_size < MIN_BLOCK_SIZE) _block_size = MIN_BLOCK_SIZE;
      }

      /// Finds the next entering edge.
      bool findEnteringEdge() {
        Cost curr, min = 0;
        int cnt = 0;
        for (EdgeIt e = _next_edge; e != INVALID; ++e) {
          if ((curr = _ns._state[e] * _ns._red_cost[e]) < min) {
            min = curr;
            _min_edge = e;
          }
          if (++cnt == _block_size) {
            if (min < 0) break;
            cnt = 0;
          }
        }
        if (min == 0) {
          for (EdgeIt e(_ns._graph); e != _next_edge; ++e) {
            if ((curr = _ns._state[e] * _ns._red_cost[e]) < min) {
              min = curr;
              _min_edge = e;
            }
            if (++cnt == _block_size) {
              if (min < 0) break;
              cnt = 0;
            }
          }
        }
        _ns._in_edge = _min_edge;
        _next_edge = ++_min_edge;
        return min < 0;
      }
    }; //class BlockSearchPivotRule

    /// \brief Implementation of the "Limited Search" pivot rule for the
    /// \ref NetworkSimplex "network simplex" algorithm.
    ///
    /// This class implements the "Limited Search" pivot rule
    /// for the \ref NetworkSimplex "network simplex" algorithm.
    class LimitedSearchPivotRule
    {
    private:

      NetworkSimplex &_ns;
      EdgeIt _next_edge, _min_edge;
      int _sample_size;

      static const int SAMPLE_SIZE_FACTOR = 15;
      static const int MIN_SAMPLE_SIZE = 10;

    public:

      /// Constructor.
      LimitedSearchPivotRule(NetworkSimplex &ns) :
        _ns(ns), _next_edge(ns._graph), _min_edge(ns._graph)
      {
        _sample_size = countEdges(_ns._graph) *
                       SAMPLE_SIZE_FACTOR / 10000;
        if (_sample_size < MIN_SAMPLE_SIZE)
          _sample_size = MIN_SAMPLE_SIZE;
      }

      /// Finds the next entering edge.
      bool findEnteringEdge() {
        Cost curr, min = 0;
        int cnt = 0;
        for (EdgeIt e = _next_edge; e != INVALID; ++e) {
          if ((curr = _ns._state[e] * _ns._red_cost[e]) < min) {
            min = curr;
            _min_edge = e;
          }
          if (curr < 0 && ++cnt == _sample_size) break;
        }
        if (min == 0) {
          for (EdgeIt e(_ns._graph); e != _next_edge; ++e) {
            if ((curr = _ns._state[e] * _ns._red_cost[e]) < min) {
              min = curr;
              _min_edge = e;
            }
            if (curr < 0 && ++cnt == _sample_size) break;
          }
        }
        _ns._in_edge = _min_edge;
        _next_edge = ++_min_edge;
        return min < 0;
      }
    }; //class LimitedSearchPivotRule

    /// \brief Implementation of the "Candidate List" pivot rule for the
    /// \ref NetworkSimplex "network simplex" algorithm.
    ///
    /// This class implements the "Candidate List" pivot rule
    /// for the \ref NetworkSimplex "network simplex" algorithm.
    class CandidateListPivotRule
    {
    private:

      NetworkSimplex &_ns;

      // The list of candidate edges.
      std::vector<Edge> _candidates;
      // The maximum length of the edge list.
      int _list_length;
      // The maximum number of minor iterations between two major
      // itarations.
      int _minor_limit;

      int _minor_count;
      EdgeIt _next_edge;

      static const int LIST_LENGTH_FACTOR = 20;
      static const int MINOR_LIMIT_FACTOR = 10;
      static const int MIN_LIST_LENGTH = 10;
      static const int MIN_MINOR_LIMIT = 2;

    public:

      /// Constructor.
      CandidateListPivotRule(NetworkSimplex &ns) :
        _ns(ns), _next_edge(ns._graph)
      {
        int edge_num = countEdges(_ns._graph);
        _minor_count = 0;
        _list_length = edge_num * LIST_LENGTH_FACTOR / 10000;
        if (_list_length < MIN_LIST_LENGTH)
          _list_length = MIN_LIST_LENGTH;
        _minor_limit = _list_length * MINOR_LIMIT_FACTOR / 100;
        if (_minor_limit < MIN_MINOR_LIMIT)
          _minor_limit = MIN_MINOR_LIMIT;
      }

      /// Finds the next entering edge.
      bool findEnteringEdge() {
        Cost min, curr;
        if (_minor_count < _minor_limit && _candidates.size() > 0) {
          // Minor iteration
          ++_minor_count;
          Edge e;
          min = 0;
          for (int i = 0; i < int(_candidates.size()); ++i) {
            e = _candidates[i];
            if ((curr = _ns._state[e] * _ns._red_cost[e]) < min) {
              min = curr;
              _ns._in_edge = e;
            }
          }
          if (min < 0) return true;
        }

        // Major iteration
        _candidates.clear();
        EdgeIt e = _next_edge;
        min = 0;
        for ( ; e != INVALID; ++e) {
          if ((curr = _ns._state[e] * _ns._red_cost[e]) < 0) {
            _candidates.push_back(e);
            if (curr < min) {
              min = curr;
              _ns._in_edge = e;
            }
            if (int(_candidates.size()) == _list_length) break;
          }
        }
        if (int(_candidates.size()) < _list_length) {
          for (e = EdgeIt(_ns._graph); e != _next_edge; ++e) {
            if ((curr = _ns._state[e] * _ns._red_cost[e]) < 0) {
              _candidates.push_back(e);
              if (curr < min) {
                min = curr;
                _ns._in_edge = e;
              }
              if (int(_candidates.size()) == _list_length) break;
            }
          }
        }
        if (_candidates.size() == 0) return false;
        _minor_count = 1;
        _next_edge = ++e;
        return true;
      }
    }; //class CandidateListPivotRule

  private:

    // State constants for edges
    enum EdgeStateEnum {
      STATE_UPPER = -1,
      STATE_TREE  =  0,
      STATE_LOWER =  1
    };

    // Constant for the combined pivot rule.
    static const int COMBINED_PIVOT_MAX_DEG = 5;

  private:

    // The directed graph the algorithm runs on
    SGraph _graph;
    // The original graph
    const Graph &_graph_ref;
    // The original lower bound map
    const LowerMap *_lower;
    // The capacity map
    SCapacityMap _capacity;
    // The cost map
    SCostMap _cost;
    // The supply map
    SSupplyMap _supply;
    bool _valid_supply;

    // Edge map of the current flow
    SCapacityMap _flow;
    // Node map of the current potentials
    SPotentialMap _potential;

    // The depth node map of the spanning tree structure
    IntNodeMap _depth;
    // The parent node map of the spanning tree structure
    NodeNodeMap _parent;
    // The pred_edge node map of the spanning tree structure
    EdgeNodeMap _pred_edge;
    // The thread node map of the spanning tree structure
    NodeNodeMap _thread;
    // The forward node map of the spanning tree structure
    BoolNodeMap _forward;
    // The state edge map
    IntEdgeMap _state;
    // The root node of the starting spanning tree
    Node _root;

    // The reduced cost map
    ReducedCostMap _red_cost;

    // Members for handling the original graph
    FlowMap *_flow_result;
    PotentialMap *_potential_result;
    bool _local_flow;
    bool _local_potential;
    NodeRefMap _node_ref;
    EdgeRefMap _edge_ref;

    // The entering edge of the current pivot iteration.
    Edge _in_edge;

    // Temporary nodes used in the current pivot iteration.
    Node join, u_in, v_in, u_out, v_out;
    Node right, first, second, last;
    Node stem, par_stem, new_stem;
    // The maximum augment amount along the found cycle in the current
    // pivot iteration.
    Capacity delta;

  public :

    /// \brief General constructor (with lower bounds).
    ///
    /// General constructor (with lower bounds).
    ///
    /// \param graph The directed graph the algorithm runs on.
    /// \param lower The lower bounds of the edges.
    /// \param capacity The capacities (upper bounds) of the edges.
    /// \param cost The cost (length) values of the edges.
    /// \param supply The supply values of the nodes (signed).
    NetworkSimplex( const Graph &graph,
                    const LowerMap &lower,
                    const CapacityMap &capacity,
                    const CostMap &cost,
                    const SupplyMap &supply ) :
      _graph(), _graph_ref(graph), _lower(&lower), _capacity(_graph),
      _cost(_graph), _supply(_graph), _flow(_graph),
      _potential(_graph), _depth(_graph), _parent(_graph),
      _pred_edge(_graph), _thread(_graph), _forward(_graph),
      _state(_graph), _red_cost(_graph, _cost, _potential),
      _flow_result(0), _potential_result(0),
      _local_flow(false), _local_potential(false),
      _node_ref(graph), _edge_ref(graph)
    {
      // Checking the sum of supply values
      Supply sum = 0;
      for (typename Graph::NodeIt n(_graph_ref); n != INVALID; ++n)
        sum += supply[n];
      if (!(_valid_supply = sum == 0)) return;

      // Copying _graph_ref to _graph
      _graph.reserveNode(countNodes(_graph_ref) + 1);
      _graph.reserveEdge(countEdges(_graph_ref) + countNodes(_graph_ref));
      copyGraph(_graph, _graph_ref)
        .edgeMap(_cost, cost)
        .nodeRef(_node_ref)
        .edgeRef(_edge_ref)
        .run();

      // Removing non-zero lower bounds
      for (typename Graph::EdgeIt e(_graph_ref); e != INVALID; ++e) {
        _capacity[_edge_ref[e]] = capacity[e] - lower[e];
      }
      for (typename Graph::NodeIt n(_graph_ref); n != INVALID; ++n) {
        Supply s = supply[n];
        for (typename Graph::InEdgeIt e(_graph_ref, n); e != INVALID; ++e)
          s += lower[e];
        for (typename Graph::OutEdgeIt e(_graph_ref, n); e != INVALID; ++e)
          s -= lower[e];
        _supply[_node_ref[n]] = s;
      }
    }

    /// \brief General constructor (without lower bounds).
    ///
    /// General constructor (without lower bounds).
    ///
    /// \param graph The directed graph the algorithm runs on.
    /// \param capacity The capacities (upper bounds) of the edges.
    /// \param cost The cost (length) values of the edges.
    /// \param supply The supply values of the nodes (signed).
    NetworkSimplex( const Graph &graph,
                    const CapacityMap &capacity,
                    const CostMap &cost,
                    const SupplyMap &supply ) :
      _graph(), _graph_ref(graph), _lower(NULL), _capacity(_graph),
      _cost(_graph), _supply(_graph), _flow(_graph),
      _potential(_graph), _depth(_graph), _parent(_graph),
      _pred_edge(_graph), _thread(_graph), _forward(_graph),
      _state(_graph), _red_cost(_graph, _cost, _potential),
      _flow_result(0), _potential_result(0),
      _local_flow(false), _local_potential(false),
      _node_ref(graph), _edge_ref(graph)
    {
      // Checking the sum of supply values
      Supply sum = 0;
      for (typename Graph::NodeIt n(_graph_ref); n != INVALID; ++n)
        sum += supply[n];
      if (!(_valid_supply = sum == 0)) return;

      // Copying _graph_ref to graph
      copyGraph(_graph, _graph_ref)
        .edgeMap(_capacity, capacity)
        .edgeMap(_cost, cost)
        .nodeMap(_supply, supply)
        .nodeRef(_node_ref)
        .edgeRef(_edge_ref)
        .run();
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
    NetworkSimplex( const Graph &graph,
                    const LowerMap &lower,
                    const CapacityMap &capacity,
                    const CostMap &cost,
                    typename Graph::Node s,
                    typename Graph::Node t,
                    typename SupplyMap::Value flow_value ) :
      _graph(), _graph_ref(graph), _lower(&lower), _capacity(_graph),
      _cost(_graph), _supply(_graph), _flow(_graph),
      _potential(_graph), _depth(_graph), _parent(_graph),
      _pred_edge(_graph), _thread(_graph), _forward(_graph),
      _state(_graph), _red_cost(_graph, _cost, _potential),
      _flow_result(0), _potential_result(0),
      _local_flow(false), _local_potential(false),
      _node_ref(graph), _edge_ref(graph)
    {
      // Copying _graph_ref to graph
      copyGraph(_graph, _graph_ref)
        .edgeMap(_cost, cost)
        .nodeRef(_node_ref)
        .edgeRef(_edge_ref)
        .run();

      // Removing non-zero lower bounds
      for (typename Graph::EdgeIt e(_graph_ref); e != INVALID; ++e) {
        _capacity[_edge_ref[e]] = capacity[e] - lower[e];
      }
      for (typename Graph::NodeIt n(_graph_ref); n != INVALID; ++n) {
        Supply sum = 0;
        if (n == s) sum =  flow_value;
        if (n == t) sum = -flow_value;
        for (typename Graph::InEdgeIt e(_graph_ref, n); e != INVALID; ++e)
          sum += lower[e];
        for (typename Graph::OutEdgeIt e(_graph_ref, n); e != INVALID; ++e)
          sum -= lower[e];
        _supply[_node_ref[n]] = sum;
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
    NetworkSimplex( const Graph &graph,
                    const CapacityMap &capacity,
                    const CostMap &cost,
                    typename Graph::Node s,
                    typename Graph::Node t,
                    typename SupplyMap::Value flow_value ) :
      _graph(), _graph_ref(graph), _lower(NULL), _capacity(_graph),
      _cost(_graph), _supply(_graph, 0), _flow(_graph),
      _potential(_graph), _depth(_graph), _parent(_graph),
      _pred_edge(_graph), _thread(_graph), _forward(_graph),
      _state(_graph), _red_cost(_graph, _cost, _potential),
      _flow_result(0), _potential_result(0),
      _local_flow(false), _local_potential(false),
      _node_ref(graph), _edge_ref(graph)
    {
      // Copying _graph_ref to graph
      copyGraph(_graph, _graph_ref)
        .edgeMap(_capacity, capacity)
        .edgeMap(_cost, cost)
        .nodeRef(_node_ref)
        .edgeRef(_edge_ref)
        .run();
      _supply[_node_ref[s]] =  flow_value;
      _supply[_node_ref[t]] = -flow_value;
      _valid_supply = true;
    }

    /// Destructor.
    ~NetworkSimplex() {
      if (_local_flow) delete _flow_result;
      if (_local_potential) delete _potential_result;
    }

    /// \brief Sets the flow map.
    ///
    /// Sets the flow map.
    ///
    /// \return \c (*this)
    NetworkSimplex& flowMap(FlowMap &map) {
      if (_local_flow) {
        delete _flow_result;
        _local_flow = false;
      }
      _flow_result = &map;
      return *this;
    }

    /// \brief Sets the potential map.
    ///
    /// Sets the potential map.
    ///
    /// \return \c (*this)
    NetworkSimplex& potentialMap(PotentialMap &map) {
      if (_local_potential) {
        delete _potential_result;
        _local_potential = false;
      }
      _potential_result = &map;
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
    /// \param pivot_rule The pivot rule that is used during the
    /// algorithm.
    ///
    /// The available pivot rules:
    ///
    /// - FIRST_ELIGIBLE_PIVOT The next eligible edge is selected in
    /// a wraparound fashion in every iteration
    /// (\ref FirstEligiblePivotRule).
    ///
    /// - BEST_ELIGIBLE_PIVOT The best eligible edge is selected in
    /// every iteration (\ref BestEligiblePivotRule).
    ///
    /// - BLOCK_SEARCH_PIVOT A specified number of edges are examined in
    /// every iteration in a wraparound fashion and the best eligible
    /// edge is selected from this block (\ref BlockSearchPivotRule).
    ///
    /// - LIMITED_SEARCH_PIVOT A specified number of eligible edges are
    /// examined in every iteration in a wraparound fashion and the best
    /// one is selected from them (\ref LimitedSearchPivotRule).
    ///
    /// - CANDIDATE_LIST_PIVOT In major iterations a candidate list is
    /// built from eligible edges and it is used for edge selection in
    /// the following minor iterations (\ref CandidateListPivotRule).
    ///
    /// - COMBINED_PIVOT This is a combined version of the two fastest
    /// pivot rules.
    /// For rather sparse graphs \ref LimitedSearchPivotRule
    /// "Limited Search" implementation is used, otherwise
    /// \ref BlockSearchPivotRule "Block Search" pivot rule is used.
    /// According to our benchmark tests this combined method is the
    /// most efficient.
    ///
    /// \return \c true if a feasible flow can be found.
    bool run(PivotRuleEnum pivot_rule = COMBINED_PIVOT) {
      return init() && start(pivot_rule);
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
      return *_flow_result;
    }

    /// \brief Returns a const reference to the node map storing the
    /// found potentials (the dual solution).
    ///
    /// Returns a const reference to the node map storing the found
    /// potentials (the dual solution).
    ///
    /// \pre \ref run() must be called before using this function.
    const PotentialMap& potentialMap() const {
      return *_potential_result;
    }

    /// \brief Returns the flow on the given edge.
    ///
    /// Returns the flow on the given edge.
    ///
    /// \pre \ref run() must be called before using this function.
    Capacity flow(const typename Graph::Edge& edge) const {
      return (*_flow_result)[edge];
    }

    /// \brief Returns the potential of the given node.
    ///
    /// Returns the potential of the given node.
    ///
    /// \pre \ref run() must be called before using this function.
    Cost potential(const typename Graph::Node& node) const {
      return (*_potential_result)[node];
    }

    /// \brief Returns the total cost of the found flow.
    ///
    /// Returns the total cost of the found flow. The complexity of the
    /// function is \f$ O(e) \f$.
    ///
    /// \pre \ref run() must be called before using this function.
    Cost totalCost() const {
      Cost c = 0;
      for (typename Graph::EdgeIt e(_graph_ref); e != INVALID; ++e)
        c += (*_flow_result)[e] * _cost[_edge_ref[e]];
      return c;
    }

    /// @}

  private:

    /// \brief Extends the underlaying graph and initializes all the
    /// node and edge maps.
    bool init() {
      if (!_valid_supply) return false;

      // Initializing result maps
      if (!_flow_result) {
        _flow_result = new FlowMap(_graph_ref);
        _local_flow = true;
      }
      if (!_potential_result) {
        _potential_result = new PotentialMap(_graph_ref);
        _local_potential = true;
      }

      // Initializing state and flow maps
      for (EdgeIt e(_graph); e != INVALID; ++e) {
        _flow[e] = 0;
        _state[e] = STATE_LOWER;
      }

      // Adding an artificial root node to the graph
      _root = _graph.addNode();
      _parent[_root] = INVALID;
      _pred_edge[_root] = INVALID;
      _depth[_root] = 0;
      _supply[_root] = 0;
      _potential[_root] = 0;

      // Adding artificial edges to the graph and initializing the node
      // maps of the spanning tree data structure
      Node last = _root;
      Edge e;
      Cost max_cost = std::numeric_limits<Cost>::max() / 4;
      for (NodeIt u(_graph); u != INVALID; ++u) {
        if (u == _root) continue;
        _thread[last] = u;
        last = u;
        _parent[u] = _root;
        _depth[u] = 1;
        if (_supply[u] >= 0) {
          e = _graph.addEdge(u, _root);
          _flow[e] = _supply[u];
          _forward[u] = true;
          _potential[u] = -max_cost;
        } else {
          e = _graph.addEdge(_root, u);
          _flow[e] = -_supply[u];
          _forward[u] = false;
          _potential[u] = max_cost;
        }
        _cost[e] = max_cost;
        _capacity[e] = std::numeric_limits<Capacity>::max();
        _state[e] = STATE_TREE;
        _pred_edge[u] = e;
      }
      _thread[last] = _root;

      return true;
    }

    /// Finds the join node.
    Node findJoinNode() {
      Node u = _graph.source(_in_edge);
      Node v = _graph.target(_in_edge);
      while (u != v) {
        if (_depth[u] == _depth[v]) {
          u = _parent[u];
          v = _parent[v];
        }
        else if (_depth[u] > _depth[v]) u = _parent[u];
        else v = _parent[v];
      }
      return u;
    }

    /// \brief Finds the leaving edge of the cycle. Returns \c true if
    /// the leaving edge is not the same as the entering edge.
    bool findLeavingEdge() {
      // Initializing first and second nodes according to the direction
      // of the cycle
      if (_state[_in_edge] == STATE_LOWER) {
        first  = _graph.source(_in_edge);
        second = _graph.target(_in_edge);
      } else {
        first  = _graph.target(_in_edge);
        second = _graph.source(_in_edge);
      }
      delta = _capacity[_in_edge];
      bool result = false;
      Capacity d;
      Edge e;

      // Searching the cycle along the path form the first node to the
      // root node
      for (Node u = first; u != join; u = _parent[u]) {
        e = _pred_edge[u];
        d = _forward[u] ? _flow[e] : _capacity[e] - _flow[e];
        if (d < delta) {
          delta = d;
          u_out = u;
          u_in = first;
          v_in = second;
          result = true;
        }
      }
      // Searching the cycle along the path form the second node to the
      // root node
      for (Node u = second; u != join; u = _parent[u]) {
        e = _pred_edge[u];
        d = _forward[u] ? _capacity[e] - _flow[e] : _flow[e];
        if (d <= delta) {
          delta = d;
          u_out = u;
          u_in = second;
          v_in = first;
          result = true;
        }
      }
      return result;
    }

    /// Changes \c flow and \c state edge maps.
    void changeFlows(bool change) {
      // Augmenting along the cycle
      if (delta > 0) {
        Capacity val = _state[_in_edge] * delta;
        _flow[_in_edge] += val;
        for (Node u = _graph.source(_in_edge); u != join; u = _parent[u]) {
          _flow[_pred_edge[u]] += _forward[u] ? -val : val;
        }
        for (Node u = _graph.target(_in_edge); u != join; u = _parent[u]) {
          _flow[_pred_edge[u]] += _forward[u] ? val : -val;
        }
      }
      // Updating the state of the entering and leaving edges
      if (change) {
        _state[_in_edge] = STATE_TREE;
        _state[_pred_edge[u_out]] =
          (_flow[_pred_edge[u_out]] == 0) ? STATE_LOWER : STATE_UPPER;
      } else {
        _state[_in_edge] = -_state[_in_edge];
      }
    }

    /// Updates \c thread and \c parent node maps.
    void updateThreadParent() {
      Node u;
      v_out = _parent[u_out];

      // Handling the case when join and v_out coincide
      bool par_first = false;
      if (join == v_out) {
        for (u = join; u != u_in && u != v_in; u = _thread[u]) ;
        if (u == v_in) {
          par_first = true;
          while (_thread[u] != u_out) u = _thread[u];
          first = u;
        }
      }

      // Finding the last successor of u_in (u) and the node after it
      // (right) according to the thread index
      for (u = u_in; _depth[_thread[u]] > _depth[u_in]; u = _thread[u]) ;
      right = _thread[u];
      if (_thread[v_in] == u_out) {
        for (last = u; _depth[last] > _depth[u_out]; last = _thread[last]) ;
        if (last == u_out) last = _thread[last];
      }
      else last = _thread[v_in];

      // Updating stem nodes
      _thread[v_in] = stem = u_in;
      par_stem = v_in;
      while (stem != u_out) {
        _thread[u] = new_stem = _parent[stem];

        // Finding the node just before the stem node (u) according to
        // the original thread index
        for (u = new_stem; _thread[u] != stem; u = _thread[u]) ;
        _thread[u] = right;

        // Changing the parent node of stem and shifting stem and
        // par_stem nodes
        _parent[stem] = par_stem;
        par_stem = stem;
        stem = new_stem;

        // Finding the last successor of stem (u) and the node after it
        // (right) according to the thread index
        for (u = stem; _depth[_thread[u]] > _depth[stem]; u = _thread[u]) ;
        right = _thread[u];
      }
      _parent[u_out] = par_stem;
      _thread[u] = last;

      if (join == v_out && par_first) {
        if (first != v_in) _thread[first] = right;
      } else {
        for (u = v_out; _thread[u] != u_out; u = _thread[u]) ;
        _thread[u] = right;
      }
    }

    /// Updates \c pred_edge and \c forward node maps.
    void updatePredEdge() {
      Node u = u_out, v;
      while (u != u_in) {
        v = _parent[u];
        _pred_edge[u] = _pred_edge[v];
        _forward[u] = !_forward[v];
        u = v;
      }
      _pred_edge[u_in] = _in_edge;
      _forward[u_in] = (u_in == _graph.source(_in_edge));
    }

    /// Updates \c depth and \c potential node maps.
    void updateDepthPotential() {
      _depth[u_in] = _depth[v_in] + 1;
      _potential[u_in] = _forward[u_in] ?
        _potential[v_in] - _cost[_pred_edge[u_in]] :
        _potential[v_in] + _cost[_pred_edge[u_in]];

      Node u = _thread[u_in], v;
      while (true) {
        v = _parent[u];
        if (v == INVALID) break;
        _depth[u] = _depth[v] + 1;
        _potential[u] = _forward[u] ?
          _potential[v] - _cost[_pred_edge[u]] :
          _potential[v] + _cost[_pred_edge[u]];
        if (_depth[u] <= _depth[v_in]) break;
        u = _thread[u];
      }
    }

    /// Executes the algorithm.
    bool start(PivotRuleEnum pivot_rule) {
      switch (pivot_rule) {
        case FIRST_ELIGIBLE_PIVOT:
          return start<FirstEligiblePivotRule>();
        case BEST_ELIGIBLE_PIVOT:
          return start<BestEligiblePivotRule>();
        case BLOCK_SEARCH_PIVOT:
          return start<BlockSearchPivotRule>();
        case LIMITED_SEARCH_PIVOT:
          return start<LimitedSearchPivotRule>();
        case CANDIDATE_LIST_PIVOT:
          return start<CandidateListPivotRule>();
        case COMBINED_PIVOT:
          if ( countEdges(_graph) / countNodes(_graph) <=
               COMBINED_PIVOT_MAX_DEG )
            return start<LimitedSearchPivotRule>();
          else
            return start<BlockSearchPivotRule>();
      }
      return false;
    }

    template<class PivotRuleImplementation>
    bool start() {
      PivotRuleImplementation pivot(*this);

      // Executing the network simplex algorithm
      while (pivot.findEnteringEdge()) {
        join = findJoinNode();
        bool change = findLeavingEdge();
        changeFlows(change);
        if (change) {
          updateThreadParent();
          updatePredEdge();
          updateDepthPotential();
        }
      }

      // Checking if the flow amount equals zero on all the artificial
      // edges
      for (InEdgeIt e(_graph, _root); e != INVALID; ++e)
        if (_flow[e] > 0) return false;
      for (OutEdgeIt e(_graph, _root); e != INVALID; ++e)
        if (_flow[e] > 0) return false;

      // Copying flow values to _flow_result
      if (_lower) {
        for (typename Graph::EdgeIt e(_graph_ref); e != INVALID; ++e)
          (*_flow_result)[e] = (*_lower)[e] + _flow[_edge_ref[e]];
      } else {
        for (typename Graph::EdgeIt e(_graph_ref); e != INVALID; ++e)
          (*_flow_result)[e] = _flow[_edge_ref[e]];
      }
      // Copying potential values to _potential_result
      for (typename Graph::NodeIt n(_graph_ref); n != INVALID; ++n)
        (*_potential_result)[n] = _potential[_node_ref[n]];

      return true;
    }

  }; //class NetworkSimplex

  ///@}

} //namespace lemon

#endif //LEMON_NETWORK_SIMPLEX_H
