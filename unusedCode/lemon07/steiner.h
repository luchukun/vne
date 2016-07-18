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

#ifndef LEMON_STEINER_H
#define LEMON_STEINER_H

///\ingroup approx
///\file
///\brief Algorithm for the 2-approximation of Steiner Tree problem.
///

#include <lemon/smart_graph.h>
#include <lemon/graph_utils.h>
#include <lemon/error.h>

#include <lemon/ugraph_adaptor.h>
#include <lemon/maps.h>

#include <lemon/dijkstra.h>
#include <lemon/prim.h>


namespace lemon {

  /// \ingroup approx
  
  /// \brief Algorithm for the 2-approximation of Steiner Tree problem
  ///
  /// The Steiner-tree problem is the next: Given a connected
  /// undirected graph, a cost function on the edges and a subset of
  /// the nodes. Construct a tree with minimum cost which covers the
  /// given subset of the nodes. The problem is NP-hard moreover
  /// it is APX-complete too.
  ///
  /// Mehlhorn's approximation algorithm is implemented in this class,
  /// which gives a 2-approximation for the Steiner-tree problem. The
  /// algorithm's time complexity is O(nlog(n)+e).
  template <typename UGraph,
            typename CostMap = typename UGraph:: template UEdgeMap<double> >
  class SteinerTree {
  public:
    
    UGRAPH_TYPEDEFS(typename UGraph);

    typedef typename CostMap::Value Value;
    
  private:

    class CompMap {
    public:
      typedef Node Key;
      typedef Edge Value;

    private:
      const UGraph& _graph;
      typename UGraph::template NodeMap<int> _comp;

    public:
      CompMap(const UGraph& graph) : _graph(graph), _comp(graph) {}

      void set(const Node& node, const Edge& edge) {
        if (edge != INVALID) {
          _comp.set(node, _comp[_graph.source(edge)]);
        } else {
          _comp.set(node, -1);
        }
      }

      int comp(const Node& node) const { return _comp[node]; }
      void comp(const Node& node, int value) { _comp.set(node, value); }
    };

    typedef typename UGraph::template NodeMap<Edge> PredMap;

    typedef ForkWriteMap<PredMap, CompMap> ForkedMap;


    struct External {
      int source, target;
      UEdge uedge;
      Value value;

      External(int s, int t, const UEdge& e, const Value& v)
        : source(s), target(t), uedge(e), value(v) {}
    };

    struct ExternalLess {
      bool operator()(const External& left, const External& right) const {
        return (left.source < right.source) || 
          (left.source == right.source && left.target < right.target);
      }
    };


    typedef typename UGraph::template NodeMap<bool> FilterMap;

    typedef typename UGraph::template UEdgeMap<bool> TreeMap;

    const UGraph& _graph;
    const CostMap& _cost;

    typename Dijkstra<UGraph, CostMap>::
    template DefPredMap<ForkedMap>::Create _dijkstra;

    PredMap* _pred;
    CompMap* _comp;
    ForkedMap* _forked;

    int _terminal_num;

    FilterMap *_filter;
    TreeMap *_tree;

    Value _value;

  public:

    /// \brief Constructor
    
    /// Constructor
    ///
    SteinerTree(const UGraph &graph, const CostMap &cost)
      : _graph(graph), _cost(cost), _dijkstra(graph, _cost), 
        _pred(0), _comp(0), _forked(0), _filter(0), _tree(0) {}

    /// \brief Initializes the internal data structures.
    ///
    /// Initializes the internal data structures.
    void init() {
      if (!_pred) _pred = new PredMap(_graph);
      if (!_comp) _comp = new CompMap(_graph);
      if (!_forked) _forked = new ForkedMap(*_pred, *_comp);
      if (!_filter) _filter = new FilterMap(_graph);
      if (!_tree) _tree = new TreeMap(_graph);
      _dijkstra.predMap(*_forked);
      _dijkstra.init();
      _terminal_num = 0;
      for (NodeIt it(_graph); it != INVALID; ++it) {
        _filter->set(it, false);
      }
    }

    /// \brief Adds a new terminal node.
    ///
    /// Adds a new terminal node to the Steiner-tree problem.
    void addTerminal(const Node& node) {
      if (!_dijkstra.reached(node)) {
        _dijkstra.addSource(node);
        _comp->comp(node, _terminal_num);
        ++_terminal_num;
      }
    }
    
    /// \brief Executes the algorithm.
    ///
    /// Executes the algorithm.
    ///
    /// \pre init() must be called and at least some nodes should be
    /// added with addTerminal() before using this function.
    ///
    /// This method constructs an approximation of the Steiner-Tree.
    void start() {
      _dijkstra.start();
      
      std::vector<External> externals;
      for (UEdgeIt it(_graph); it != INVALID; ++it) {
        Node s = _graph.source(it);
        Node t = _graph.target(it);
        if (_comp->comp(s) == _comp->comp(t)) continue;

        Value cost = _dijkstra.dist(s) + _dijkstra.dist(t) + _cost[it];

        if (_comp->comp(s) < _comp->comp(t)) {
          externals.push_back(External(_comp->comp(s), _comp->comp(t), 
                                       it, cost));
        } else {
          externals.push_back(External(_comp->comp(t), _comp->comp(s), 
                                       it, cost));
        }
      }
      std::sort(externals.begin(), externals.end(), ExternalLess());

      SmartUGraph aux_graph;
      std::vector<SmartUGraph::Node> aux_nodes;

      for (int i = 0; i < _terminal_num; ++i) {
        aux_nodes.push_back(aux_graph.addNode());
      }

      SmartUGraph::UEdgeMap<Value> aux_cost(aux_graph);
      SmartUGraph::UEdgeMap<UEdge> cross(aux_graph);
      {
        int i = 0;
        while (i < int(externals.size())) {
          int sn = externals[i].source;
          int tn = externals[i].target;
          Value ev = externals[i].value;
          UEdge ee = externals[i].uedge;
          ++i;
          while (i < int(externals.size()) && 
                 sn == externals[i].source && tn == externals[i].target) {
            if (externals[i].value < ev) {
              ev = externals[i].value;
              ee = externals[i].uedge;
            }
            ++i;
          }
          SmartUGraph::UEdge ne = 
            aux_graph.addEdge(aux_nodes[sn], aux_nodes[tn]);
          aux_cost.set(ne, ev);
          cross.set(ne, ee);
        }
      }

      std::vector<SmartUGraph::UEdge> aux_tree_edges;
      BackInserterBoolMap<std::vector<SmartUGraph::UEdge> >
        aux_tree_map(aux_tree_edges);
      prim(aux_graph, aux_cost, aux_tree_map);

      for (std::vector<SmartUGraph::UEdge>::iterator 
             it = aux_tree_edges.begin(); it != aux_tree_edges.end(); ++it) {
        Node node;
        node = _graph.source(cross[*it]);
        while (node != INVALID && !(*_filter)[node]) {
          _filter->set(node, true);
          node = (*_pred)[node] != INVALID ? 
            _graph.source((*_pred)[node]) : INVALID;
        }
        node = _graph.target(cross[*it]);
        while (node != INVALID && !(*_filter)[node]) {
          _filter->set(node, true);
          node = (*_pred)[node] != INVALID ? 
            _graph.source((*_pred)[node]) : INVALID;
        }
      }

      _value = prim(nodeSubUGraphAdaptor(_graph, *_filter), _cost, *_tree);
            
    }

    /// \brief Checks if an edge is in the Steiner-tree or not.
    ///
    /// Checks if an edge is in the Steiner-tree or not.
    /// \param e is the edge that will be checked
    /// \return \c true if e is in the Steiner-tree, \c false otherwise
    bool tree(UEdge e){
      return (*_tree)[e];
    }

    /// \brief Checks if the node is in the Steiner-tree or not.
    ///
    /// Checks if a node is in the Steiner-tree or not.
    /// \param n is the node that will be checked
    /// \return \c true if n is in the Steiner-tree, \c false otherwise
    bool tree(Node n){
      return (*_filter)[n];
    }

    /// \brief Checks if the node is a Steiner-node.
    ///
    /// Checks if the node is a Steiner-node (i.e. a tree node but
    /// not terminal).
    /// \param n is the node that will be checked
    /// \return \c true if n is a Steiner-node, \c false otherwise
    bool steiner(Node n){
      return (*_filter)[n] && (*_pred)[n] != INVALID;
    }

    /// \brief Checks if the node is a terminal.
    ///
    /// Checks if the node is a terminal.
    /// \param n is the node that will be checked
    /// \return \c true if n is a terminal, \c false otherwise
    bool terminal(Node n){
      return _dijkstra.reached(n) && (*_pred)[n] == INVALID;
    }
    
    /// \brief The total cost of the tree
    ///
    /// The total cost of the constructed tree. The calculated value does
    /// not exceed the double of the optimal value.
    Value treeValue() const {
      return _value;
    }
    
  };

} //END OF NAMESPACE LEMON

#endif
