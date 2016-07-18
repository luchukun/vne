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

#ifndef LEMON_PLANARITY_H
#define LEMON_PLANARITY_H

/// \ingroup planar
/// \file
/// \brief Planarity checking, embedding, drawing and coloring

#include <vector>
#include <list>

#include <lemon/dfs.h>
#include <lemon/bfs.h>
#include <lemon/radix_sort.h>
#include <lemon/maps.h>
#include <lemon/path.h>
#include <lemon/iterable_maps.h>
#include <lemon/edge_set.h>
#include <lemon/bucket_heap.h>
#include <lemon/ugraph_adaptor.h>
#include <lemon/color.h>


namespace lemon {

  namespace _planarity_bits {

    template <typename UGraph>
    struct PlanarityVisitor : DfsVisitor<UGraph> {

      typedef typename UGraph::Node Node;
      typedef typename UGraph::Edge Edge;

      typedef typename UGraph::template NodeMap<Edge> PredMap;

      typedef typename UGraph::template UEdgeMap<bool> TreeMap;

      typedef typename UGraph::template NodeMap<int> OrderMap;
      typedef std::vector<Node> OrderList;

      typedef typename UGraph::template NodeMap<int> LowMap;
      typedef typename UGraph::template NodeMap<int> AncestorMap;

      PlanarityVisitor(const UGraph& ugraph,
		       PredMap& pred_map, TreeMap& tree_map,
		       OrderMap& order_map, OrderList& order_list,
		       AncestorMap& ancestor_map, LowMap& low_map)
	: _ugraph(ugraph), _pred_map(pred_map), _tree_map(tree_map),
	  _order_map(order_map), _order_list(order_list),
	  _ancestor_map(ancestor_map), _low_map(low_map) {}
      
      void reach(const Node& node) {
	_order_map[node] = _order_list.size();
	_low_map[node] = _order_list.size();
	_ancestor_map[node] = _order_list.size();
	_order_list.push_back(node);
      }

      void discover(const Edge& edge) {
	Node source = _ugraph.source(edge);
	Node target = _ugraph.target(edge);

	_tree_map[edge] = true;
	_pred_map[target] = edge;
      }

      void examine(const Edge& edge) {
	Node source = _ugraph.source(edge);
	Node target = _ugraph.target(edge);
	
	if (_order_map[target] < _order_map[source] && !_tree_map[edge]) {
	  if (_low_map[source] > _order_map[target]) {
	    _low_map[source] = _order_map[target];
	  }
	  if (_ancestor_map[source] > _order_map[target]) {
	    _ancestor_map[source] = _order_map[target];
	  }
	}
      }

      void backtrack(const Edge& edge) {
	Node source = _ugraph.source(edge);
	Node target = _ugraph.target(edge);

	if (_low_map[source] > _low_map[target]) {
	  _low_map[source] = _low_map[target];
	}
      }

      const UGraph& _ugraph;
      PredMap& _pred_map;
      TreeMap& _tree_map;
      OrderMap& _order_map;
      OrderList& _order_list;
      AncestorMap& _ancestor_map;
      LowMap& _low_map;
    };

    template <typename UGraph, bool embedding = true>
    struct NodeDataNode {
      int prev, next;
      int visited;
      typename UGraph::Edge first;
      bool inverted;
    };

    template <typename UGraph>
    struct NodeDataNode<UGraph, false> {
      int prev, next;
      int visited;
    };

    template <typename UGraph>
    struct ChildListNode {
      typedef typename UGraph::Node Node;
      Node first;
      Node prev, next;
    };

    template <typename UGraph>
    struct EdgeListNode {
      typename UGraph::Edge prev, next;
    };

  }

  /// \ingroup planar
  ///
  /// \brief Planarity checking of an undirected simple graph
  ///
  /// This class implements the Boyer-Myrvold algorithm for planarity
  /// checking of an undirected graph. This class is a simplified
  /// version of the PlanarEmbedding algorithm class, and it does
  /// provide neither embedding nor kuratowski subdivisons.
  template <typename UGraph>
  class PlanarityChecking {
  private:
    
    UGRAPH_TYPEDEFS(typename UGraph);
      
    const UGraph& _ugraph;

  private:

    typedef typename UGraph::template NodeMap<Edge> PredMap;
    
    typedef typename UGraph::template UEdgeMap<bool> TreeMap;

    typedef typename UGraph::template NodeMap<int> OrderMap;
    typedef std::vector<Node> OrderList;

    typedef typename UGraph::template NodeMap<int> LowMap;
    typedef typename UGraph::template NodeMap<int> AncestorMap;

    typedef _planarity_bits::NodeDataNode<UGraph> NodeDataNode;
    typedef std::vector<NodeDataNode> NodeData;
    
    typedef _planarity_bits::ChildListNode<UGraph> ChildListNode;
    typedef typename UGraph::template NodeMap<ChildListNode> ChildLists;

    typedef typename UGraph::template NodeMap<std::list<int> > MergeRoots;
 
    typedef typename UGraph::template NodeMap<bool> EmbedEdge;

  public:

    /// \brief Constructor
    ///
    /// \warining The graph should be simple, i.e. parallel and loop edge
    /// free.
    PlanarityChecking(const UGraph& ugraph) : _ugraph(ugraph) {}

    /// \brief Runs the algorithm.
    ///
    /// Runs the algorithm.  
    /// \return %True when the graph is planar.
    bool run() {
      typedef _planarity_bits::PlanarityVisitor<UGraph> Visitor;

      PredMap pred_map(_ugraph, INVALID);
      TreeMap tree_map(_ugraph, false);

      OrderMap order_map(_ugraph, -1);
      OrderList order_list;

      AncestorMap ancestor_map(_ugraph, -1);
      LowMap low_map(_ugraph, -1);

      Visitor visitor(_ugraph, pred_map, tree_map,
		      order_map, order_list, ancestor_map, low_map);
      DfsVisit<UGraph, Visitor> visit(_ugraph, visitor);
      visit.run();

      ChildLists child_lists(_ugraph);
      createChildLists(tree_map, order_map, low_map, child_lists);

      NodeData node_data(2 * order_list.size());
      
      EmbedEdge embed_edge(_ugraph, false);

      MergeRoots merge_roots(_ugraph);
      
      for (int i = order_list.size() - 1; i >= 0; --i) {

	Node node = order_list[i];

	Node source = node;
	for (OutEdgeIt e(_ugraph, node); e != INVALID; ++e) {
	  Node target = _ugraph.target(e);
	  
	  if (order_map[source] < order_map[target] && tree_map[e]) {
	    initFace(target, node_data, order_map, order_list);
	  }
	}
	
	for (OutEdgeIt e(_ugraph, node); e != INVALID; ++e) {
	  Node target = _ugraph.target(e);
	  
	  if (order_map[source] < order_map[target] && !tree_map[e]) {
	    embed_edge[target] = true;
	    walkUp(target, source, i, pred_map, low_map,
		   order_map, order_list, node_data, merge_roots);
	  }
	}

	for (typename MergeRoots::Value::iterator it = 
	       merge_roots[node].begin(); it != merge_roots[node].end(); ++it) {
	  int rn = *it;
	  walkDown(rn, i, node_data, order_list, child_lists, 
		   ancestor_map, low_map, embed_edge, merge_roots);
	}
	merge_roots[node].clear();
	
	for (OutEdgeIt e(_ugraph, node); e != INVALID; ++e) {
	  Node target = _ugraph.target(e);
	  
	  if (order_map[source] < order_map[target] && !tree_map[e]) {
	    if (embed_edge[target]) {
	      return false;
	    }
	  }
	}
      }

      return true;
    }
    
  private:

    void createChildLists(const TreeMap& tree_map, const OrderMap& order_map,
			  const LowMap& low_map, ChildLists& child_lists) {

      for (NodeIt n(_ugraph); n != INVALID; ++n) {
	Node source = n;
	
	std::vector<Node> targets;  
	for (OutEdgeIt e(_ugraph, n); e != INVALID; ++e) {
	  Node target = _ugraph.target(e);

	  if (order_map[source] < order_map[target] && tree_map[e]) {
	    targets.push_back(target);
	  }
	}	

	if (targets.size() == 0) {
	  child_lists[source].first = INVALID;
	} else if (targets.size() == 1) {
	  child_lists[source].first = targets[0];
	  child_lists[targets[0]].prev = INVALID;
	  child_lists[targets[0]].next = INVALID;
	} else {
	  radixSort(targets.begin(), targets.end(), mapFunctor(low_map));
	  for (int i = 1; i < int(targets.size()); ++i) {
	    child_lists[targets[i]].prev = targets[i - 1];
	    child_lists[targets[i - 1]].next = targets[i];
	  }
	  child_lists[targets.back()].next = INVALID; 
	  child_lists[targets.front()].prev = INVALID;
	  child_lists[source].first = targets.front();
	}
      }
    }

    void walkUp(const Node& node, Node root, int rorder,
		const PredMap& pred_map, const LowMap& low_map,
		const OrderMap& order_map, const OrderList& order_list,
		NodeData& node_data, MergeRoots& merge_roots) {

      int na, nb;
      bool da, db;
      
      na = nb = order_map[node];
      da = true; db = false;
      
      while (true) {
	
	if (node_data[na].visited == rorder) break;
	if (node_data[nb].visited == rorder) break;

	node_data[na].visited = rorder;
	node_data[nb].visited = rorder;

	int rn = -1;

	if (na >= int(order_list.size())) {
	  rn = na;
	} else if (nb >= int(order_list.size())) {
	  rn = nb;
	}

	if (rn == -1) {
	  int nn;
	  
	  nn = da ? node_data[na].prev : node_data[na].next;
	  da = node_data[nn].prev != na;
	  na = nn;
	  
	  nn = db ? node_data[nb].prev : node_data[nb].next;
	  db = node_data[nn].prev != nb;
	  nb = nn;

	} else {

	  Node rep = order_list[rn - order_list.size()];
	  Node parent = _ugraph.source(pred_map[rep]);

	  if (low_map[rep] < rorder) {
	    merge_roots[parent].push_back(rn);
	  } else {
	    merge_roots[parent].push_front(rn);
	  }
	  
	  if (parent != root) {  
	    na = nb = order_map[parent];
	    da = true; db = false;
	  } else {
	    break;
	  }
	}	
      }
    }

    void walkDown(int rn, int rorder, NodeData& node_data,
		  OrderList& order_list, ChildLists& child_lists,
		  AncestorMap& ancestor_map, LowMap& low_map,
		  EmbedEdge& embed_edge, MergeRoots& merge_roots) {

      std::vector<std::pair<int, bool> > merge_stack;

      for (int di = 0; di < 2; ++di) {
	bool rd = di == 0;
	int pn = rn;
	int n = rd ? node_data[rn].next : node_data[rn].prev;
	
	while (n != rn) {
	  
	  Node node = order_list[n];
	  
	  if (embed_edge[node]) {

	    // Merging components on the critical path
	    while (!merge_stack.empty()) {

	      // Component root
	      int cn = merge_stack.back().first;
	      bool cd = merge_stack.back().second;
	      merge_stack.pop_back();

	      // Parent of component
	      int dn = merge_stack.back().first;
	      bool dd = merge_stack.back().second;
	      merge_stack.pop_back();

	      Node parent = order_list[dn];

	      // Erasing from merge_roots
	      merge_roots[parent].pop_front();
	      
	      Node child = order_list[cn - order_list.size()];

	      // Erasing from child_lists
	      if (child_lists[child].prev != INVALID) {
		child_lists[child_lists[child].prev].next =
		  child_lists[child].next;
	      } else {
		child_lists[parent].first = child_lists[child].next;
	      }
	      
	      if (child_lists[child].next != INVALID) {
		child_lists[child_lists[child].next].prev =
		  child_lists[child].prev;
	      }
	      
	      // Merging external faces
	      {
		int en = cn;
		cn = cd ? node_data[cn].prev : node_data[cn].next;
		cd = node_data[cn].next == en;

	      }

	      if (cd) node_data[cn].next = dn; else node_data[cn].prev = dn;
	      if (dd) node_data[dn].prev = cn; else node_data[dn].next = cn;

	    }

	    bool d = pn == node_data[n].prev;

	    if (node_data[n].prev == node_data[n].next && 
		node_data[n].inverted) {
	      d = !d;
	    }

	    // Embedding edge into external face
	    if (rd) node_data[rn].next = n; else node_data[rn].prev = n;
	    if (d) node_data[n].prev = rn; else node_data[n].next = rn;
	    pn = rn;

	    embed_edge[order_list[n]] = false;
	  }

	  if (!merge_roots[node].empty()) {

	    bool d = pn == node_data[n].prev;

	    merge_stack.push_back(std::make_pair(n, d));

	    int rn = merge_roots[node].front();
	    
	    int xn = node_data[rn].next;
	    Node xnode = order_list[xn];
	    
	    int yn = node_data[rn].prev;
	    Node ynode = order_list[yn];
	    
	    bool rd;
	    if (!external(xnode, rorder, child_lists, ancestor_map, low_map)) {
	      rd = true;
	    } else if (!external(ynode, rorder, child_lists, 
				 ancestor_map, low_map)) {
	      rd = false;
	    } else if (pertinent(xnode, embed_edge, merge_roots)) {
	      rd = true;
	    } else {
	      rd = false;
	    }
	    
	    merge_stack.push_back(std::make_pair(rn, rd));
	    
	    pn = rn;
	    n = rd ? xn : yn;	      
	    	    
	  } else if (!external(node, rorder, child_lists,
			       ancestor_map, low_map)) {
	    int nn = (node_data[n].next != pn ? 
		      node_data[n].next : node_data[n].prev);

	    bool nd = n == node_data[nn].prev;

	    if (nd) node_data[nn].prev = pn;
	    else node_data[nn].next = pn; 

	    if (n == node_data[pn].prev) node_data[pn].prev = nn;
	    else node_data[pn].next = nn;

	    node_data[nn].inverted = 
	      (node_data[nn].prev == node_data[nn].next && nd != rd);

	    n = nn;
	  }
	  else break;
	  
	}

	if (!merge_stack.empty() || n == rn) {
	  break;
	}
      }
    }

    void initFace(const Node& node, NodeData& node_data, 
		  const OrderMap& order_map, const OrderList& order_list) {
      int n = order_map[node];
      int rn = n + order_list.size();
      
      node_data[n].next = node_data[n].prev = rn;
      node_data[rn].next = node_data[rn].prev = n;
      
      node_data[n].visited = order_list.size();
      node_data[rn].visited = order_list.size();
      
    }

    bool external(const Node& node, int rorder,
		  ChildLists& child_lists, AncestorMap& ancestor_map, 
		  LowMap& low_map) {
      Node child = child_lists[node].first;

      if (child != INVALID) {
	if (low_map[child] < rorder) return true;
      }

      if (ancestor_map[node] < rorder) return true;

      return false;
    }

    bool pertinent(const Node& node, const EmbedEdge& embed_edge,
		   const MergeRoots& merge_roots) {
      return !merge_roots[node].empty() || embed_edge[node];
    }

  };

  /// \ingroup planar
  ///
  /// \brief Planar embedding of an undirected simple graph
  ///
  /// This class implements the Boyer-Myrvold algorithm for planar
  /// embedding of an undirected graph. The planar embeding is an
  /// ordering of the outgoing edges in each node, which is a possible
  /// configuration to draw the graph in the plane. If there is not
  /// such ordering then the graph contains a \f$ K_5 \f$ (full graph
  /// with 5 nodes) or an \f$ K_{3,3} \f$ (complete bipartite graph on
  /// 3 ANode and 3 BNode) subdivision.
  ///
  /// The current implementation calculates an embedding or an
  /// Kuratowski subdivision if the graph is not planar. The running
  /// time of the algorithm is \f$ O(n) \f$.
  template <typename UGraph>
  class PlanarEmbedding {
  private:
    
    UGRAPH_TYPEDEFS(typename UGraph);
      
    const UGraph& _ugraph;
    typename UGraph::template EdgeMap<Edge> _embedding;

    typename UGraph::template UEdgeMap<bool> _kuratowski;

  private:

    typedef typename UGraph::template NodeMap<Edge> PredMap;
    
    typedef typename UGraph::template UEdgeMap<bool> TreeMap;

    typedef typename UGraph::template NodeMap<int> OrderMap;
    typedef std::vector<Node> OrderList;

    typedef typename UGraph::template NodeMap<int> LowMap;
    typedef typename UGraph::template NodeMap<int> AncestorMap;

    typedef _planarity_bits::NodeDataNode<UGraph> NodeDataNode;
    typedef std::vector<NodeDataNode> NodeData;
    
    typedef _planarity_bits::ChildListNode<UGraph> ChildListNode;
    typedef typename UGraph::template NodeMap<ChildListNode> ChildLists;

    typedef typename UGraph::template NodeMap<std::list<int> > MergeRoots;
 
    typedef typename UGraph::template NodeMap<Edge> EmbedEdge;

    typedef _planarity_bits::EdgeListNode<UGraph> EdgeListNode;
    typedef typename UGraph::template EdgeMap<EdgeListNode> EdgeLists;

    typedef typename UGraph::template NodeMap<bool> FlipMap;

    typedef typename UGraph::template NodeMap<int> TypeMap;

    enum IsolatorNodeType {
      HIGHX = 6, LOWX = 7,
      HIGHY = 8, LOWY = 9,
      ROOT = 10, PERTINENT = 11,
      INTERNAL = 12
    };

  public:

    /// \brief The map for store of embedding
    typedef typename UGraph::template EdgeMap<Edge> EmbeddingMap;

    /// \brief Constructor
    ///
    /// \warining The graph should be simple, i.e. parallel and loop edge
    /// free.
    PlanarEmbedding(const UGraph& ugraph)
      : _ugraph(ugraph), _embedding(_ugraph), _kuratowski(ugraph, false) {}

    /// \brief Runs the algorithm.
    ///
    /// Runs the algorithm.  
    /// \param kuratowski If the parameter is false, then the
    /// algorithm does not calculate the isolate Kuratowski
    /// subdivisions.
    ///\return %True when the graph is planar.
    bool run(bool kuratowski = true) {
      typedef _planarity_bits::PlanarityVisitor<UGraph> Visitor;

      PredMap pred_map(_ugraph, INVALID);
      TreeMap tree_map(_ugraph, false);

      OrderMap order_map(_ugraph, -1);
      OrderList order_list;

      AncestorMap ancestor_map(_ugraph, -1);
      LowMap low_map(_ugraph, -1);

      Visitor visitor(_ugraph, pred_map, tree_map,
		      order_map, order_list, ancestor_map, low_map);
      DfsVisit<UGraph, Visitor> visit(_ugraph, visitor);
      visit.run();

      ChildLists child_lists(_ugraph);
      createChildLists(tree_map, order_map, low_map, child_lists);

      NodeData node_data(2 * order_list.size());
      
      EmbedEdge embed_edge(_ugraph, INVALID);

      MergeRoots merge_roots(_ugraph);

      EdgeLists edge_lists(_ugraph);

      FlipMap flip_map(_ugraph, false);

      for (int i = order_list.size() - 1; i >= 0; --i) {

	Node node = order_list[i];

	node_data[i].first = INVALID;
	
	Node source = node;
	for (OutEdgeIt e(_ugraph, node); e != INVALID; ++e) {
	  Node target = _ugraph.target(e);
	  
	  if (order_map[source] < order_map[target] && tree_map[e]) {
	    initFace(target, edge_lists, node_data,
		      pred_map, order_map, order_list);
	  }
	}
	
	for (OutEdgeIt e(_ugraph, node); e != INVALID; ++e) {
	  Node target = _ugraph.target(e);
	  
	  if (order_map[source] < order_map[target] && !tree_map[e]) {
	    embed_edge[target] = e;
	    walkUp(target, source, i, pred_map, low_map,
		   order_map, order_list, node_data, merge_roots);
	  }
	}

	for (typename MergeRoots::Value::iterator it = 
	       merge_roots[node].begin(); it != merge_roots[node].end(); ++it) {
	  int rn = *it;
	  walkDown(rn, i, node_data, edge_lists, flip_map, order_list, 
		   child_lists, ancestor_map, low_map, embed_edge, merge_roots);
	}
	merge_roots[node].clear();
	
	for (OutEdgeIt e(_ugraph, node); e != INVALID; ++e) {
	  Node target = _ugraph.target(e);
	  
	  if (order_map[source] < order_map[target] && !tree_map[e]) {
	    if (embed_edge[target] != INVALID) {
	      if (kuratowski) {
		isolateKuratowski(e, node_data, edge_lists, flip_map,
				  order_map, order_list, pred_map, child_lists,
				  ancestor_map, low_map, 
				  embed_edge, merge_roots);
	      }
	      return false;
	    }
	  }
	}
      }

      for (int i = 0; i < int(order_list.size()); ++i) {

	mergeRemainingFaces(order_list[i], node_data, order_list, order_map,
			    child_lists, edge_lists);
	storeEmbedding(order_list[i], node_data, order_map, pred_map,
		       edge_lists, flip_map);
      }

      return true;
    }

    /// \brief Gives back the successor of an edge
    ///
    /// Gives back the successor of an edge. This function makes
    /// possible to query the cyclic order of the outgoing edges from
    /// a node.
    Edge next(const Edge& edge) const {
      return _embedding[edge];
    }

    /// \brief Gives back the calculated embedding map
    ///
    /// The returned map contains the successor of each edge in the
    /// graph.
    const EmbeddingMap& embedding() const {
      return _embedding;
    }

    /// \brief Gives back true when the undirected edge is in the
    /// kuratowski subdivision
    ///
    /// Gives back true when the undirected edge is in the kuratowski
    /// subdivision
    bool kuratowski(const UEdge& uedge) {
      return _kuratowski[uedge];
    }

  private:

    void createChildLists(const TreeMap& tree_map, const OrderMap& order_map,
			  const LowMap& low_map, ChildLists& child_lists) {

      for (NodeIt n(_ugraph); n != INVALID; ++n) {
	Node source = n;
	
	std::vector<Node> targets;  
	for (OutEdgeIt e(_ugraph, n); e != INVALID; ++e) {
	  Node target = _ugraph.target(e);

	  if (order_map[source] < order_map[target] && tree_map[e]) {
	    targets.push_back(target);
	  }
	}	

	if (targets.size() == 0) {
	  child_lists[source].first = INVALID;
	} else if (targets.size() == 1) {
	  child_lists[source].first = targets[0];
	  child_lists[targets[0]].prev = INVALID;
	  child_lists[targets[0]].next = INVALID;
	} else {
	  radixSort(targets.begin(), targets.end(), mapFunctor(low_map));
	  for (int i = 1; i < int(targets.size()); ++i) {
	    child_lists[targets[i]].prev = targets[i - 1];
	    child_lists[targets[i - 1]].next = targets[i];
	  }
	  child_lists[targets.back()].next = INVALID; 
	  child_lists[targets.front()].prev = INVALID;
	  child_lists[source].first = targets.front();
	}
      }
    }

    void walkUp(const Node& node, Node root, int rorder,
		const PredMap& pred_map, const LowMap& low_map,
		const OrderMap& order_map, const OrderList& order_list,
		NodeData& node_data, MergeRoots& merge_roots) {

      int na, nb;
      bool da, db;
      
      na = nb = order_map[node];
      da = true; db = false;
      
      while (true) {
	
	if (node_data[na].visited == rorder) break;
	if (node_data[nb].visited == rorder) break;

	node_data[na].visited = rorder;
	node_data[nb].visited = rorder;

	int rn = -1;

	if (na >= int(order_list.size())) {
	  rn = na;
	} else if (nb >= int(order_list.size())) {
	  rn = nb;
	}

	if (rn == -1) {
	  int nn;
	  
	  nn = da ? node_data[na].prev : node_data[na].next;
	  da = node_data[nn].prev != na;
	  na = nn;
	  
	  nn = db ? node_data[nb].prev : node_data[nb].next;
	  db = node_data[nn].prev != nb;
	  nb = nn;

	} else {

	  Node rep = order_list[rn - order_list.size()];
	  Node parent = _ugraph.source(pred_map[rep]);

	  if (low_map[rep] < rorder) {
	    merge_roots[parent].push_back(rn);
	  } else {
	    merge_roots[parent].push_front(rn);
	  }
	  
	  if (parent != root) {  
	    na = nb = order_map[parent];
	    da = true; db = false;
	  } else {
	    break;
	  }
	}	
      }
    }

    void walkDown(int rn, int rorder, NodeData& node_data,
		  EdgeLists& edge_lists, FlipMap& flip_map, 
		  OrderList& order_list, ChildLists& child_lists,
		  AncestorMap& ancestor_map, LowMap& low_map,
		  EmbedEdge& embed_edge, MergeRoots& merge_roots) {

      std::vector<std::pair<int, bool> > merge_stack;

      for (int di = 0; di < 2; ++di) {
	bool rd = di == 0;
	int pn = rn;
	int n = rd ? node_data[rn].next : node_data[rn].prev;
	
	while (n != rn) {
	  
	  Node node = order_list[n];
	  
	  if (embed_edge[node] != INVALID) {

	    // Merging components on the critical path
	    while (!merge_stack.empty()) {

	      // Component root
	      int cn = merge_stack.back().first;
	      bool cd = merge_stack.back().second;
	      merge_stack.pop_back();

	      // Parent of component
	      int dn = merge_stack.back().first;
	      bool dd = merge_stack.back().second;
	      merge_stack.pop_back();

	      Node parent = order_list[dn];

	      // Erasing from merge_roots
	      merge_roots[parent].pop_front();
	      
	      Node child = order_list[cn - order_list.size()];

	      // Erasing from child_lists
	      if (child_lists[child].prev != INVALID) {
		child_lists[child_lists[child].prev].next =
		  child_lists[child].next;
	      } else {
		child_lists[parent].first = child_lists[child].next;
	      }
	      
	      if (child_lists[child].next != INVALID) {
		child_lists[child_lists[child].next].prev =
		  child_lists[child].prev;
	      }

	      // Merging edges + flipping
	      Edge de = node_data[dn].first;
	      Edge ce = node_data[cn].first;

	      flip_map[order_list[cn - order_list.size()]] = cd != dd;
	      if (cd != dd) {
		std::swap(edge_lists[ce].prev, edge_lists[ce].next);
		ce = edge_lists[ce].prev;
		std::swap(edge_lists[ce].prev, edge_lists[ce].next);
	      }

	      {
		Edge dne = edge_lists[de].next; 
		Edge cne = edge_lists[ce].next; 

		edge_lists[de].next = cne;
		edge_lists[ce].next = dne;
	      
		edge_lists[dne].prev = ce;
		edge_lists[cne].prev = de;
	      }
	      	      
	      if (dd) {
		node_data[dn].first = ce;
	      }
	      
	      // Merging external faces
	      {
		int en = cn;
		cn = cd ? node_data[cn].prev : node_data[cn].next;
		cd = node_data[cn].next == en;

 		if (node_data[cn].prev == node_data[cn].next && 
		    node_data[cn].inverted) {
 		  cd = !cd;
 		}
	      }

	      if (cd) node_data[cn].next = dn; else node_data[cn].prev = dn;
	      if (dd) node_data[dn].prev = cn; else node_data[dn].next = cn;

	    }

	    bool d = pn == node_data[n].prev;

	    if (node_data[n].prev == node_data[n].next && 
		node_data[n].inverted) {
	      d = !d;
	    }

	    // Add new edge
	    {
	      Edge edge = embed_edge[node];
	      Edge re = node_data[rn].first;

	      edge_lists[edge_lists[re].next].prev = edge;
	      edge_lists[edge].next = edge_lists[re].next;
	      edge_lists[edge].prev = re;
	      edge_lists[re].next = edge;

	      if (!rd) {
		node_data[rn].first = edge;
	      }

	      Edge rev = _ugraph.oppositeEdge(edge);
	      Edge e = node_data[n].first;

	      edge_lists[edge_lists[e].next].prev = rev;
	      edge_lists[rev].next = edge_lists[e].next;
	      edge_lists[rev].prev = e;
	      edge_lists[e].next = rev;

	      if (d) {
		node_data[n].first = rev;
	      }
	      
	    }

	    // Embedding edge into external face
	    if (rd) node_data[rn].next = n; else node_data[rn].prev = n;
	    if (d) node_data[n].prev = rn; else node_data[n].next = rn;
	    pn = rn;

	    embed_edge[order_list[n]] = INVALID;
	  }

	  if (!merge_roots[node].empty()) {

	    bool d = pn == node_data[n].prev;
	    if (node_data[n].prev == node_data[n].next && 
		node_data[n].inverted) {
	      d = !d;
	    }

	    merge_stack.push_back(std::make_pair(n, d));

	    int rn = merge_roots[node].front();
	    
	    int xn = node_data[rn].next;
	    Node xnode = order_list[xn];
	    
	    int yn = node_data[rn].prev;
	    Node ynode = order_list[yn];
	    
	    bool rd;
	    if (!external(xnode, rorder, child_lists, ancestor_map, low_map)) {
	      rd = true;
	    } else if (!external(ynode, rorder, child_lists, 
				 ancestor_map, low_map)) {
	      rd = false;
	    } else if (pertinent(xnode, embed_edge, merge_roots)) {
	      rd = true;
	    } else {
	      rd = false;
	    }
	    
	    merge_stack.push_back(std::make_pair(rn, rd));
	    
	    pn = rn;
	    n = rd ? xn : yn;	      
	    	    
	  } else if (!external(node, rorder, child_lists,
			       ancestor_map, low_map)) {
	    int nn = (node_data[n].next != pn ? 
		      node_data[n].next : node_data[n].prev);

	    bool nd = n == node_data[nn].prev;

	    if (nd) node_data[nn].prev = pn;
	    else node_data[nn].next = pn; 

	    if (n == node_data[pn].prev) node_data[pn].prev = nn;
	    else node_data[pn].next = nn;

	    node_data[nn].inverted = 
	      (node_data[nn].prev == node_data[nn].next && nd != rd);

	    n = nn;
	  }
	  else break;
	  
	}

	if (!merge_stack.empty() || n == rn) {
	  break;
	}
      }
    }

    void initFace(const Node& node, EdgeLists& edge_lists,
		   NodeData& node_data, const PredMap& pred_map,
		   const OrderMap& order_map, const OrderList& order_list) {
      int n = order_map[node];
      int rn = n + order_list.size();
      
      node_data[n].next = node_data[n].prev = rn;
      node_data[rn].next = node_data[rn].prev = n;
      
      node_data[n].visited = order_list.size();
      node_data[rn].visited = order_list.size();

      node_data[n].inverted = false;
      node_data[rn].inverted = false;

      Edge edge = pred_map[node];
      Edge rev = _ugraph.oppositeEdge(edge);

      node_data[rn].first = edge;
      node_data[n].first = rev;

      edge_lists[edge].prev = edge;
      edge_lists[edge].next = edge;

      edge_lists[rev].prev = rev;
      edge_lists[rev].next = rev;

    }

    void mergeRemainingFaces(const Node& node, NodeData& node_data,
			     OrderList& order_list, OrderMap& order_map,
			     ChildLists& child_lists, EdgeLists& edge_lists) {
      while (child_lists[node].first != INVALID) {
	int dd = order_map[node];
	Node child = child_lists[node].first; 
	int cd = order_map[child] + order_list.size();
	child_lists[node].first = child_lists[child].next;

	Edge de = node_data[dd].first;
	Edge ce = node_data[cd].first;

	if (de != INVALID) {
	  Edge dne = edge_lists[de].next; 
	  Edge cne = edge_lists[ce].next; 
	  
	  edge_lists[de].next = cne;
	  edge_lists[ce].next = dne;
	  
	  edge_lists[dne].prev = ce;
	  edge_lists[cne].prev = de;
	}
	
	node_data[dd].first = ce;

      }
    }

    void storeEmbedding(const Node& node, NodeData& node_data,
			OrderMap& order_map, PredMap& pred_map,
			EdgeLists& edge_lists, FlipMap& flip_map) {

      if (node_data[order_map[node]].first == INVALID) return;

      if (pred_map[node] != INVALID) {
	Node source = _ugraph.source(pred_map[node]);
	flip_map[node] = flip_map[node] != flip_map[source];
      }
      
      Edge first = node_data[order_map[node]].first;
      Edge prev = first;

      Edge edge = flip_map[node] ?
	edge_lists[prev].prev : edge_lists[prev].next;

      _embedding[prev] = edge;
      
      while (edge != first) {
	Edge next = edge_lists[edge].prev == prev ?
	  edge_lists[edge].next : edge_lists[edge].prev;
	prev = edge; edge = next;
	_embedding[prev] = edge;
      }
    }


    bool external(const Node& node, int rorder,
		  ChildLists& child_lists, AncestorMap& ancestor_map, 
		  LowMap& low_map) {
      Node child = child_lists[node].first;

      if (child != INVALID) {
	if (low_map[child] < rorder) return true;
      }

      if (ancestor_map[node] < rorder) return true;

      return false;
    }

    bool pertinent(const Node& node, const EmbedEdge& embed_edge,
		   const MergeRoots& merge_roots) {
      return !merge_roots[node].empty() || embed_edge[node] != INVALID;
    }

    int lowPoint(const Node& node, OrderMap& order_map, ChildLists& child_lists,
		 AncestorMap& ancestor_map, LowMap& low_map) {
      int low_point;

      Node child = child_lists[node].first;

      if (child != INVALID) {
	low_point = low_map[child];
      } else {
	low_point = order_map[node];
      }

      if (low_point > ancestor_map[node]) {
	low_point = ancestor_map[node];
      }

      return low_point;
    }

    int findComponentRoot(Node root, Node node, ChildLists& child_lists, 
			  OrderMap& order_map, OrderList& order_list) {

      int order = order_map[root];
      int norder = order_map[node];

      Node child = child_lists[root].first;
      while (child != INVALID) {
	int corder = order_map[child];
	if (corder > order && corder < norder) {
	  order = corder;
	}
	child = child_lists[child].next;
      }
      return order + order_list.size();
    }

    Node findPertinent(Node node, OrderMap& order_map, NodeData& node_data,
		       EmbedEdge& embed_edge, MergeRoots& merge_roots) {
      Node wnode =_ugraph.target(node_data[order_map[node]].first);
      while (!pertinent(wnode, embed_edge, merge_roots)) {
	wnode = _ugraph.target(node_data[order_map[wnode]].first);
      }
      return wnode;
    }


    Node findExternal(Node node, int rorder, OrderMap& order_map, 
		      ChildLists& child_lists, AncestorMap& ancestor_map,
		      LowMap& low_map, NodeData& node_data) {
      Node wnode =_ugraph.target(node_data[order_map[node]].first);
      while (!external(wnode, rorder, child_lists, ancestor_map, low_map)) {
	wnode = _ugraph.target(node_data[order_map[wnode]].first);
      }
      return wnode;
    }

    void markCommonPath(Node node, int rorder, Node& wnode, Node& znode, 
			OrderList& order_list, OrderMap& order_map, 
			NodeData& node_data, EdgeLists& edge_lists, 
			EmbedEdge& embed_edge, MergeRoots& merge_roots, 
			ChildLists& child_lists, AncestorMap& ancestor_map, 
			LowMap& low_map) {
      
      Node cnode = node;
      Node pred = INVALID;
      
      while (true) {

	bool pert = pertinent(cnode, embed_edge, merge_roots);
	bool ext = external(cnode, rorder, child_lists, ancestor_map, low_map);

	if (pert && ext) {
	  if (!merge_roots[cnode].empty()) {
	    int cn = merge_roots[cnode].back();
	    
	    if (low_map[order_list[cn - order_list.size()]] < rorder) {
	      Edge edge = node_data[cn].first;
	      _kuratowski.set(edge, true);
	      
	      pred = cnode;
	      cnode = _ugraph.target(edge);
	    
	      continue;
	    }
	  }
	  wnode = znode = cnode;
	  return;

	} else if (pert) {
	  wnode = cnode;
	  
	  while (!external(cnode, rorder, child_lists, ancestor_map, low_map)) {
	    Edge edge = node_data[order_map[cnode]].first;
	  
	    if (_ugraph.target(edge) == pred) {
	      edge = edge_lists[edge].next;
	    }
	    _kuratowski.set(edge, true);
	    
	    Node next = _ugraph.target(edge);
	    pred = cnode; cnode = next;
	  }
	  
	  znode = cnode;
	  return;

	} else if (ext) {
	  znode = cnode;
	  
	  while (!pertinent(cnode, embed_edge, merge_roots)) {
	    Edge edge = node_data[order_map[cnode]].first;
	  
	    if (_ugraph.target(edge) == pred) {
	      edge = edge_lists[edge].next;
	    }
	    _kuratowski.set(edge, true);
	    
	    Node next = _ugraph.target(edge);
	    pred = cnode; cnode = next;
	  }
	  
	  wnode = cnode;
	  return;
	  
	} else {
	  Edge edge = node_data[order_map[cnode]].first;
	  
	  if (_ugraph.target(edge) == pred) {
	    edge = edge_lists[edge].next;
	  }
	  _kuratowski.set(edge, true);

	  Node next = _ugraph.target(edge);
	  pred = cnode; cnode = next;
	}
	
      }

    }

    void orientComponent(Node root, int rn, OrderMap& order_map,
			 PredMap& pred_map, NodeData& node_data, 
			 EdgeLists& edge_lists, FlipMap& flip_map, 
			 TypeMap& type_map) {
      node_data[order_map[root]].first = node_data[rn].first;
      type_map[root] = 1;

      std::vector<Node> st, qu;

      st.push_back(root);
      while (!st.empty()) {
	Node node = st.back();
	st.pop_back();
	qu.push_back(node);
	
	Edge edge = node_data[order_map[node]].first;
	
	if (type_map[_ugraph.target(edge)] == 0) {
	  st.push_back(_ugraph.target(edge));
	  type_map[_ugraph.target(edge)] = 1;
	} 
	
	Edge last = edge, pred = edge;
	edge = edge_lists[edge].next;
	while (edge != last) {

	  if (type_map[_ugraph.target(edge)] == 0) {
	    st.push_back(_ugraph.target(edge));
	    type_map[_ugraph.target(edge)] = 1;
	  } 
	  
	  Edge next = edge_lists[edge].next != pred ? 
	    edge_lists[edge].next : edge_lists[edge].prev;
	  pred = edge; edge = next;
	}

      }

      type_map[root] = 2;
      flip_map[root] = false;

      for (int i = 1; i < int(qu.size()); ++i) {

	Node node = qu[i];

	while (type_map[node] != 2) {
	  st.push_back(node);
	  type_map[node] = 2;
	  node = _ugraph.source(pred_map[node]);
	}

	bool flip = flip_map[node];

	while (!st.empty()) {
	  node = st.back();
	  st.pop_back();
	  
	  flip_map[node] = flip != flip_map[node];
	  flip = flip_map[node];

	  if (flip) {
	    Edge edge = node_data[order_map[node]].first;
	    std::swap(edge_lists[edge].prev, edge_lists[edge].next);
	    edge = edge_lists[edge].prev;
	    std::swap(edge_lists[edge].prev, edge_lists[edge].next);
	    node_data[order_map[node]].first = edge;
	  }
	}
      }

      for (int i = 0; i < int(qu.size()); ++i) {

	Edge edge = node_data[order_map[qu[i]]].first;
	Edge last = edge, pred = edge;

	edge = edge_lists[edge].next;
	while (edge != last) {

	  if (edge_lists[edge].next == pred) {
	    std::swap(edge_lists[edge].next, edge_lists[edge].prev);
	  } 
	  pred = edge; edge = edge_lists[edge].next;
	}
	
      }
    }

    void setFaceFlags(Node root, Node wnode, Node ynode, Node xnode,
		      OrderMap& order_map, NodeData& node_data, 
		      TypeMap& type_map) {
      Node node = _ugraph.target(node_data[order_map[root]].first);

      while (node != ynode) {
	type_map[node] = HIGHY;
	node = _ugraph.target(node_data[order_map[node]].first);
      }

      while (node != wnode) {
	type_map[node] = LOWY;
	node = _ugraph.target(node_data[order_map[node]].first);
      }
      
      node = _ugraph.target(node_data[order_map[wnode]].first);
      
      while (node != xnode) {
	type_map[node] = LOWX;
	node = _ugraph.target(node_data[order_map[node]].first);
      }
      type_map[node] = LOWX;

      node = _ugraph.target(node_data[order_map[xnode]].first);
      while (node != root) {
	type_map[node] = HIGHX;
	node = _ugraph.target(node_data[order_map[node]].first);
      }

      type_map[wnode] = PERTINENT;
      type_map[root] = ROOT;
    }

    void findInternalPath(std::vector<Edge>& ipath,
			  Node wnode, Node root, TypeMap& type_map, 
			  OrderMap& order_map, NodeData& node_data, 
			  EdgeLists& edge_lists) {
      std::vector<Edge> st;

      Node node = wnode;
      
      while (node != root) {
	Edge edge = edge_lists[node_data[order_map[node]].first].next;
	st.push_back(edge);
	node = _ugraph.target(edge);
      }
      
      while (true) {
	Edge edge = st.back();
	if (type_map[_ugraph.target(edge)] == LOWX ||
	    type_map[_ugraph.target(edge)] == HIGHX) {
	  break;
	}
	if (type_map[_ugraph.target(edge)] == 2) {
	  type_map[_ugraph.target(edge)] = 3;
	  
	  edge = edge_lists[_ugraph.oppositeEdge(edge)].next;
	  st.push_back(edge);
	} else {
	  st.pop_back();
	  edge = edge_lists[edge].next;
	  
	  while (_ugraph.oppositeEdge(edge) == st.back()) {
	    edge = st.back();
	    st.pop_back();
	    edge = edge_lists[edge].next;
	  }
	  st.push_back(edge);
	}
      }
      
      for (int i = 0; i < int(st.size()); ++i) {
	if (type_map[_ugraph.target(st[i])] != LOWY &&
	    type_map[_ugraph.target(st[i])] != HIGHY) {
	  for (; i < int(st.size()); ++i) {
	    ipath.push_back(st[i]);
	  }
	}
      }
    }

    void setInternalFlags(std::vector<Edge>& ipath, TypeMap& type_map) {
      for (int i = 1; i < int(ipath.size()); ++i) {
	type_map[_ugraph.source(ipath[i])] = INTERNAL;
      }
    }

    void findPilePath(std::vector<Edge>& ppath,
		      Node root, TypeMap& type_map, OrderMap& order_map, 
		      NodeData& node_data, EdgeLists& edge_lists) {
      std::vector<Edge> st;

      st.push_back(_ugraph.oppositeEdge(node_data[order_map[root]].first));
      st.push_back(node_data[order_map[root]].first);
      
      while (st.size() > 1) {
	Edge edge = st.back();
	if (type_map[_ugraph.target(edge)] == INTERNAL) {
	  break;
	}
	if (type_map[_ugraph.target(edge)] == 3) {
	  type_map[_ugraph.target(edge)] = 4;
	  
	  edge = edge_lists[_ugraph.oppositeEdge(edge)].next;
	  st.push_back(edge);
	} else {
	  st.pop_back();
	  edge = edge_lists[edge].next;
	  
	  while (!st.empty() && _ugraph.oppositeEdge(edge) == st.back()) {
	    edge = st.back();
	    st.pop_back();
	    edge = edge_lists[edge].next;
	  }
	  st.push_back(edge);
	}
      }
      
      for (int i = 1; i < int(st.size()); ++i) {
	ppath.push_back(st[i]);
      }
    }


    int markExternalPath(Node node, OrderMap& order_map,
			 ChildLists& child_lists, PredMap& pred_map,
			 AncestorMap& ancestor_map, LowMap& low_map) {
      int lp = lowPoint(node, order_map, child_lists,
			ancestor_map, low_map);
      
      if (ancestor_map[node] != lp) {
	node = child_lists[node].first;
	_kuratowski[pred_map[node]] = true;

	while (ancestor_map[node] != lp) {
	  for (OutEdgeIt e(_ugraph, node); e != INVALID; ++e) {
	    Node tnode = _ugraph.target(e); 
	    if (order_map[tnode] > order_map[node] && low_map[tnode] == lp) {
	      node = tnode;
	      _kuratowski[e] = true;
	      break;
	    }
	  }
	}
      }

      for (OutEdgeIt e(_ugraph, node); e != INVALID; ++e) {
	if (order_map[_ugraph.target(e)] == lp) {
	  _kuratowski[e] = true;
	  break;
	}
      }
      
      return lp;
    }

    void markPertinentPath(Node node, OrderMap& order_map, 
			   NodeData& node_data, EdgeLists& edge_lists,
			   EmbedEdge& embed_edge, MergeRoots& merge_roots) {
      while (embed_edge[node] == INVALID) {
	int n = merge_roots[node].front();
	Edge edge = node_data[n].first;

	_kuratowski.set(edge, true);
	
	Node pred = node;
	node = _ugraph.target(edge);
	while (!pertinent(node, embed_edge, merge_roots)) {
	  edge = node_data[order_map[node]].first;
	  if (_ugraph.target(edge) == pred) {
	    edge = edge_lists[edge].next;
	  }
	  _kuratowski.set(edge, true);
	  pred = node;
	  node = _ugraph.target(edge);
	}
      }
      _kuratowski.set(embed_edge[node], true);
    } 

    void markPredPath(Node node, Node snode, PredMap& pred_map) {
      while (node != snode) {
	_kuratowski.set(pred_map[node], true);
	node = _ugraph.source(pred_map[node]);
      }
    }

    void markFacePath(Node ynode, Node xnode, 
		      OrderMap& order_map, NodeData& node_data) {
      Edge edge = node_data[order_map[ynode]].first;
      Node node = _ugraph.target(edge);
      _kuratowski.set(edge, true);
	
      while (node != xnode) {
	edge = node_data[order_map[node]].first;
	_kuratowski.set(edge, true);
	node = _ugraph.target(edge);
      }
    }

    void markInternalPath(std::vector<Edge>& path) {
      for (int i = 0; i < int(path.size()); ++i) {
	_kuratowski.set(path[i], true);
      }
    }

    void markPilePath(std::vector<Edge>& path) {
      for (int i = 0; i < int(path.size()); ++i) {
	_kuratowski.set(path[i], true);
      }
    }

    void isolateKuratowski(Edge edge, NodeData& node_data, 
			   EdgeLists& edge_lists, FlipMap& flip_map,
			   OrderMap& order_map, OrderList& order_list, 
			   PredMap& pred_map, ChildLists& child_lists,
			   AncestorMap& ancestor_map, LowMap& low_map, 
			   EmbedEdge& embed_edge, MergeRoots& merge_roots) {

      Node root = _ugraph.source(edge);
      Node enode = _ugraph.target(edge);

      int rorder = order_map[root];

      TypeMap type_map(_ugraph, 0);

      int rn = findComponentRoot(root, enode, child_lists, 
				 order_map, order_list);

      Node xnode = order_list[node_data[rn].next];
      Node ynode = order_list[node_data[rn].prev];

      // Minor-A
      {
	while (!merge_roots[xnode].empty() || !merge_roots[ynode].empty()) {
	  
	  if (!merge_roots[xnode].empty()) {
	    root = xnode;
	    rn = merge_roots[xnode].front();
	  } else {
	    root = ynode;
	    rn = merge_roots[ynode].front();
	  }
	  
	  xnode = order_list[node_data[rn].next];
	  ynode = order_list[node_data[rn].prev];
	}
	
	if (root != _ugraph.source(edge)) {
	  orientComponent(root, rn, order_map, pred_map, 
			  node_data, edge_lists, flip_map, type_map);
	  markFacePath(root, root, order_map, node_data);
	  int xlp = markExternalPath(xnode, order_map, child_lists, 
				     pred_map, ancestor_map, low_map);
	  int ylp = markExternalPath(ynode, order_map, child_lists, 
				     pred_map, ancestor_map, low_map);
	  markPredPath(root, order_list[xlp < ylp ? xlp : ylp], pred_map);
	  Node lwnode = findPertinent(ynode, order_map, node_data,
				     embed_edge, merge_roots);
	  
	  markPertinentPath(lwnode, order_map, node_data, edge_lists,
			    embed_edge, merge_roots);
	  
	  return;
	}
      }
      
      orientComponent(root, rn, order_map, pred_map, 
		      node_data, edge_lists, flip_map, type_map);

      Node wnode = findPertinent(ynode, order_map, node_data,
				 embed_edge, merge_roots);
      setFaceFlags(root, wnode, ynode, xnode, order_map, node_data, type_map);

      
      //Minor-B
      if (!merge_roots[wnode].empty()) {
	int cn = merge_roots[wnode].back();
	Node rep = order_list[cn - order_list.size()];
	if (low_map[rep] < rorder) {
	  markFacePath(root, root, order_map, node_data);
	  int xlp = markExternalPath(xnode, order_map, child_lists, 
				     pred_map, ancestor_map, low_map);
	  int ylp = markExternalPath(ynode, order_map, child_lists, 
				     pred_map, ancestor_map, low_map);

	  Node lwnode, lznode;
	  markCommonPath(wnode, rorder, lwnode, lznode, order_list, 
			 order_map, node_data, edge_lists, embed_edge, 
			 merge_roots, child_lists, ancestor_map, low_map);
	  	  
	  markPertinentPath(lwnode, order_map, node_data, edge_lists,
			    embed_edge, merge_roots);
	  int zlp = markExternalPath(lznode, order_map, child_lists, 
				     pred_map, ancestor_map, low_map);

	  int minlp = xlp < ylp ? xlp : ylp;
	  if (zlp < minlp) minlp = zlp;

	  int maxlp = xlp > ylp ? xlp : ylp;
	  if (zlp > maxlp) maxlp = zlp;
	  
	  markPredPath(order_list[maxlp], order_list[minlp], pred_map);
	  
	  return;
	}
      }

      Node pxnode, pynode;
      std::vector<Edge> ipath;
      findInternalPath(ipath, wnode, root, type_map, order_map,
		       node_data, edge_lists);
      setInternalFlags(ipath, type_map);
      pynode = _ugraph.source(ipath.front());
      pxnode = _ugraph.target(ipath.back());

      wnode = findPertinent(pynode, order_map, node_data,
			    embed_edge, merge_roots);
      
      // Minor-C
      {
	if (type_map[_ugraph.source(ipath.front())] == HIGHY) {
	  if (type_map[_ugraph.target(ipath.back())] == HIGHX) {
	    markFacePath(xnode, pxnode, order_map, node_data);
	  }
	  markFacePath(root, xnode, order_map, node_data);
	  markPertinentPath(wnode, order_map, node_data, edge_lists,
			    embed_edge, merge_roots);
	  markInternalPath(ipath);
	  int xlp = markExternalPath(xnode, order_map, child_lists, 
				     pred_map, ancestor_map, low_map);
	  int ylp = markExternalPath(ynode, order_map, child_lists, 
				     pred_map, ancestor_map, low_map);
	  markPredPath(root, order_list[xlp < ylp ? xlp : ylp], pred_map);
	  return;
	}

	if (type_map[_ugraph.target(ipath.back())] == HIGHX) {
	  markFacePath(ynode, root, order_map, node_data);
	  markPertinentPath(wnode, order_map, node_data, edge_lists,
			    embed_edge, merge_roots);
	  markInternalPath(ipath);
	  int xlp = markExternalPath(xnode, order_map, child_lists, 
				     pred_map, ancestor_map, low_map);
	  int ylp = markExternalPath(ynode, order_map, child_lists, 
				     pred_map, ancestor_map, low_map);
	  markPredPath(root, order_list[xlp < ylp ? xlp : ylp], pred_map);
	  return;
	}
      }

      std::vector<Edge> ppath;
      findPilePath(ppath, root, type_map, order_map, node_data, edge_lists);
      
      // Minor-D
      if (!ppath.empty()) {
	markFacePath(ynode, xnode, order_map, node_data);
	markPertinentPath(wnode, order_map, node_data, edge_lists,
			  embed_edge, merge_roots);
	markPilePath(ppath);
	markInternalPath(ipath);
	int xlp = markExternalPath(xnode, order_map, child_lists, 
				   pred_map, ancestor_map, low_map);
	int ylp = markExternalPath(ynode, order_map, child_lists, 
				   pred_map, ancestor_map, low_map);
	markPredPath(root, order_list[xlp < ylp ? xlp : ylp], pred_map);
	return;
      }

      // Minor-E*
      {

	if (!external(wnode, rorder, child_lists, ancestor_map, low_map)) {
	  Node znode = findExternal(pynode, rorder, order_map, 
				    child_lists, ancestor_map,
				    low_map, node_data);
	  
	  if (type_map[znode] == LOWY) {
	    markFacePath(root, xnode, order_map, node_data);
	    markPertinentPath(wnode, order_map, node_data, edge_lists,
			      embed_edge, merge_roots);
	    markInternalPath(ipath);
	    int xlp = markExternalPath(xnode, order_map, child_lists, 
				       pred_map, ancestor_map, low_map);
	    int zlp = markExternalPath(znode, order_map, child_lists, 
				       pred_map, ancestor_map, low_map);
	    markPredPath(root, order_list[xlp < zlp ? xlp : zlp], pred_map);
	  } else {
	    markFacePath(ynode, root, order_map, node_data);
	    markPertinentPath(wnode, order_map, node_data, edge_lists,
			      embed_edge, merge_roots);
	    markInternalPath(ipath);
	    int ylp = markExternalPath(ynode, order_map, child_lists, 
				       pred_map, ancestor_map, low_map);
	    int zlp = markExternalPath(znode, order_map, child_lists, 
				       pred_map, ancestor_map, low_map);
	    markPredPath(root, order_list[ylp < zlp ? ylp : zlp], pred_map);
	  }
	  return;
	}

	int xlp = markExternalPath(xnode, order_map, child_lists, 
				   pred_map, ancestor_map, low_map);
	int ylp = markExternalPath(ynode, order_map, child_lists, 
				   pred_map, ancestor_map, low_map);
	int wlp = markExternalPath(wnode, order_map, child_lists, 
				   pred_map, ancestor_map, low_map);

	if (wlp > xlp && wlp > ylp) {
	  markFacePath(root, root, order_map, node_data);
	  markPredPath(root, order_list[xlp < ylp ? xlp : ylp], pred_map);
	  return;
	}

	markInternalPath(ipath);
	markPertinentPath(wnode, order_map, node_data, edge_lists,
			  embed_edge, merge_roots);

	if (xlp > ylp && xlp > wlp) {
	  markFacePath(root, pynode, order_map, node_data);
	  markFacePath(wnode, xnode, order_map, node_data);
	  markPredPath(root, order_list[ylp < wlp ? ylp : wlp], pred_map);
	  return;
	}

	if (ylp > xlp && ylp > wlp) {
	  markFacePath(pxnode, root, order_map, node_data);
	  markFacePath(ynode, wnode, order_map, node_data);
	  markPredPath(root, order_list[xlp < wlp ? xlp : wlp], pred_map);
	  return;
	}

	if (pynode != ynode) {
	  markFacePath(pxnode, wnode, order_map, node_data);

	  int minlp = xlp < ylp ? xlp : ylp;
	  if (wlp < minlp) minlp = wlp;

	  int maxlp = xlp > ylp ? xlp : ylp;
	  if (wlp > maxlp) maxlp = wlp;
	  
	  markPredPath(order_list[maxlp], order_list[minlp], pred_map);
	  return;
	}

	if (pxnode != xnode) {
	  markFacePath(wnode, pynode, order_map, node_data);

	  int minlp = xlp < ylp ? xlp : ylp;
	  if (wlp < minlp) minlp = wlp;

	  int maxlp = xlp > ylp ? xlp : ylp;
	  if (wlp > maxlp) maxlp = wlp;
	  
	  markPredPath(order_list[maxlp], order_list[minlp], pred_map);
	  return;
	}

	markFacePath(root, root, order_map, node_data);
	int minlp = xlp < ylp ? xlp : ylp;
	if (wlp < minlp) minlp = wlp;
	markPredPath(root, order_list[minlp], pred_map);
	return;
      }
      
    }

  };

  namespace _planarity_bits {

    template <typename UGraph, typename EmbeddingMap>
    void makeConnected(UGraph& ugraph, EmbeddingMap& embedding) {
      DfsVisitor<UGraph> null_visitor;
      DfsVisit<UGraph, DfsVisitor<UGraph> > dfs(ugraph, null_visitor);
      dfs.init();

      typename UGraph::Node u = INVALID;
      for (typename UGraph::NodeIt n(ugraph); n != INVALID; ++n) {
	if (!dfs.reached(n)) {
	  dfs.addSource(n);
	  dfs.start();
	  if (u == INVALID) {
	    u = n;
	  } else {
	    typename UGraph::Node v = n;
	    
	    typename UGraph::Edge ue = typename UGraph::OutEdgeIt(ugraph, u);
	    typename UGraph::Edge ve = typename UGraph::OutEdgeIt(ugraph, v);

	    typename UGraph::Edge e = ugraph.direct(ugraph.addEdge(u, v), true);
	    
	    if (ue != INVALID) {
	      embedding[e] = embedding[ue];
	      embedding[ue] = e;
	    } else {
	      embedding[e] = e;
	    }

	    if (ve != INVALID) {
	      embedding[ugraph.oppositeEdge(e)] = embedding[ve];
	      embedding[ve] = ugraph.oppositeEdge(e);
	    } else {
	      embedding[ugraph.oppositeEdge(e)] = ugraph.oppositeEdge(e);
	    }
	  }
	}
      }
    }

    template <typename UGraph, typename EmbeddingMap>
    void makeBiNodeConnected(UGraph& ugraph, EmbeddingMap& embedding) {
      typename UGraph::template EdgeMap<bool> processed(ugraph);

      std::vector<typename UGraph::Edge> edges;
      for (typename UGraph::EdgeIt e(ugraph); e != INVALID; ++e) {
	edges.push_back(e);
      }

      IterableBoolMap<UGraph, typename UGraph::Node> visited(ugraph, false);
      
      for (int i = 0; i < int(edges.size()); ++i) {
	typename UGraph::Edge pp = edges[i];
	if (processed[pp]) continue;

	typename UGraph::Edge e = embedding[ugraph.oppositeEdge(pp)];
	processed[e] = true;
	visited.set(ugraph.source(e), true);
	
	typename UGraph::Edge p = e, l = e;
	e = embedding[ugraph.oppositeEdge(e)];
	
	while (e != l) {
	  processed[e] = true;

	  if (visited[ugraph.source(e)]) {
	    
	    typename UGraph::Edge n = 
	      ugraph.direct(ugraph.addEdge(ugraph.source(p), 
					   ugraph.target(e)), true);
	    embedding[n] = p;
	    embedding[ugraph.oppositeEdge(pp)] = n;

	    embedding[ugraph.oppositeEdge(n)] = 
	      embedding[ugraph.oppositeEdge(e)];
	    embedding[ugraph.oppositeEdge(e)] = 
	      ugraph.oppositeEdge(n);
	    
	    p = n;
	    e = embedding[ugraph.oppositeEdge(n)];
	  } else {
	    visited.set(ugraph.source(e), true);
	    pp = p;
	    p = e;
	    e = embedding[ugraph.oppositeEdge(e)];
	  }
	}
	visited.setAll(false);
      }
    }
    
    
    template <typename UGraph, typename EmbeddingMap>
    void makeMaxPlanar(UGraph& ugraph, EmbeddingMap& embedding) {
      
      typename UGraph::template NodeMap<int> degree(ugraph);

      for (typename UGraph::NodeIt n(ugraph); n != INVALID; ++n) {
	degree[n] = countIncEdges(ugraph, n);
      }

      typename UGraph::template EdgeMap<bool> processed(ugraph);
      IterableBoolMap<UGraph, typename UGraph::Node> visited(ugraph, false);

      std::vector<typename UGraph::Edge> edges;
      for (typename UGraph::EdgeIt e(ugraph); e != INVALID; ++e) {
	edges.push_back(e);
      }

      for (int i = 0; i < int(edges.size()); ++i) {
	typename UGraph::Edge e = edges[i];

	if (processed[e]) continue;
	processed[e] = true;

	typename UGraph::Edge mine = e;
	int mind = degree[ugraph.source(e)];

	int face_size = 1;	

	typename UGraph::Edge l = e;
	e = embedding[ugraph.oppositeEdge(e)];
	while (l != e) {
	  processed[e] = true;

	  ++face_size;

	  if (degree[ugraph.source(e)] < mind) {
	    mine = e;
	    mind = degree[ugraph.source(e)];
	  }
	  
	  e = embedding[ugraph.oppositeEdge(e)];	  
	}
	
	if (face_size < 4) {
	  continue;
	}

	typename UGraph::Node s = ugraph.source(mine);
	for (typename UGraph::OutEdgeIt e(ugraph, s); e != INVALID; ++e) {
	  visited.set(ugraph.target(e), true); 
	}

	typename UGraph::Edge oppe = INVALID;

	e = embedding[ugraph.oppositeEdge(mine)];
	e = embedding[ugraph.oppositeEdge(e)];
	while (ugraph.target(e) != s) {
	  if (visited[ugraph.source(e)]) {
	    oppe = e;
	    break;
	  }
	  e = embedding[ugraph.oppositeEdge(e)];
	}
	visited.setAll(false);
	
	if (oppe == INVALID) {

	  e = embedding[ugraph.oppositeEdge(mine)];
	  typename UGraph::Edge pn = mine, p = e;

	  e = embedding[ugraph.oppositeEdge(e)];
	  while (ugraph.target(e) != s) {
	    typename UGraph::Edge n = 
	      ugraph.direct(ugraph.addEdge(s, ugraph.source(e)), true);

	    embedding[n] = pn;
	    embedding[ugraph.oppositeEdge(n)] = e;
	    embedding[ugraph.oppositeEdge(p)] = ugraph.oppositeEdge(n);

	    pn = n;
	    
	    p = e;
	    e = embedding[ugraph.oppositeEdge(e)];
	  }

	  embedding[ugraph.oppositeEdge(e)] = pn;

	} else {

	  mine = embedding[ugraph.oppositeEdge(mine)];
	  s = ugraph.source(mine);
	  oppe = embedding[ugraph.oppositeEdge(oppe)];
	  typename UGraph::Node t = ugraph.source(oppe);
	  
	  typename UGraph::Edge ce = ugraph.direct(ugraph.addEdge(s, t), true);
	  embedding[ce] = mine;
	  embedding[ugraph.oppositeEdge(ce)] = oppe;

	  typename UGraph::Edge pn = ce, p = oppe;	  
	  e = embedding[ugraph.oppositeEdge(oppe)];
	  while (ugraph.target(e) != s) {
	    typename UGraph::Edge n = 
	      ugraph.direct(ugraph.addEdge(s, ugraph.source(e)), true);

	    embedding[n] = pn;
	    embedding[ugraph.oppositeEdge(n)] = e;
	    embedding[ugraph.oppositeEdge(p)] = ugraph.oppositeEdge(n);

	    pn = n;
	    
	    p = e;
	    e = embedding[ugraph.oppositeEdge(e)];
	    
	  }
	  embedding[ugraph.oppositeEdge(e)] = pn;

	  pn = ugraph.oppositeEdge(ce), p = mine;	  
	  e = embedding[ugraph.oppositeEdge(mine)];
	  while (ugraph.target(e) != t) {
	    typename UGraph::Edge n = 
	      ugraph.direct(ugraph.addEdge(t, ugraph.source(e)), true);

	    embedding[n] = pn;
	    embedding[ugraph.oppositeEdge(n)] = e;
	    embedding[ugraph.oppositeEdge(p)] = ugraph.oppositeEdge(n);

	    pn = n;
	    
	    p = e;
	    e = embedding[ugraph.oppositeEdge(e)];
	    
	  }
	  embedding[ugraph.oppositeEdge(e)] = pn;
	}
      }
    }

  }

  /// \ingroup planar
  ///
  /// \brief Schnyder's planar drawing algorithms
  ///
  /// The planar drawing algorithm calculates location for each node
  /// in the plane, which coordinates satisfies that if each edge is
  /// represented with a straight line then the edges will not
  /// intersect each other.
  ///
  /// Scnyder's algorithm embeds the graph on \c (n-2,n-2) size grid,
  /// ie. each node will be located in the \c [0,n-2]x[0,n-2] square.
  /// The time complexity of the algorithm is O(n).
  template <typename UGraph>
  class PlanarDrawing {
  public:

    UGRAPH_TYPEDEFS(typename UGraph);

    /// \brief The point type for store coordinates
    typedef dim2::Point<int> Point;
    /// \brief The map type for store coordinates
    typedef typename UGraph::template NodeMap<Point> PointMap;


    /// \brief Constructor
    ///
    /// Constructor
    /// \pre The ugraph should be simple, ie. loop and parallel edge free. 
    PlanarDrawing(const UGraph& ugraph)
      : _ugraph(ugraph), _point_map(ugraph) {}

  private:

    template <typename AuxUGraph, typename AuxEmbeddingMap>
    void drawing(const AuxUGraph& ugraph, 
		 const AuxEmbeddingMap& next, 
		 PointMap& point_map) {
      UGRAPH_TYPEDEFS(typename AuxUGraph);

      typename AuxUGraph::template EdgeMap<Edge> prev(ugraph);

      for (NodeIt n(ugraph); n != INVALID; ++n) {
	Edge e = OutEdgeIt(ugraph, n);
	
	Edge p = e, l = e;
	
	e = next[e];
	while (e != l) {
	  prev[e] = p;
	  p = e;
	  e = next[e];
	}
	prev[e] = p;
      }

      Node anode, bnode, cnode;

      {
	Edge e = EdgeIt(ugraph);
	anode = ugraph.source(e);
	bnode = ugraph.target(e);
	cnode = ugraph.target(next[ugraph.oppositeEdge(e)]);
      }
      
      IterableBoolMap<AuxUGraph, Node> proper(ugraph, false);
      typename AuxUGraph::template NodeMap<int> conn(ugraph, -1);

      conn[anode] = conn[bnode] = -2;      
      {
	for (OutEdgeIt e(ugraph, anode); e != INVALID; ++e) {
	  Node m = ugraph.target(e);
	  if (conn[m] == -1) {
	    conn[m] = 1;
	  }
	}
	conn[cnode] = 2;

	for (OutEdgeIt e(ugraph, bnode); e != INVALID; ++e) {
	  Node m = ugraph.target(e);
	  if (conn[m] == -1) {
	    conn[m] = 1;
	  } else if (conn[m] != -2) {
	    conn[m] += 1;	    
	    Edge pe = ugraph.oppositeEdge(e);
	    if (conn[ugraph.target(next[pe])] == -2) {
	      conn[m] -= 1;
	    }
	    if (conn[ugraph.target(prev[pe])] == -2) {
	      conn[m] -= 1;
	    }

	    proper.set(m, conn[m] == 1);
	  }
	}
      }
      

      typename AuxUGraph::template EdgeMap<int> angle(ugraph, -1);

      while (proper.trueNum() != 0) {
	Node n = typename IterableBoolMap<AuxUGraph, Node>::TrueIt(proper);
	proper.set(n, false);
	conn[n] = -2;

	for (OutEdgeIt e(ugraph, n); e != INVALID; ++e) {
	  Node m = ugraph.target(e);
	  if (conn[m] == -1) {
	    conn[m] = 1;
	  } else if (conn[m] != -2) {
	    conn[m] += 1;	    
	    Edge pe = ugraph.oppositeEdge(e);
	    if (conn[ugraph.target(next[pe])] == -2) {
	      conn[m] -= 1;
	    }
	    if (conn[ugraph.target(prev[pe])] == -2) {
	      conn[m] -= 1;
	    }

	    proper.set(m, conn[m] == 1);
	  }
	}

	{
	  Edge e = OutEdgeIt(ugraph, n);
	  Edge p = e, l = e;
	  
	  e = next[e];
	  while (e != l) {
	    
	    if (conn[ugraph.target(e)] == -2 && conn[ugraph.target(p)] == -2) {
	      Edge f = e;
	      angle[f] = 0;
	      f = next[ugraph.oppositeEdge(f)];
	      angle[f] = 1;
	      f = next[ugraph.oppositeEdge(f)];
	      angle[f] = 2;
	    }
	    
	    p = e;
	    e = next[e];
	  }
	  
	  if (conn[ugraph.target(e)] == -2 && conn[ugraph.target(p)] == -2) {
	    Edge f = e;
	    angle[f] = 0;
	    f = next[ugraph.oppositeEdge(f)];
	    angle[f] = 1;
	    f = next[ugraph.oppositeEdge(f)];
	    angle[f] = 2;
	  }
	}
      }

      typename AuxUGraph::template NodeMap<Node> apred(ugraph, INVALID);
      typename AuxUGraph::template NodeMap<Node> bpred(ugraph, INVALID);
      typename AuxUGraph::template NodeMap<Node> cpred(ugraph, INVALID);

      typename AuxUGraph::template NodeMap<int> apredid(ugraph, -1);
      typename AuxUGraph::template NodeMap<int> bpredid(ugraph, -1);
      typename AuxUGraph::template NodeMap<int> cpredid(ugraph, -1);

      for (EdgeIt e(ugraph); e != INVALID; ++e) {
	if (angle[e] == angle[next[e]]) {
	  switch (angle[e]) {
	  case 2:
	    apred[ugraph.target(e)] = ugraph.source(e);
	    apredid[ugraph.target(e)] = ugraph.id(ugraph.source(e));
	    break;
	  case 1:
	    bpred[ugraph.target(e)] = ugraph.source(e);
	    bpredid[ugraph.target(e)] = ugraph.id(ugraph.source(e));
	    break;
	  case 0:
	    cpred[ugraph.target(e)] = ugraph.source(e);
	    cpredid[ugraph.target(e)] = ugraph.id(ugraph.source(e));
	    break;
	  }
	}
      }

      cpred[anode] = INVALID;
      cpred[bnode] = INVALID;

      std::vector<Node> aorder, border, corder; 

      {
	typename AuxUGraph::template NodeMap<bool> processed(ugraph, false);
	std::vector<Node> st;
	for (NodeIt n(ugraph); n != INVALID; ++n) {
	  if (!processed[n] && n != bnode && n != cnode) {
	    st.push_back(n);
	    processed[n] = true;
	    Node m = apred[n];
	    while (m != INVALID && !processed[m]) {
	      st.push_back(m);
	      processed[m] = true;
	      m = apred[m];
	    }
	    while (!st.empty()) {
	      aorder.push_back(st.back());
	      st.pop_back();
	    }
	  }
	}
      }

      {
	typename AuxUGraph::template NodeMap<bool> processed(ugraph, false);
	std::vector<Node> st;
	for (NodeIt n(ugraph); n != INVALID; ++n) {
	  if (!processed[n] && n != cnode && n != anode) {
	    st.push_back(n);
	    processed[n] = true;
	    Node m = bpred[n];
	    while (m != INVALID && !processed[m]) {
	      st.push_back(m);
	      processed[m] = true;
	      m = bpred[m];
	    }
	    while (!st.empty()) {
	      border.push_back(st.back());
	      st.pop_back();
	    }
	  }
	}
      }

      {
	typename AuxUGraph::template NodeMap<bool> processed(ugraph, false);
	std::vector<Node> st;
	for (NodeIt n(ugraph); n != INVALID; ++n) {
	  if (!processed[n] && n != anode && n != bnode) {
	    st.push_back(n);
	    processed[n] = true;
	    Node m = cpred[n];
	    while (m != INVALID && !processed[m]) {
	      st.push_back(m);
	      processed[m] = true;
	      m = cpred[m];
	    }
	    while (!st.empty()) {
	      corder.push_back(st.back());
	      st.pop_back();
	    }
	  }
	}
      }

      typename AuxUGraph::template NodeMap<int> atree(ugraph, 0);
      for (int i = aorder.size() - 1; i >= 0; --i) {
	Node n = aorder[i];
	atree[n] = 1;
	for (OutEdgeIt e(ugraph, n); e != INVALID; ++e) {
	  if (apred[ugraph.target(e)] == n) {
	    atree[n] += atree[ugraph.target(e)];
	  }
	} 
      }

      typename AuxUGraph::template NodeMap<int> btree(ugraph, 0);
      for (int i = border.size() - 1; i >= 0; --i) {
	Node n = border[i];
	btree[n] = 1;
	for (OutEdgeIt e(ugraph, n); e != INVALID; ++e) {
	  if (bpred[ugraph.target(e)] == n) {
	    btree[n] += btree[ugraph.target(e)];
	  }
	} 
      }
      
      typename AuxUGraph::template NodeMap<int> apath(ugraph, 0);
      apath[bnode] = apath[cnode] = 1;
      typename AuxUGraph::template NodeMap<int> apath_btree(ugraph, 0);
      apath_btree[bnode] = btree[bnode];
      for (int i = 1; i < int(aorder.size()); ++i) {
	Node n = aorder[i];
	apath[n] = apath[apred[n]] + 1;
	apath_btree[n] = btree[n] + apath_btree[apred[n]];
      }

      typename AuxUGraph::template NodeMap<int> bpath_atree(ugraph, 0);
      bpath_atree[anode] = atree[anode];
      for (int i = 1; i < int(border.size()); ++i) {
	Node n = border[i];
	bpath_atree[n] = atree[n] + bpath_atree[bpred[n]];
      }

      typename AuxUGraph::template NodeMap<int> cpath(ugraph, 0);
      cpath[anode] = cpath[bnode] = 1;
      typename AuxUGraph::template NodeMap<int> cpath_atree(ugraph, 0);
      cpath_atree[anode] = atree[anode];
      typename AuxUGraph::template NodeMap<int> cpath_btree(ugraph, 0);
      cpath_btree[bnode] = btree[bnode];
      for (int i = 1; i < int(corder.size()); ++i) {
	Node n = corder[i];
	cpath[n] = cpath[cpred[n]] + 1;
	cpath_atree[n] = atree[n] + cpath_atree[cpred[n]];
	cpath_btree[n] = btree[n] + cpath_btree[cpred[n]];
      }

      typename AuxUGraph::template NodeMap<int> third(ugraph);
      for (NodeIt n(ugraph); n != INVALID; ++n) {
	point_map[n].x = 
	  bpath_atree[n] + cpath_atree[n] - atree[n] - cpath[n] + 1;
	point_map[n].y = 
	  cpath_btree[n] + apath_btree[n] - btree[n] - apath[n] + 1;
      }
      
    }

  public:

    /// \brief Calculates the node locations
    ///
    /// This function calculates the node locations.
    bool run() {
      PlanarEmbedding<UGraph> pe(_ugraph);
      if (!pe.run()) return false;
      
      run(pe);
      return true;
    }

    /// \brief Calculates the node locations according to a
    /// combinatorical embedding
    ///
    /// This function calculates the node locations. The \c embedding
    /// parameter should contain a valid combinatorical embedding, ie.
    /// a valid cyclic order of the edges.
    template <typename EmbeddingMap>
    void run(const EmbeddingMap& embedding) {
      typedef SmartUEdgeSet<UGraph> AuxUGraph; 

      if (3 * countNodes(_ugraph) - 6 == countUEdges(_ugraph)) {
	drawing(_ugraph, embedding, _point_map);
	return;
      }

      AuxUGraph aux_ugraph(_ugraph);
      typename AuxUGraph::template EdgeMap<typename AuxUGraph::Edge> 
	aux_embedding(aux_ugraph);

      {

	typename UGraph::template UEdgeMap<typename AuxUGraph::UEdge> 
	  ref(_ugraph);
	
	for (UEdgeIt e(_ugraph); e != INVALID; ++e) {
	  ref[e] = aux_ugraph.addEdge(_ugraph.source(e), _ugraph.target(e));
	}

	for (UEdgeIt e(_ugraph); e != INVALID; ++e) {
	  Edge ee = embedding[_ugraph.direct(e, true)];
	  aux_embedding[aux_ugraph.direct(ref[e], true)] = 
	    aux_ugraph.direct(ref[ee], _ugraph.direction(ee));
	  ee = embedding[_ugraph.direct(e, false)];
	  aux_embedding[aux_ugraph.direct(ref[e], false)] = 
	    aux_ugraph.direct(ref[ee], _ugraph.direction(ee));
	}
      }
      _planarity_bits::makeConnected(aux_ugraph, aux_embedding);
      _planarity_bits::makeBiNodeConnected(aux_ugraph, aux_embedding);
      _planarity_bits::makeMaxPlanar(aux_ugraph, aux_embedding);
      drawing(aux_ugraph, aux_embedding, _point_map);
    }

    /// \brief The coordinate of the given node
    ///
    /// The coordinate of the given node.
    Point operator[](const Node& node) {
      return _point_map[node];
    }

    /// \brief Returns the grid embedding in a \e NodeMap.
    ///
    /// Returns the grid embedding in a \e NodeMap of \c dim2::Point<int> .
    const PointMap& coords() const {
      return _point_map;
    }

  private:
    
    const UGraph& _ugraph;
    PointMap _point_map;
    
  };

  namespace _planarity_bits {

    template <typename ColorMap>
    class KempeFilter {
    public:
      typedef typename ColorMap::Key Key;
      typedef bool Value;

      KempeFilter(const ColorMap& color_map, 
		  const typename ColorMap::Value& first,
		  const typename ColorMap::Value& second)
	: _color_map(color_map), _first(first), _second(second) {}

      Value operator[](const Key& key) const {
	return _color_map[key] == _first || _color_map[key] == _second;
      }

    private:
      const ColorMap& _color_map;
      typename ColorMap::Value _first, _second;
    };
  }

  /// \ingroup planar
  ///
  /// \brief Coloring planar graphs
  ///
  /// The graph coloring problem is the coloring of the graph nodes
  /// such way that there are not adjacent nodes with the same
  /// color. The planar graphs can be always colored with four colors,
  /// it is proved by Appel and Haken and their proofs provide a
  /// quadratic time algorithm for four coloring, but it could not be
  /// used to implement efficient algorithm. The five and six coloring
  /// can be made in linear time, but in this class the five coloring
  /// has quadratic worst case time complexity. The two coloring (if
  /// possible) is solvable with a graph search algorithm and it is
  /// implemented in \ref bipartitePartitions() function in Lemon. To
  /// decide whether the planar graph is three colorable is
  /// NP-complete.
  ///
  /// This class contains member functions for calculate colorings
  /// with five and six colors. The six coloring algorithm is a simple
  /// greedy coloring on the backward minimum outgoing order of nodes.
  /// This order can be computed such way, that in each phase the node
  /// with least outgoing edges to unprocessed nodes is choosen. This
  /// order guarantees that at the time of coloring of a node it has
  /// at most five already colored adjacents. The five coloring
  /// algorithm works in the same way, but if the greedy approach
  /// fails to color with five color, ie. the node has five already
  /// different colored neighbours, it swaps the colors in one
  /// connected two colored set with the Kempe recoloring method.
  template <typename UGraph>
  class PlanarColoring {
  public:

    UGRAPH_TYPEDEFS(typename UGraph);

    /// \brief The map type for store color indexes
    typedef typename UGraph::template NodeMap<int> IndexMap;
    /// \brief The map type for store colors
    typedef ComposeMap<Palette, IndexMap> ColorMap;

    /// \brief Constructor
    ///
    /// Constructor
    /// \pre The ugraph should be simple, ie. loop and parallel edge free. 
    PlanarColoring(const UGraph& ugraph)
      : _ugraph(ugraph), _color_map(ugraph), _palette(0) {
      _palette.add(Color(1,0,0));
      _palette.add(Color(0,1,0));
      _palette.add(Color(0,0,1));
      _palette.add(Color(1,1,0));
      _palette.add(Color(1,0,1));
      _palette.add(Color(0,1,1));
    }

    /// \brief Returns the \e NodeMap of color indexes
    ///
    /// Returns the \e NodeMap of color indexes. The values are in the
    /// range \c [0..4] or \c [0..5] according to the five coloring or 
    /// six coloring.
    IndexMap colorIndexMap() const {
      return _color_map;
    }

    /// \brief Returns the \e NodeMap of colors
    ///
    /// Returns the \e NodeMap of colors. The values are five or six
    /// distinct \ref lemon::Color "colors".
    ColorMap colorMap() const {
      return composeMap(_palette, _color_map);
    }

    /// \brief Returns the color index of the node
    ///
    /// Returns the color index of the node. The values are in the
    /// range \c [0..4] or \c [0..5] according to the five coloring or
    /// six coloring.
    int colorIndex(const Node& node) const {
      return _color_map[node];
    }

    /// \brief Returns the color of the node
    ///
    /// Returns the color of the node. The values are five or six
    /// distinct \ref lemon::Color "colors".
    Color color(const Node& node) const {
      return _palette[_color_map[node]];
    }
    

    /// \brief Calculates a coloring with at most six colors
    ///
    /// This function calculates a coloring with at most six colors. The time
    /// complexity of this variant is linear in the size of the graph.
    /// \return %True when the algorithm could color the graph with six color.
    /// If the algorithm fails, then the graph could not be planar.
    bool runSixColoring() {
      
      typename UGraph::template NodeMap<int> heap_index(_ugraph, -1);
      BucketHeap<typename UGraph::template NodeMap<int> > heap(heap_index);
      
      for (NodeIt n(_ugraph); n != INVALID; ++n) {
	_color_map[n] = -2;
	heap.push(n, countOutEdges(_ugraph, n));
      }
      
      std::vector<Node> order;
      
      while (!heap.empty()) {
	Node n = heap.top();
	heap.pop();
	_color_map[n] = -1;
	order.push_back(n);
	for (OutEdgeIt e(_ugraph, n); e != INVALID; ++e) {
	  Node t = _ugraph.runningNode(e); 
	  if (_color_map[t] == -2) {
	    heap.decrease(t, heap[t] - 1);
	  }
	}
      }

      for (int i = order.size() - 1; i >= 0; --i) {
	std::vector<bool> forbidden(6, false);
	for (OutEdgeIt e(_ugraph, order[i]); e != INVALID; ++e) {
	  Node t = _ugraph.runningNode(e); 
	  if (_color_map[t] != -1) {
	    forbidden[_color_map[t]] = true;
	  }
	}
       	for (int k = 0; k < 6; ++k) {
	  if (!forbidden[k]) {
	    _color_map[order[i]] = k;
	    break;
	  }
	}
	if (_color_map[order[i]] == -1) {
	  return false;
	}
      }
      return true;
    }

  private:

    bool recolor(const Node& u, const Node& v) {
      int ucolor = _color_map[u];
      int vcolor = _color_map[v];
      typedef _planarity_bits::KempeFilter<IndexMap> KempeFilter;
      KempeFilter filter(_color_map, ucolor, vcolor);

      typedef NodeSubUGraphAdaptor<const UGraph, const KempeFilter> KempeUGraph;
      KempeUGraph kempe_ugraph(_ugraph, filter);
      
      std::vector<Node> comp;
      Bfs<KempeUGraph> bfs(kempe_ugraph);
      bfs.init();
      bfs.addSource(u);
      while (!bfs.emptyQueue()) {
	Node n = bfs.nextNode();
	if (n == v) return false;
	comp.push_back(n);
	bfs.processNextNode();
      }

      int scolor = ucolor + vcolor;
      for (int i = 0; i < static_cast<int>(comp.size()); ++i) {
	_color_map[comp[i]] = scolor - _color_map[comp[i]]; 
      }

      return true;
    }

    template <typename EmbeddingMap>
    void kempeRecoloring(const Node& node, const EmbeddingMap& embedding) {
      std::vector<Node> nodes;
      nodes.reserve(4);

      for (Edge e = OutEdgeIt(_ugraph, node); e != INVALID; e = embedding[e]) {
	Node t = _ugraph.target(e); 
	if (_color_map[t] != -1) {
	  nodes.push_back(t);
	  if (nodes.size() == 4) break;
	}
      }
      
      int color = _color_map[nodes[0]];
      if (recolor(nodes[0], nodes[2])) {
	_color_map[node] = color;
      } else {
	color = _color_map[nodes[1]];
	recolor(nodes[1], nodes[3]);
	_color_map[node] = color;
      }
    }

  public:

    /// \brief Calculates a coloring with at most five colors
    ///
    /// This function calculates a coloring with at most five
    /// colors. The worst case time complexity of this variant is
    /// quadratic in the size of the graph.
    template <typename EmbeddingMap>
    void runFiveColoring(const EmbeddingMap& embedding) {
      
      typename UGraph::template NodeMap<int> heap_index(_ugraph, -1);
      BucketHeap<typename UGraph::template NodeMap<int> > heap(heap_index);
      
      for (NodeIt n(_ugraph); n != INVALID; ++n) {
	_color_map[n] = -2;
	heap.push(n, countOutEdges(_ugraph, n));
      }
      
      std::vector<Node> order;
      
      while (!heap.empty()) {
	Node n = heap.top();
	heap.pop();
	_color_map[n] = -1;
	order.push_back(n);
	for (OutEdgeIt e(_ugraph, n); e != INVALID; ++e) {
	  Node t = _ugraph.runningNode(e); 
	  if (_color_map[t] == -2) {
	    heap.decrease(t, heap[t] - 1);
	  }
	}
      }

      for (int i = order.size() - 1; i >= 0; --i) {
	std::vector<bool> forbidden(5, false);
	for (OutEdgeIt e(_ugraph, order[i]); e != INVALID; ++e) {
	  Node t = _ugraph.runningNode(e); 
	  if (_color_map[t] != -1) {
	    forbidden[_color_map[t]] = true;
	  }
	}
	for (int k = 0; k < 5; ++k) {
	  if (!forbidden[k]) {
	    _color_map[order[i]] = k;
	    break;
	  }
	}
	if (_color_map[order[i]] == -1) {
	  kempeRecoloring(order[i], embedding);
	}
      }
    }

    /// \brief Calculates a coloring with at most five colors
    ///
    /// This function calculates a coloring with at most five
    /// colors. The worst case time complexity of this variant is
    /// quadratic in the size of the graph, but it most cases it does
    /// not have to use Kempe recoloring method, in this case it is
    /// equivalent with the runSixColoring() algorithm.
    /// \return %True when the graph is planar.
    bool runFiveColoring() {
      PlanarEmbedding<UGraph> pe(_ugraph);
      if (!pe.run()) return false;
      
      runFiveColoring(pe.embeddingMap());
      return true;
    }

  private:
    
    const UGraph& _ugraph;
    IndexMap _color_map;
    Palette _palette;
  };

}

#endif
