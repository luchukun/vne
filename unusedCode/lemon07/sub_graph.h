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

#ifndef LEMON_SUB_GRAPH_H
#define LEMON_SUB_GRAPH_H

#include <lemon/graph_adaptor.h>
#include <lemon/bits/graph_adaptor_extender.h>
#include <lemon/bits/default_map.h>

/// \ingroup semi_adaptors
/// \file
/// \brief Subgraphs.
///
/// Graphs with filtered edge and node set.

namespace lemon {

  /// \brief Base for the SubGraph.
  ///
  /// Base for the SubGraph.
  template <typename _Graph>
  class SubGraphBase : public GraphAdaptorBase<const _Graph> {
  public:
    typedef _Graph Graph;
    typedef SubGraphBase<_Graph> SubGraph;
    typedef GraphAdaptorBase<const _Graph> Parent;
    typedef Parent Base;

    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;


  protected:

    class NodesImpl;
    class EdgesImpl;

    SubGraphBase() {}

    void construct(const Graph& _graph, NodesImpl& _nodes, EdgesImpl& _edges) {
      Parent::setGraph(_graph);
      nodes = &_nodes;
      edges = &_edges;
      firstNode = INVALID;

      Node node;
      Parent::first(node);
      while (node != INVALID) {
	(*nodes)[node].prev = node;
	(*nodes)[node].firstIn = INVALID;
	(*nodes)[node].firstOut = INVALID;
	Parent::next(node);
      }

      Edge edge;
      Parent::first(edge);
      while (edge != INVALID) {
	(*edges)[edge].prevOut = edge;
	Parent::next(edge);
      }
    }

  public:

    void first(Node& node) const {
      node = firstNode;
    }
    void next(Node& node) const {
      node = (*nodes)[node].next;
    }

    void first(Edge& edge) const {
      Node node = firstNode;
      while (node != INVALID && (*nodes)[node].firstOut == INVALID) {
	node = (*nodes)[node].next;
      }
      if (node == INVALID) {
	edge = INVALID;
      } else {
	edge = (*nodes)[node].firstOut;
      }
    }
    void next(Edge& edge) const {
      if ((*edges)[edge].nextOut != INVALID) {
	edge = (*edges)[edge].nextOut;
      } else {
	Node node = (*nodes)[source(edge)].next;
	while (node != INVALID && (*nodes)[node].firstOut == INVALID) {
	  node = (*nodes)[node].next;
	}
	if (node == INVALID) {
	  edge = INVALID;
	} else {
	  edge = (*nodes)[node].firstOut;
	}
      }
    }

    void firstOut(Edge& edge, const Node& node) const {
      edge = (*nodes)[node].firstOut;
    }
    void nextOut(Edge& edge) const {
      edge = (*edges)[edge].nextOut;
    }

    void firstIn(Edge& edge, const Node& node) const {
      edge = (*nodes)[node].firstIn;
    }
    void nextIn(Edge& edge) const {
      edge = (*edges)[edge].nextIn;
    }

    /// \brief Returns true when the given node is hidden.
    ///
    /// Returns true when the given node is hidden.
    bool hidden(const Node& node) const {
      return (*nodes)[node].prev == node;
    }

    /// \brief Hide the given node in the sub-graph.
    ///
    /// Hide the given node in the sub graph. It just lace out from
    /// the linked lists the given node. If there are incoming or outgoing
    /// edges into or from this node then all of these will be hidden.
    void hide(const Node& node) {
      if (hidden(node)) return;
      Edge edge;
      firstOut(edge, node);
      while (edge != INVALID) {
	hide(edge);
	firstOut(edge, node);
      }

      firstOut(edge, node);
      while (edge != INVALID) {
	hide(edge);
	firstOut(edge, node);
      }
      if ((*nodes)[node].prev != INVALID) {
	(*nodes)[(*nodes)[node].prev].next = (*nodes)[node].next;
      } else {
	firstNode = (*nodes)[node].next;
      }
      if ((*nodes)[node].next != INVALID) {
	(*nodes)[(*nodes)[node].next].prev = (*nodes)[node].prev;
      }
      (*nodes)[node].prev = node;
      (*nodes)[node].firstIn = INVALID;
      (*nodes)[node].firstOut = INVALID;
    }

    /// \brief Unhide the given node in the sub-graph.
    ///
    /// Unhide the given node in the sub graph. It just lace in the given
    /// node into the linked lists.
    void unHide(const Node& node) {
      if (!hidden(node)) return;
      (*nodes)[node].next = firstNode;
      (*nodes)[node].prev = INVALID;
      if ((*nodes)[node].next != INVALID) {
	(*nodes)[(*nodes)[node].next].prev = node;
      }
      firstNode = node;
    }

    /// \brief Returns true when the given edge is hidden.
    ///
    /// Returns true when the given edge is hidden.
    bool hidden(const Edge& edge) const {
      return (*edges)[edge].prevOut == edge;
    }

    /// \brief Hide the given edge in the sub-graph.
    ///
    /// Hide the given edge in the sub graph. It just lace out from
    /// the linked lists the given edge.
    void hide(const Edge& edge) {
      if (hidden(edge)) return;
      if ((*edges)[edge].prevOut == edge) return;
      if ((*edges)[edge].prevOut != INVALID) {
	(*edges)[(*edges)[edge].prevOut].nextOut = (*edges)[edge].nextOut;
      } else {
	(*nodes)[source(edge)].firstOut = (*edges)[edge].nextOut;
      }
      if ((*edges)[edge].nextOut != INVALID) {
	(*edges)[(*edges)[edge].nextOut].prevOut = (*edges)[edge].prevOut;
      }

      if ((*edges)[edge].prevIn != INVALID) {
	(*edges)[(*edges)[edge].prevIn].nextIn = (*edges)[edge].nextIn;
      } else {
	(*nodes)[target(edge)].firstIn = (*edges)[edge].nextIn;
      }
      if ((*edges)[edge].nextIn != INVALID) {
	(*edges)[(*edges)[edge].nextIn].prevIn = (*edges)[edge].prevIn;
      }
      (*edges)[edge].next = edge;
    }

    /// \brief Unhide the given edge in the sub-graph.
    ///
    /// Unhide the given edge in the sub graph. It just lace in the given
    /// edge into the linked lists. If the source or the target of the
    /// node is hidden then it will unhide it.
    void unHide(const Edge& edge) {
      if (!hidden(edge)) return;

      Node node;

      node = Parent::source(edge);
      unHide(node);
      (*edges)[edge].nextOut = (*nodes)[node].firstOut;
      (*edges)[edge].prevOut = INVALID;
      if ((*edges)[edge].nextOut != INVALID) {
	(*edges)[(*edges)[edge].nextOut].prevOut = edge;
      }
      (*nodes)[node].firstOut = edge;

      node = Parent::target(edge);
      unHide(node);
      (*edges)[edge].nextIn = (*nodes)[node].firstIn;
      (*edges)[edge].prevIn = INVALID;
      if ((*edges)[edge].nextIn != INVALID) {
	(*edges)[(*edges)[edge].nextIn].prevIn = edge;
      }
      (*nodes)[node].firstIn = edge;      
    }
    
    typedef False NodeNumTag;
    typedef False EdgeNumTag;

  protected:
    struct NodeT {
      Node prev, next;
      Edge firstIn, firstOut;
    };
    class NodesImpl : public DefaultMap<Graph, Node, NodeT> {
      friend class SubGraphBase;
    public:
      typedef DefaultMap<Graph, Node, NodeT> Parent;

      NodesImpl(SubGraph& _adaptor, const Graph& _graph) 
	: Parent(_graph), adaptor(_adaptor) {}

      virtual ~NodesImpl() {}

      virtual void build() {
	Parent::build();
	Node node;
	adaptor.Base::first(node);
	while (node != INVALID) {
	  Parent::operator[](node).prev = node;
	  Parent::operator[](node).firstIn = INVALID;
	  Parent::operator[](node).firstOut = INVALID;
	  adaptor.Base::next(node);
	}
      }

      virtual void clear() {
	adaptor.firstNode = INVALID;
	Parent::clear();
      }

      virtual void add(const Node& node) {
	Parent::add(node);
	Parent::operator[](node).prev = node;
	Parent::operator[](node).firstIn = INVALID;
	Parent::operator[](node).firstOut = INVALID;
      }

      virtual void add(const std::vector<Node>& nodes) {
	Parent::add(nodes);
	for (int i = 0; i < (int)nodes.size(); ++i) {
	  Parent::operator[](nodes[i]).prev = nodes[i];
	  Parent::operator[](nodes[i]).firstIn = INVALID;
	  Parent::operator[](nodes[i]).firstOut = INVALID;
	}
      } 

      virtual void erase(const Node& node) {
	adaptor.hide(node);
	Parent::erase(node);
      }

      virtual void erase(const std::vector<Node>& nodes) {
	for (int i = 0; i < (int)nodes.size(); ++i) {
	  adaptor.hide(nodes[i]);
	}
	Parent::erase(nodes);
      }

    private:
      SubGraph& adaptor;
    };

    struct EdgeT {
      Edge prevOut, nextOut;
      Edge prevIn, nextIn;
    };
    class EdgesImpl : public DefaultMap<Graph, Edge, EdgeT> {
      friend class SubGraphBase;
    public:
      typedef DefaultMap<Graph, Edge, EdgeT> Parent;

      EdgesImpl(SubGraph& _adaptor, const Graph& _graph) 
	: Parent(_graph), adaptor(_adaptor) {}

      virtual ~EdgesImpl() {}

      virtual void build() {
	Parent::build();
	Edge edge;
	adaptor.Base::first(edge);
	while (edge != INVALID) {
	  Parent::operator[](edge).prevOut = edge;
	  adaptor.Base::next(edge);
	}
      }

      virtual void clear() {
	Node node;
	adaptor.first(node);
	while (node != INVALID) {
	  (*adaptor.nodes).firstIn = INVALID;
	  (*adaptor.nodes).firstOut = INVALID;
	  adaptor.next(node);
	}
	Parent::clear();
      }

      virtual void add(const Edge& edge) {
	Parent::add(edge);
	Parent::operator[](edge).prevOut = edge;
      }

      virtual void add(const std::vector<Edge>& edges) {
	Parent::add(edges);
	for (int i = 0; i < (int)edges.size(); ++i) {
	  Parent::operator[](edges[i]).prevOut = edges[i];
	}
      }

      virtual void erase(const Edge& edge) {
	adaptor.hide(edge);
	Parent::erase(edge);
      }

      virtual void erase(const std::vector<Edge>& edges) {
	for (int i = 0; i < (int)edges.size(); ++i) {
	  adaptor.hide(edges[i]);
	}
	Parent::erase(edges);
      }

    private:
      SubGraph& adaptor;
    };

    NodesImpl* nodes;
    EdgesImpl* edges;
    Node firstNode;
  };

  /// \ingroup semi_adaptors
  ///
  /// \brief Graph which uses a subset of another graph's nodes and edges.
  ///
  /// Graph which uses a subset of another graph's nodes and edges. This class
  /// is an alternative to the SubGraphAdaptor which is created for the
  /// same reason. The main difference between the two class that it
  /// makes linked lists on the unhidden nodes and edges what cause that
  /// on sparse subgraphs the algorithms can be more efficient and some times
  /// provide better time complexity. On other way this implemetation is
  /// less efficient in most case when the subgraph filters out only
  /// a few nodes or edges.
  /// 
  /// \see SubGraphAdaptor
  /// \see EdgeSubGraphBase
  template <typename Graph>
  class SubGraph 
    : public GraphAdaptorExtender< SubGraphBase<Graph> > {
  public:
    typedef GraphAdaptorExtender< SubGraphBase<Graph> > Parent;
  public:
    /// \brief Constructor for sub-graph.
    ///
    /// Constructor for sub-graph. Initially all the edges and nodes
    /// are hidden in the graph.
    SubGraph(const Graph& _graph) 
      : Parent(), nodes(*this, _graph), edges(*this, _graph) { 
      Parent::construct(_graph, nodes, edges);
    }
  private:
    typename Parent::NodesImpl nodes;
    typename Parent::EdgesImpl edges;
  };

  /// \brief Base for the EdgeSubGraph.
  ///
  /// Base for the EdgeSubGraph.
  template <typename _Graph>
  class EdgeSubGraphBase : public GraphAdaptorBase<const _Graph> {
  public:
    typedef _Graph Graph;
    typedef EdgeSubGraphBase<_Graph> SubGraph;
    typedef GraphAdaptorBase<const _Graph> Parent;
    typedef Parent Base;

    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;


  protected:

    class NodesImpl;
    class EdgesImpl;

    EdgeSubGraphBase() {}

    void construct(const Graph& _graph, NodesImpl& _nodes, EdgesImpl& _edges) {
      Parent::setGraph(_graph);
      nodes = &_nodes;
      edges = &_edges;

      Node node;
      Parent::first(node);
      while (node != INVALID) {
	(*nodes)[node].firstIn = INVALID;
	(*nodes)[node].firstOut = INVALID;
	Parent::next(node);
      }

      Edge edge;
      Parent::first(edge);
      while (edge != INVALID) {
	(*edges)[edge].prevOut = edge;
	Parent::next(edge);
      }
    }

  public:

    void first(Node& node) const {
      Parent::first(node);
    }
    void next(Node& node) const {
      Parent::next(node);
    }

    void first(Edge& edge) const {
      Node node;
      Parent::first(node);
      while (node != INVALID && (*nodes)[node].firstOut == INVALID) {
	Parent::next(node);
      }
      if (node == INVALID) {
	edge = INVALID;
      } else {
	edge = (*nodes)[node].firstOut;
      }
    }
    void next(Edge& edge) const {
      if ((*edges)[edge].nextOut != INVALID) {
	edge = (*edges)[edge].nextOut;
      } else {
	Node node = source(edge);
	Parent::next(node);
	while (node != INVALID && (*nodes)[node].firstOut == INVALID) {
	  Parent::next(node);
	}
	if (node == INVALID) {
	  edge = INVALID;
	} else {
	  edge = (*nodes)[node].firstOut;
	}
      }
    }

    void firstOut(Edge& edge, const Node& node) const {
      edge = (*nodes)[node].firstOut;
    }
    void nextOut(Edge& edge) const {
      edge = (*edges)[edge].nextOut;
    }

    void firstIn(Edge& edge, const Node& node) const {
      edge = (*nodes)[node].firstIn;
    }
    void nextIn(Edge& edge) const {
      edge = (*edges)[edge].nextIn;
    }

    /// \brief Returns true when the given edge is hidden.
    ///
    /// Returns true when the given edge is hidden.
    bool hidden(const Edge& edge) const {
      return (*edges)[edge].prevOut == edge;
    }

    /// \brief Hide the given edge in the sub-graph.
    ///
    /// Hide the given edge in the sub graph. It just lace out from
    /// the linked lists the given edge.
    void hide(const Edge& edge) {
      if (hidden(edge)) return;
      if ((*edges)[edge].prevOut != INVALID) {
	(*edges)[(*edges)[edge].prevOut].nextOut = (*edges)[edge].nextOut;
      } else {
	(*nodes)[source(edge)].firstOut = (*edges)[edge].nextOut;
      }
      if ((*edges)[edge].nextOut != INVALID) {
	(*edges)[(*edges)[edge].nextOut].prevOut = (*edges)[edge].prevOut;
      }

      if ((*edges)[edge].prevIn != INVALID) {
	(*edges)[(*edges)[edge].prevIn].nextIn = (*edges)[edge].nextIn;
      } else {
	(*nodes)[target(edge)].firstIn = (*edges)[edge].nextIn;
      }
      if ((*edges)[edge].nextIn != INVALID) {
	(*edges)[(*edges)[edge].nextIn].prevIn = (*edges)[edge].prevIn;
      }
      (*edges)[edge].prevOut = edge;
    }

    /// \brief Unhide the given edge in the sub-graph.
    ///
    /// Unhide the given edge in the sub graph. It just lace in the given
    /// edge into the linked lists.
    void unHide(const Edge& edge) {
      if (!hidden(edge)) return;
      Node node;

      node = Parent::source(edge);
      (*edges)[edge].nextOut = (*nodes)[node].firstOut;
      (*edges)[edge].prevOut = INVALID;
      if ((*edges)[edge].nextOut != INVALID) {
	(*edges)[(*edges)[edge].nextOut].prevOut = edge;
      }
      (*nodes)[node].firstOut = edge;

      node = Parent::target(edge);
      (*edges)[edge].nextIn = (*nodes)[node].firstIn;
      (*edges)[edge].prevIn = INVALID;
      if ((*edges)[edge].nextIn != INVALID) {
	(*edges)[(*edges)[edge].nextIn].prevIn = edge;
      }
      (*nodes)[node].firstIn = edge;      
    }
    
  protected:
    struct NodeT {
      Edge firstIn, firstOut;
    };
    class NodesImpl : public DefaultMap<Graph, Node, NodeT> {
      friend class EdgeSubGraphBase;
    public:
      typedef DefaultMap<Graph, Node, NodeT> Parent;

      NodesImpl(SubGraph& _adaptor, const Graph& _graph) 
	: Parent(_graph), adaptor(_adaptor) {}

      virtual ~NodesImpl() {}

      virtual void build() {
	Parent::build();
	Node node;
	adaptor.Base::first(node);
	while (node != INVALID) {
	  Parent::operator[](node).firstIn = INVALID;
	  Parent::operator[](node).firstOut = INVALID;
	  adaptor.Base::next(node);
	}
      }

      virtual void add(const Node& node) {
	Parent::add(node);
	Parent::operator[](node).firstIn = INVALID;
	Parent::operator[](node).firstOut = INVALID;
      }

      virtual void add(const std::vector<Node>& nodes) {
        Parent::add(nodes);
        for (int i = 0; i < (int)nodes.size(); ++i) {
          Parent::operator[](nodes[i]).firstIn = INVALID;
          Parent::operator[](nodes[i]).firstOut = INVALID;
        }
      }

    private:
      SubGraph& adaptor;
    };

    struct EdgeT {
      Edge prevOut, nextOut;
      Edge prevIn, nextIn;
    };
    class EdgesImpl : public DefaultMap<Graph, Edge, EdgeT> {
      friend class EdgeSubGraphBase;
    public:
      typedef DefaultMap<Graph, Edge, EdgeT> Parent;

      EdgesImpl(SubGraph& _adaptor, const Graph& _graph) 
	: Parent(_graph), adaptor(_adaptor) {}

      virtual ~EdgesImpl() {}

      virtual void build() {
	Parent::build();
	Edge edge;
	adaptor.Base::first(edge);
	while (edge != INVALID) {
	  Parent::operator[](edge).prevOut = edge;
	  adaptor.Base::next(edge);
	}
      }

      virtual void clear() {
	Node node;
	adaptor.Base::first(node);
	while (node != INVALID) {
	  (*adaptor.nodes)[node].firstIn = INVALID;
	  (*adaptor.nodes)[node].firstOut = INVALID;
	  adaptor.Base::next(node);
	}
	Parent::clear();
      }

      virtual void add(const Edge& edge) {
	Parent::add(edge);
	Parent::operator[](edge).prevOut = edge;
      }

      virtual void add(const std::vector<Edge>& edges) {
	Parent::add(edges);
	for (int i = 0; i < (int)edges.size(); ++i) {
	  Parent::operator[](edges[i]).prevOut = edges[i];
	}
      }

      virtual void erase(const Edge& edge) {
	adaptor.hide(edge);
	Parent::erase(edge);
      }

      virtual void erase(const std::vector<Edge>& edges) {
	for (int i = 0; i < (int)edges.size(); ++i) {
	  adaptor.hide(edges[i]);
	}
	Parent::erase(edges);
      }

    private:
      SubGraph& adaptor;
    };

    NodesImpl* nodes;
    EdgesImpl* edges;
  };

  /// \ingroup semi_adaptors
  ///
  /// \brief Graph which uses a subset of another graph's edges.
  ///
  /// Graph which uses a subset of another graph's edges. This class
  /// is an alternative to the EdgeSubGraphAdaptor which is created for the
  /// same reason. The main difference between the two class that it
  /// makes linked lists on the unhidden edges what cause that
  /// on sparse subgraphs the algorithms can be more efficient and some times
  /// provide better time complexity. On other way this implemetation is
  /// less efficient in most case when the subgraph filters out only
  /// a few edges.
  /// 
  /// \see EdgeSubGraphAdaptor
  /// \see EdgeSubGraphBase
  template <typename Graph>
  class EdgeSubGraph 
    : public GraphAdaptorExtender< EdgeSubGraphBase<Graph> > {
  public:
    typedef GraphAdaptorExtender< EdgeSubGraphBase<Graph> > Parent;
  public:
    /// \brief Constructor for sub-graph.
    ///
    /// Constructor for sub-graph. Initially all the edges are hidden in the 
    /// graph.
    EdgeSubGraph(const Graph& _graph) 
      : Parent(), nodes(*this, _graph), edges(*this, _graph) { 
      Parent::construct(_graph, nodes, edges);
    }
  private:
    typename Parent::NodesImpl nodes;
    typename Parent::EdgesImpl edges;
  };


//   template<typename Graph, typename Number, 
// 	   typename CapacityMap, typename FlowMap>
//   class ResGraph
//     : public IterableGraphExtender<EdgeSubGraphBase<
//     UGraphAdaptor<Graph> > > {
//   public:
//     typedef IterableGraphExtender<EdgeSubGraphBase<
//       UGraphAdaptor<Graph> > > Parent;

//   protected:
//     UGraphAdaptor<Graph> u;

//     const CapacityMap* capacity;
//     FlowMap* flow;

//     typename Parent::NodesImpl nodes;
//     typename Parent::EdgesImpl edges;

//     void setCapacityMap(const CapacityMap& _capacity) {
//       capacity=&_capacity;
//     }

//     void setFlowMap(FlowMap& _flow) {
//       flow=&_flow;
//     }

//   public:

//     typedef typename UGraphAdaptor<Graph>::Node Node;
//     typedef typename UGraphAdaptor<Graph>::Edge Edge;
//     typedef typename UGraphAdaptor<Graph>::UEdge UEdge;

//     ResGraphAdaptor(Graph& _graph, 
// 		    const CapacityMap& _capacity, FlowMap& _flow) 
//       : Parent(), u(_graph), capacity(&_capacity), flow(&_flow),
// 	nodes(*this, _graph), edges(*this, _graph) {
//       Parent::construct(u, nodes, edges);
//       setFlowMap(_flow);
//       setCapacityMap(_capacity);
//       typename Graph::Edge edge; 
//       for (_graph.first(edge); edge != INVALID; _graph.next(edge)) {
// 	if ((*flow)[edge] != (*capacity)[edge]) {
// 	  Parent::unHide(direct(edge, true));
// 	}
// 	if ((*flow)[edge] != 0) {
// 	  Parent::unHide(direct(edge, false));
// 	}
//       }
//     }

//     void augment(const Edge& e, Number a) {
//       if (direction(e)) {
// 	flow->set(e, (*flow)[e]+a);
//       } else { 
// 	flow->set(e, (*flow)[e]-a);
//       }
//       if ((*flow)[e] == (*capacity)[e]) {
// 	Parent::hide(e);
//       } else {
// 	Parent::unHide(e);
//       }
//       if ((*flow)[e] == 0) {
// 	Parent::hide(oppositeEdge(e));
//       } else {
// 	Parent::unHide(oppositeEdge(e));
//       }
//     }

//     Number resCap(const Edge& e) {
//       if (direction(e)) { 
// 	return (*capacity)[e]-(*flow)[e]; 
//       } else { 
// 	return (*flow)[e];
//       }
//     }
    
//     bool direction(const Edge& edge) const {
//       return Parent::getGraph().direction(edge);
//     }

//     Edge direct(const UEdge& edge, bool direction) const {
//       return Parent::getGraph().direct(edge, direction);
//     }

//     Edge direct(const UEdge& edge, const Node& node) const {
//       return Parent::getGraph().direct(edge, node);
//     }

//     Edge oppositeEdge(const Edge& edge) const {
//       return Parent::getGraph().oppositeEdge(edge);
//     }

//     /// \brief Residual capacity map.
//     ///
//     /// In generic residual graphs the residual capacity can be obtained 
//     /// as a map. 
//     class ResCap {
//     protected:
//       const ResGraphAdaptor* res_graph;
//     public:
//       typedef Number Value;
//       typedef Edge Key;
//       ResCap(const ResGraphAdaptor& _res_graph) 
// 	: res_graph(&_res_graph) {}
//       Number operator[](const Edge& e) const {
// 	return res_graph->resCap(e);
//       }
//     };
//   };

}

#endif
