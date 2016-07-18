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

#ifndef LEMON_BITS_EDGE_SET_EXTENDER_H
#define LEMON_BITS_EDGE_SET_EXTENDER_H

#include <lemon/bits/invalid.h>
#include <lemon/error.h>

#include <lemon/bits/default_map.h>

///\ingroup graphbits
///\file
///\brief Extenders for the edge set types
namespace lemon {

  /// \ingroup graphbits
  ///
  /// \brief Extender for the EdgeSets
  template <typename Base>
  class EdgeSetExtender : public Base {
  public:

    typedef Base Parent;
    typedef EdgeSetExtender Graph;

    // Base extensions

    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;

    int maxId(Node) const {
      return Parent::maxNodeId();
    }

    int maxId(Edge) const {
      return Parent::maxEdgeId();
    }

    Node fromId(int id, Node) const {
      return Parent::nodeFromId(id);
    }

    Edge fromId(int id, Edge) const {
      return Parent::edgeFromId(id);
    }

    Node oppositeNode(const Node &n, const Edge &e) const {
      if (n == Parent::source(e))
	return Parent::target(e);
      else if(n==Parent::target(e))
	return Parent::source(e);
      else
	return INVALID;
    }


    // Alteration notifier extensions

    /// The edge observer registry.
    typedef AlterationNotifier<EdgeSetExtender, Edge> EdgeNotifier;

  protected:

    mutable EdgeNotifier edge_notifier;

  public:

    using Parent::notifier;

    /// \brief Gives back the edge alteration notifier.
    ///
    /// Gives back the edge alteration notifier.
    EdgeNotifier& notifier(Edge) const {
      return edge_notifier;
    }

    // Iterable extensions

    class NodeIt : public Node { 
      const Graph* graph;
    public:

      NodeIt() {}

      NodeIt(Invalid i) : Node(i) { }

      explicit NodeIt(const Graph& _graph) : graph(&_graph) {
	_graph.first(static_cast<Node&>(*this));
      }

      NodeIt(const Graph& _graph, const Node& node) 
	: Node(node), graph(&_graph) {}

      NodeIt& operator++() { 
	graph->next(*this);
	return *this; 
      }

    };


    class EdgeIt : public Edge { 
      const Graph* graph;
    public:

      EdgeIt() { }

      EdgeIt(Invalid i) : Edge(i) { }

      explicit EdgeIt(const Graph& _graph) : graph(&_graph) {
	_graph.first(static_cast<Edge&>(*this));
      }

      EdgeIt(const Graph& _graph, const Edge& e) : 
	Edge(e), graph(&_graph) { }

      EdgeIt& operator++() { 
	graph->next(*this);
	return *this; 
      }

    };


    class OutEdgeIt : public Edge { 
      const Graph* graph;
    public:

      OutEdgeIt() { }

      OutEdgeIt(Invalid i) : Edge(i) { }

      OutEdgeIt(const Graph& _graph, const Node& node) 
	: graph(&_graph) {
	_graph.firstOut(*this, node);
      }

      OutEdgeIt(const Graph& _graph, const Edge& edge) 
	: Edge(edge), graph(&_graph) {}

      OutEdgeIt& operator++() { 
	graph->nextOut(*this);
	return *this; 
      }

    };


    class InEdgeIt : public Edge { 
      const Graph* graph;
    public:

      InEdgeIt() { }

      InEdgeIt(Invalid i) : Edge(i) { }

      InEdgeIt(const Graph& _graph, const Node& node) 
	: graph(&_graph) {
	_graph.firstIn(*this, node);
      }

      InEdgeIt(const Graph& _graph, const Edge& edge) : 
	Edge(edge), graph(&_graph) {}

      InEdgeIt& operator++() { 
	graph->nextIn(*this);
	return *this; 
      }

    };

    /// \brief Base node of the iterator
    ///
    /// Returns the base node (ie. the source in this case) of the iterator
    Node baseNode(const OutEdgeIt &e) const {
      return Parent::source(static_cast<const Edge&>(e));
    }
    /// \brief Running node of the iterator
    ///
    /// Returns the running node (ie. the target in this case) of the
    /// iterator
    Node runningNode(const OutEdgeIt &e) const {
      return Parent::target(static_cast<const Edge&>(e));
    }

    /// \brief Base node of the iterator
    ///
    /// Returns the base node (ie. the target in this case) of the iterator
    Node baseNode(const InEdgeIt &e) const {
      return Parent::target(static_cast<const Edge&>(e));
    }
    /// \brief Running node of the iterator
    ///
    /// Returns the running node (ie. the source in this case) of the
    /// iterator
    Node runningNode(const InEdgeIt &e) const {
      return Parent::source(static_cast<const Edge&>(e));
    }

    using Parent::first;

    // Mappable extension
    
    template <typename _Value>
    class EdgeMap 
      : public MapExtender<DefaultMap<Graph, Edge, _Value> > {
    public:
      typedef EdgeSetExtender Graph;
      typedef MapExtender<DefaultMap<Graph, Edge, _Value> > Parent;

      explicit EdgeMap(const Graph& _g) 
	: Parent(_g) {}
      EdgeMap(const Graph& _g, const _Value& _v) 
	: Parent(_g, _v) {}

      EdgeMap& operator=(const EdgeMap& cmap) {
	return operator=<EdgeMap>(cmap);
      }

      template <typename CMap>
      EdgeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
	return *this;
      }

    };


    // Alteration extension

    Edge addEdge(const Node& from, const Node& to) {
      Edge edge = Parent::addEdge(from, to);
      notifier(Edge()).add(edge);
      return edge;
    }
    
    void clear() {
      notifier(Edge()).clear();
      Parent::clear();
    }

    void erase(const Edge& edge) {
      notifier(Edge()).erase(edge);
      Parent::erase(edge);
    }

    EdgeSetExtender() {
      edge_notifier.setContainer(*this);
    }

    ~EdgeSetExtender() {
      edge_notifier.clear();
    }

  };


  /// \ingroup graphbits
  ///
  /// \brief Extender for the UEdgeSets
  template <typename Base>
  class UEdgeSetExtender : public Base {

  public:

    typedef Base Parent;
    typedef UEdgeSetExtender Graph;

    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;
    typedef typename Parent::UEdge UEdge;


    int maxId(Node) const {
      return Parent::maxNodeId();
    }

    int maxId(Edge) const {
      return Parent::maxEdgeId();
    }

    int maxId(UEdge) const {
      return Parent::maxUEdgeId();
    }

    Node fromId(int id, Node) const {
      return Parent::nodeFromId(id);
    }

    Edge fromId(int id, Edge) const {
      return Parent::edgeFromId(id);
    }

    UEdge fromId(int id, UEdge) const {
      return Parent::uEdgeFromId(id);
    }

    Node oppositeNode(const Node &n, const UEdge &e) const {
      if( n == Parent::source(e))
	return Parent::target(e);
      else if( n == Parent::target(e))
	return Parent::source(e);
      else
	return INVALID;
    }

    Edge oppositeEdge(const Edge &e) const {
      return Parent::direct(e, !Parent::direction(e));
    }

    using Parent::direct;
    Edge direct(const UEdge &ue, const Node &s) const {
      return Parent::direct(ue, Parent::source(ue) == s);
    }

    typedef AlterationNotifier<UEdgeSetExtender, Edge> EdgeNotifier;
    typedef AlterationNotifier<UEdgeSetExtender, UEdge> UEdgeNotifier;


  protected:

    mutable EdgeNotifier edge_notifier;
    mutable UEdgeNotifier uedge_notifier;

  public:

    using Parent::notifier;
    
    EdgeNotifier& notifier(Edge) const {
      return edge_notifier;
    }

    UEdgeNotifier& notifier(UEdge) const {
      return uedge_notifier;
    }


    class NodeIt : public Node { 
      const Graph* graph;
    public:

      NodeIt() {}

      NodeIt(Invalid i) : Node(i) { }

      explicit NodeIt(const Graph& _graph) : graph(&_graph) {
	_graph.first(static_cast<Node&>(*this));
      }

      NodeIt(const Graph& _graph, const Node& node) 
	: Node(node), graph(&_graph) {}

      NodeIt& operator++() { 
	graph->next(*this);
	return *this; 
      }

    };


    class EdgeIt : public Edge { 
      const Graph* graph;
    public:

      EdgeIt() { }

      EdgeIt(Invalid i) : Edge(i) { }

      explicit EdgeIt(const Graph& _graph) : graph(&_graph) {
	_graph.first(static_cast<Edge&>(*this));
      }

      EdgeIt(const Graph& _graph, const Edge& e) : 
	Edge(e), graph(&_graph) { }

      EdgeIt& operator++() { 
	graph->next(*this);
	return *this; 
      }

    };


    class OutEdgeIt : public Edge { 
      const Graph* graph;
    public:

      OutEdgeIt() { }

      OutEdgeIt(Invalid i) : Edge(i) { }

      OutEdgeIt(const Graph& _graph, const Node& node) 
	: graph(&_graph) {
	_graph.firstOut(*this, node);
      }

      OutEdgeIt(const Graph& _graph, const Edge& edge) 
	: Edge(edge), graph(&_graph) {}

      OutEdgeIt& operator++() { 
	graph->nextOut(*this);
	return *this; 
      }

    };


    class InEdgeIt : public Edge { 
      const Graph* graph;
    public:

      InEdgeIt() { }

      InEdgeIt(Invalid i) : Edge(i) { }

      InEdgeIt(const Graph& _graph, const Node& node) 
	: graph(&_graph) {
	_graph.firstIn(*this, node);
      }

      InEdgeIt(const Graph& _graph, const Edge& edge) : 
	Edge(edge), graph(&_graph) {}

      InEdgeIt& operator++() { 
	graph->nextIn(*this);
	return *this; 
      }

    };


    class UEdgeIt : public Parent::UEdge { 
      const Graph* graph;
    public:

      UEdgeIt() { }

      UEdgeIt(Invalid i) : UEdge(i) { }

      explicit UEdgeIt(const Graph& _graph) : graph(&_graph) {
	_graph.first(static_cast<UEdge&>(*this));
      }

      UEdgeIt(const Graph& _graph, const UEdge& e) : 
	UEdge(e), graph(&_graph) { }

      UEdgeIt& operator++() { 
	graph->next(*this);
	return *this; 
      }

    };

    class IncEdgeIt : public Parent::UEdge {
      friend class UEdgeSetExtender;
      const Graph* graph;
      bool direction;
    public:

      IncEdgeIt() { }

      IncEdgeIt(Invalid i) : UEdge(i), direction(false) { }

      IncEdgeIt(const Graph& _graph, const Node &n) : graph(&_graph) {
	_graph.firstInc(*this, direction, n);
      }

      IncEdgeIt(const Graph& _graph, const UEdge &ue, const Node &n)
	: graph(&_graph), UEdge(ue) {
	direction = (_graph.source(ue) == n);
      }

      IncEdgeIt& operator++() {
	graph->nextInc(*this, direction);
	return *this; 
      }
    };

    /// \brief Base node of the iterator
    ///
    /// Returns the base node (ie. the source in this case) of the iterator
    Node baseNode(const OutEdgeIt &e) const {
      return Parent::source(static_cast<const Edge&>(e));
    }
    /// \brief Running node of the iterator
    ///
    /// Returns the running node (ie. the target in this case) of the
    /// iterator
    Node runningNode(const OutEdgeIt &e) const {
      return Parent::target(static_cast<const Edge&>(e));
    }

    /// \brief Base node of the iterator
    ///
    /// Returns the base node (ie. the target in this case) of the iterator
    Node baseNode(const InEdgeIt &e) const {
      return Parent::target(static_cast<const Edge&>(e));
    }
    /// \brief Running node of the iterator
    ///
    /// Returns the running node (ie. the source in this case) of the
    /// iterator
    Node runningNode(const InEdgeIt &e) const {
      return Parent::source(static_cast<const Edge&>(e));
    }

    /// Base node of the iterator
    ///
    /// Returns the base node of the iterator
    Node baseNode(const IncEdgeIt &e) const {
      return e.direction ? source(e) : target(e);
    }
    /// Running node of the iterator
    ///
    /// Returns the running node of the iterator
    Node runningNode(const IncEdgeIt &e) const {
      return e.direction ? target(e) : source(e);
    }


    template <typename _Value>
    class EdgeMap 
      : public MapExtender<DefaultMap<Graph, Edge, _Value> > {
    public:
      typedef UEdgeSetExtender Graph;
      typedef MapExtender<DefaultMap<Graph, Edge, _Value> > Parent;

      EdgeMap(const Graph& _g) 
	: Parent(_g) {}
      EdgeMap(const Graph& _g, const _Value& _v) 
	: Parent(_g, _v) {}

      EdgeMap& operator=(const EdgeMap& cmap) {
	return operator=<EdgeMap>(cmap);
      }

      template <typename CMap>
      EdgeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
	return *this;
      }

    };


    template <typename _Value>
    class UEdgeMap 
      : public MapExtender<DefaultMap<Graph, UEdge, _Value> > {
    public:
      typedef UEdgeSetExtender Graph;
      typedef MapExtender<DefaultMap<Graph, UEdge, _Value> > Parent;

      UEdgeMap(const Graph& _g) 
	: Parent(_g) {}

      UEdgeMap(const Graph& _g, const _Value& _v) 
	: Parent(_g, _v) {}

      UEdgeMap& operator=(const UEdgeMap& cmap) {
	return operator=<UEdgeMap>(cmap);
      }

      template <typename CMap>
      UEdgeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
	return *this;
      }

    };


    // Alteration extension

    UEdge addEdge(const Node& from, const Node& to) {
      UEdge uedge = Parent::addEdge(from, to);
      notifier(UEdge()).add(uedge);
      std::vector<Edge> edges;
      edges.push_back(Parent::direct(uedge, true));
      edges.push_back(Parent::direct(uedge, false));
      notifier(Edge()).add(edges);
      return uedge;
    }
    
    void clear() {
      notifier(Edge()).clear();
      notifier(UEdge()).clear();
      Parent::clear();
    }

    void erase(const UEdge& uedge) {
      std::vector<Edge> edges;
      edges.push_back(Parent::direct(uedge, true));
      edges.push_back(Parent::direct(uedge, false));
      notifier(Edge()).erase(edges);
      notifier(UEdge()).erase(uedge);
      Parent::erase(uedge);
    }


    UEdgeSetExtender() {
      edge_notifier.setContainer(*this);
      uedge_notifier.setContainer(*this);
    }

    ~UEdgeSetExtender() {
      uedge_notifier.clear();
      edge_notifier.clear();
    }
    
  };

}

#endif
