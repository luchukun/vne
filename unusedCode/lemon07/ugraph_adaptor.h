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

#ifndef LEMON_UGRAPH_ADAPTOR_H
#define LEMON_UGRAPH_ADAPTOR_H

///\ingroup graph_adaptors
///\file
///\brief Several graph adaptors.
///
///This file contains several useful ugraph adaptor functions.
///
///\author Balazs Dezso

#include <lemon/bits/invalid.h>
#include <lemon/maps.h>

#include <lemon/bits/graph_adaptor_extender.h>

#include <lemon/bits/traits.h>

#include <iostream>

namespace lemon {

  /// \brief Base type for the Graph Adaptors
  ///
  /// This is the base type for most of LEMON graph adaptors. 
  /// This class implements a trivial graph adaptor i.e. it only wraps the 
  /// functions and types of the graph. The purpose of this class is to 
  /// make easier implementing graph adaptors. E.g. if an adaptor is 
  /// considered which differs from the wrapped graph only in some of its 
  /// functions or types, then it can be derived from GraphAdaptor, and only 
  /// the differences should be implemented.
  ///
  /// \author Balazs Dezso 
  template<typename _UGraph>
  class UGraphAdaptorBase {
  public:
    typedef _UGraph Graph;
    typedef Graph ParentGraph;

  protected:
    Graph* graph;

    UGraphAdaptorBase() : graph(0) {}

    void setGraph(Graph& _graph) { graph=&_graph; }

  public:
    UGraphAdaptorBase(Graph& _graph) : graph(&_graph) {}
 
    typedef typename Graph::Node Node;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::UEdge UEdge;
   
    void first(Node& i) const { graph->first(i); }
    void first(Edge& i) const { graph->first(i); }
    void first(UEdge& i) const { graph->first(i); }
    void firstIn(Edge& i, const Node& n) const { graph->firstIn(i, n); }
    void firstOut(Edge& i, const Node& n ) const { graph->firstOut(i, n); }
    void firstInc(UEdge &i, bool &d, const Node &n) const {
      graph->firstInc(i, d, n);
    }

    void next(Node& i) const { graph->next(i); }
    void next(Edge& i) const { graph->next(i); }
    void next(UEdge& i) const { graph->next(i); }
    void nextIn(Edge& i) const { graph->nextIn(i); }
    void nextOut(Edge& i) const { graph->nextOut(i); }
    void nextInc(UEdge &i, bool &d) const { graph->nextInc(i, d); }

    Node source(const UEdge& e) const { return graph->source(e); }
    Node target(const UEdge& e) const { return graph->target(e); }

    Node source(const Edge& e) const { return graph->source(e); }
    Node target(const Edge& e) const { return graph->target(e); }

    typedef NodeNumTagIndicator<Graph> NodeNumTag;
    int nodeNum() const { return graph->nodeNum(); }
    
    typedef EdgeNumTagIndicator<Graph> EdgeNumTag;
    int edgeNum() const { return graph->edgeNum(); }
    int uEdgeNum() const { return graph->uEdgeNum(); }

    typedef FindEdgeTagIndicator<Graph> FindEdgeTag;
    Edge findEdge(const Node& u, const Node& v, 
		  const Edge& prev = INVALID) {
      return graph->findEdge(u, v, prev);
    }
    UEdge findUEdge(const Node& u, const Node& v, 
                    const UEdge& prev = INVALID) {
      return graph->findUEdge(u, v, prev);
    }
  
    Node addNode() const { return graph->addNode(); }
    UEdge addEdge(const Node& u, const Node& v) const { 
      return graph->addEdge(u, v); 
    }

    void erase(const Node& i) const { graph->erase(i); }
    void erase(const UEdge& i) const { graph->erase(i); }
  
    void clear() const { graph->clear(); }
    
    bool direction(const Edge& e) const { return graph->direction(e); }
    Edge direct(const UEdge& e, bool d) const { return graph->direct(e, d); }

    int id(const Node& v) const { return graph->id(v); }
    int id(const Edge& e) const { return graph->id(e); }
    int id(const UEdge& e) const { return graph->id(e); }

    Node fromNodeId(int ix) const {
      return graph->fromNodeId(ix);
    }

    Edge fromEdgeId(int ix) const {
      return graph->fromEdgeId(ix);
    }

    UEdge fromUEdgeId(int ix) const {
      return graph->fromUEdgeId(ix);
    }

    int maxNodeId() const {
      return graph->maxNodeId();
    }

    int maxEdgeId() const {
      return graph->maxEdgeId();
    }

    int maxUEdgeId() const {
      return graph->maxEdgeId();
    }

    typedef typename ItemSetTraits<Graph, Node>::ItemNotifier NodeNotifier;

    NodeNotifier& notifier(Node) const {
      return graph->notifier(Node());
    } 

    typedef typename ItemSetTraits<Graph, Edge>::ItemNotifier EdgeNotifier;

    EdgeNotifier& notifier(Edge) const {
      return graph->notifier(Edge());
    } 

    typedef typename ItemSetTraits<Graph, UEdge>::ItemNotifier UEdgeNotifier;

    UEdgeNotifier& notifier(UEdge) const {
      return graph->notifier(UEdge());
    } 

    template <typename _Value>
    class NodeMap : public Graph::template NodeMap<_Value> {
    public:
      typedef typename Graph::template NodeMap<_Value> Parent;
      explicit NodeMap(const UGraphAdaptorBase<Graph>& ga) 
	: Parent(*ga.graph) {}
      NodeMap(const UGraphAdaptorBase<Graph>& ga, const _Value& value)
	: Parent(*ga.graph, value) {}

      NodeMap& operator=(const NodeMap& cmap) {
	return operator=<NodeMap>(cmap);
      }

      template <typename CMap>
      NodeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
        return *this;
      }

    };

    template <typename _Value>
    class EdgeMap : public Graph::template EdgeMap<_Value> {
    public:
      typedef typename Graph::template EdgeMap<_Value> Parent;
      explicit EdgeMap(const UGraphAdaptorBase<Graph>& ga) 
	: Parent(*ga.graph) {}
      EdgeMap(const UGraphAdaptorBase<Graph>& ga, const _Value& value)
	: Parent(*ga.graph, value) {}

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
    class UEdgeMap : public Graph::template UEdgeMap<_Value> {
    public:
      typedef typename Graph::template UEdgeMap<_Value> Parent;
      explicit UEdgeMap(const UGraphAdaptorBase<Graph>& ga) 
	: Parent(*ga.graph) {}
      UEdgeMap(const UGraphAdaptorBase<Graph>& ga, const _Value& value)
	: Parent(*ga.graph, value) {}

      UEdgeMap& operator=(const UEdgeMap& cmap) {
	return operator=<UEdgeMap>(cmap);
      }

      template <typename CMap>
      UEdgeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
        return *this;
      }
    };

  };

  /// \ingroup graph_adaptors
  ///
  /// \brief Trivial undirected graph adaptor
  ///
  /// This class is an adaptor which does not change the adapted undirected
  /// graph. It can be used only to test the undirected graph adaptors.
  template <typename _UGraph>
  class UGraphAdaptor 
    : public UGraphAdaptorExtender< UGraphAdaptorBase<_UGraph> > { 
  public:
    typedef _UGraph Graph;
    typedef UGraphAdaptorExtender<UGraphAdaptorBase<_UGraph> > Parent;
  protected:
    UGraphAdaptor() : Parent() {}

  public:
    explicit UGraphAdaptor(Graph& _graph) { setGraph(_graph); }
  };

  template <typename _UGraph, typename NodeFilterMap, 
	    typename UEdgeFilterMap, bool checked = true>
  class SubUGraphAdaptorBase : public UGraphAdaptorBase<_UGraph> {
  public:
    typedef _UGraph Graph;
    typedef SubUGraphAdaptorBase Adaptor;
    typedef UGraphAdaptorBase<_UGraph> Parent;
  protected:

    NodeFilterMap* node_filter_map;
    UEdgeFilterMap* uedge_filter_map;

    SubUGraphAdaptorBase() 
      : Parent(), node_filter_map(0), uedge_filter_map(0) { }

    void setNodeFilterMap(NodeFilterMap& _node_filter_map) {
      node_filter_map=&_node_filter_map;
    }
    void setUEdgeFilterMap(UEdgeFilterMap& _uedge_filter_map) {
      uedge_filter_map=&_uedge_filter_map;
    }

  public:

    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;
    typedef typename Parent::UEdge UEdge;

    void first(Node& i) const { 
      Parent::first(i); 
      while (i!=INVALID && !(*node_filter_map)[i]) Parent::next(i); 
    }

    void first(Edge& i) const { 
      Parent::first(i); 
      while (i!=INVALID && (!(*uedge_filter_map)[i] 
	     || !(*node_filter_map)[Parent::source(i)]
	     || !(*node_filter_map)[Parent::target(i)])) Parent::next(i); 
    }

    void first(UEdge& i) const { 
      Parent::first(i); 
      while (i!=INVALID && (!(*uedge_filter_map)[i] 
	     || !(*node_filter_map)[Parent::source(i)]
	     || !(*node_filter_map)[Parent::target(i)])) Parent::next(i); 
    }

    void firstIn(Edge& i, const Node& n) const { 
      Parent::firstIn(i, n); 
      while (i!=INVALID && (!(*uedge_filter_map)[i] 
	     || !(*node_filter_map)[Parent::source(i)])) Parent::nextIn(i); 
    }

    void firstOut(Edge& i, const Node& n) const { 
      Parent::firstOut(i, n); 
      while (i!=INVALID && (!(*uedge_filter_map)[i] 
	     || !(*node_filter_map)[Parent::target(i)])) Parent::nextOut(i); 
    }

    void firstInc(UEdge& i, bool& d, const Node& n) const { 
      Parent::firstInc(i, d, n); 
      while (i!=INVALID && (!(*uedge_filter_map)[i] 
            || !(*node_filter_map)[Parent::source(i)]
            || !(*node_filter_map)[Parent::target(i)])) Parent::nextInc(i, d); 
    }

    void next(Node& i) const { 
      Parent::next(i); 
      while (i!=INVALID && !(*node_filter_map)[i]) Parent::next(i); 
    }

    void next(Edge& i) const { 
      Parent::next(i); 
      while (i!=INVALID && (!(*uedge_filter_map)[i] 
	     || !(*node_filter_map)[Parent::source(i)]
	     || !(*node_filter_map)[Parent::target(i)])) Parent::next(i); 
    }

    void next(UEdge& i) const { 
      Parent::next(i); 
      while (i!=INVALID && (!(*uedge_filter_map)[i] 
	     || !(*node_filter_map)[Parent::source(i)]
	     || !(*node_filter_map)[Parent::target(i)])) Parent::next(i); 
    }

    void nextIn(Edge& i) const { 
      Parent::nextIn(i); 
      while (i!=INVALID && (!(*uedge_filter_map)[i] 
	     || !(*node_filter_map)[Parent::source(i)])) Parent::nextIn(i); 
    }

    void nextOut(Edge& i) const { 
      Parent::nextOut(i); 
      while (i!=INVALID && (!(*uedge_filter_map)[i] 
	     || !(*node_filter_map)[Parent::target(i)])) Parent::nextOut(i); 
    }

    void nextInc(UEdge& i, bool& d) const { 
      Parent::nextInc(i, d); 
      while (i!=INVALID && (!(*uedge_filter_map)[i]
            || !(*node_filter_map)[Parent::source(i)]
            || !(*node_filter_map)[Parent::target(i)])) Parent::nextInc(i, d); 
    }

    /// \brief Hide the given node in the graph.
    ///
    /// This function hides \c n in the graph, i.e. the iteration 
    /// jumps over it. This is done by simply setting the value of \c n  
    /// to be false in the corresponding node-map.
    void hide(const Node& n) const { node_filter_map->set(n, false); }

    /// \brief Hide the given undirected edge in the graph.
    ///
    /// This function hides \c e in the graph, i.e. the iteration 
    /// jumps over it. This is done by simply setting the value of \c e  
    /// to be false in the corresponding uedge-map.
    void hide(const UEdge& e) const { uedge_filter_map->set(e, false); }

    /// \brief Unhide the given node in the graph.
    ///
    /// The value of \c n is set to be true in the node-map which stores 
    /// hide information. If \c n was hidden previuosly, then it is shown 
    /// again
     void unHide(const Node& n) const { node_filter_map->set(n, true); }

    /// \brief Hide the given undirected edge in the graph.
    ///
    /// The value of \c e is set to be true in the uedge-map which stores 
    /// hide information. If \c e was hidden previuosly, then it is shown 
    /// again
    void unHide(const UEdge& e) const { uedge_filter_map->set(e, true); }

    /// \brief Returns true if \c n is hidden.
    ///
    /// Returns true if \c n is hidden.
    bool hidden(const Node& n) const { return !(*node_filter_map)[n]; }

    /// \brief Returns true if \c e is hidden.
    ///
    /// Returns true if \c e is hidden.
    bool hidden(const UEdge& e) const { return !(*uedge_filter_map)[e]; }

    typedef False NodeNumTag;
    typedef False EdgeNumTag;

    typedef FindEdgeTagIndicator<Graph> FindEdgeTag;
    Edge findEdge(const Node& u, const Node& v, 
		  const Edge& prev = INVALID) {
      if (!(*node_filter_map)[u] || !(*node_filter_map)[v]) {
        return INVALID;
      }
      Edge edge = Parent::findEdge(u, v, prev);
      while (edge != INVALID && !(*uedge_filter_map)[edge]) {
        edge = Parent::findEdge(u, v, edge);
      }
      return edge;
    }
    UEdge findUEdge(const Node& u, const Node& v, 
		  const UEdge& prev = INVALID) {
      if (!(*node_filter_map)[u] || !(*node_filter_map)[v]) {
        return INVALID;
      }
      UEdge uedge = Parent::findUEdge(u, v, prev);
      while (uedge != INVALID && !(*uedge_filter_map)[uedge]) {
        uedge = Parent::findUEdge(u, v, uedge);
      }
      return uedge;
    }

    template <typename _Value>
    class NodeMap 
      : public SubMapExtender<Adaptor, 
                              typename Parent::template NodeMap<_Value> > 
    {
    public:
      typedef Adaptor Graph;
      typedef SubMapExtender<Adaptor, typename Parent::
                             template NodeMap<_Value> > Parent;
    
      NodeMap(const Graph& g) 
	: Parent(g) {}
      NodeMap(const Graph& g, const _Value& v) 
	: Parent(g, v) {}
    
      NodeMap& operator=(const NodeMap& cmap) {
	return operator=<NodeMap>(cmap);
      }
    
      template <typename CMap>
      NodeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
	return *this;
      }
    };

    template <typename _Value>
    class EdgeMap 
      : public SubMapExtender<Adaptor, 
                              typename Parent::template EdgeMap<_Value> > 
    {
    public:
      typedef Adaptor Graph;
      typedef SubMapExtender<Adaptor, typename Parent::
                             template EdgeMap<_Value> > Parent;
    
      EdgeMap(const Graph& g) 
	: Parent(g) {}
      EdgeMap(const Graph& g, const _Value& v) 
	: Parent(g, v) {}
    
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
      : public SubMapExtender<Adaptor, 
                              typename Parent::template UEdgeMap<_Value> > 
    {
    public:
      typedef Adaptor Graph;
      typedef SubMapExtender<Adaptor, typename Parent::
                             template UEdgeMap<_Value> > Parent;
    
      UEdgeMap(const Graph& g) 
	: Parent(g) {}
      UEdgeMap(const Graph& g, const _Value& v) 
	: Parent(g, v) {}
    
      UEdgeMap& operator=(const UEdgeMap& cmap) {
	return operator=<UEdgeMap>(cmap);
      }
    
      template <typename CMap>
      UEdgeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
	return *this;
      }
    };

  };

  template <typename _UGraph, typename NodeFilterMap, typename UEdgeFilterMap>
  class SubUGraphAdaptorBase<_UGraph, NodeFilterMap, UEdgeFilterMap, false> 
    : public UGraphAdaptorBase<_UGraph> {
  public:
    typedef _UGraph Graph;
    typedef SubUGraphAdaptorBase Adaptor;
    typedef UGraphAdaptorBase<_UGraph> Parent;
  protected:
    NodeFilterMap* node_filter_map;
    UEdgeFilterMap* uedge_filter_map;
    SubUGraphAdaptorBase() : Parent(), 
			    node_filter_map(0), uedge_filter_map(0) { }

    void setNodeFilterMap(NodeFilterMap& _node_filter_map) {
      node_filter_map=&_node_filter_map;
    }
    void setUEdgeFilterMap(UEdgeFilterMap& _uedge_filter_map) {
      uedge_filter_map=&_uedge_filter_map;
    }

  public:

    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;
    typedef typename Parent::UEdge UEdge;

    void first(Node& i) const { 
      Parent::first(i); 
      while (i!=INVALID && !(*node_filter_map)[i]) Parent::next(i); 
    }

    void first(Edge& i) const { 
      Parent::first(i); 
      while (i!=INVALID && !(*uedge_filter_map)[i]) Parent::next(i); 
    }

    void first(UEdge& i) const { 
      Parent::first(i); 
      while (i!=INVALID && !(*uedge_filter_map)[i]) Parent::next(i); 
    }

    void firstIn(Edge& i, const Node& n) const { 
      Parent::firstIn(i, n); 
      while (i!=INVALID && !(*uedge_filter_map)[i]) Parent::nextIn(i); 
    }

    void firstOut(Edge& i, const Node& n) const { 
      Parent::firstOut(i, n); 
      while (i!=INVALID && !(*uedge_filter_map)[i]) Parent::nextOut(i); 
    }

    void firstInc(UEdge& i, bool& d, const Node& n) const { 
      Parent::firstInc(i, d, n); 
      while (i!=INVALID && !(*uedge_filter_map)[i]) Parent::nextInc(i, d); 
    }

    void next(Node& i) const { 
      Parent::next(i); 
      while (i!=INVALID && !(*node_filter_map)[i]) Parent::next(i); 
    }
    void next(Edge& i) const { 
      Parent::next(i); 
      while (i!=INVALID && !(*uedge_filter_map)[i]) Parent::next(i); 
    }
    void next(UEdge& i) const { 
      Parent::next(i); 
      while (i!=INVALID && !(*uedge_filter_map)[i]) Parent::next(i); 
    }
    void nextIn(Edge& i) const { 
      Parent::nextIn(i); 
      while (i!=INVALID && !(*uedge_filter_map)[i]) Parent::nextIn(i); 
    }

    void nextOut(Edge& i) const { 
      Parent::nextOut(i); 
      while (i!=INVALID && !(*uedge_filter_map)[i]) Parent::nextOut(i); 
    }
    void nextInc(UEdge& i, bool& d) const { 
      Parent::nextInc(i, d); 
      while (i!=INVALID && !(*uedge_filter_map)[i]) Parent::nextInc(i, d); 
    }

    /// \brief Hide the given node in the graph.
    ///
    /// This function hides \c n in the graph, i.e. the iteration 
    /// jumps over it. This is done by simply setting the value of \c n  
    /// to be false in the corresponding node-map.
    void hide(const Node& n) const { node_filter_map->set(n, false); }

    /// \brief Hide the given undirected edge in the graph.
    ///
    /// This function hides \c e in the graph, i.e. the iteration 
    /// jumps over it. This is done by simply setting the value of \c e  
    /// to be false in the corresponding uedge-map.
    void hide(const UEdge& e) const { uedge_filter_map->set(e, false); }

    /// \brief Unhide the given node in the graph.
    ///
    /// The value of \c n is set to be true in the node-map which stores 
    /// hide information. If \c n was hidden previuosly, then it is shown 
    /// again
     void unHide(const Node& n) const { node_filter_map->set(n, true); }

    /// \brief Hide the given undirected edge in the graph.
    ///
    /// The value of \c e is set to be true in the uedge-map which stores 
    /// hide information. If \c e was hidden previuosly, then it is shown 
    /// again
    void unHide(const UEdge& e) const { uedge_filter_map->set(e, true); }

    /// \brief Returns true if \c n is hidden.
    ///
    /// Returns true if \c n is hidden.
    bool hidden(const Node& n) const { return !(*node_filter_map)[n]; }

    /// \brief Returns true if \c e is hidden.
    ///
    /// Returns true if \c e is hidden.
    bool hidden(const UEdge& e) const { return !(*uedge_filter_map)[e]; }

    typedef False NodeNumTag;
    typedef False EdgeNumTag;

    typedef FindEdgeTagIndicator<Graph> FindEdgeTag;
    Edge findEdge(const Node& u, const Node& v, 
		  const Edge& prev = INVALID) {
      Edge edge = Parent::findEdge(u, v, prev);
      while (edge != INVALID && !(*uedge_filter_map)[edge]) {
        edge = Parent::findEdge(u, v, edge);
      }
      return edge;
    }
    UEdge findUEdge(const Node& u, const Node& v, 
		  const UEdge& prev = INVALID) {
      UEdge uedge = Parent::findUEdge(u, v, prev);
      while (uedge != INVALID && !(*uedge_filter_map)[uedge]) {
        uedge = Parent::findUEdge(u, v, uedge);
      }
      return uedge;
    }

    template <typename _Value>
    class NodeMap 
      : public SubMapExtender<Adaptor, 
                              typename Parent::template NodeMap<_Value> > 
    {
    public:
      typedef Adaptor Graph;
      typedef SubMapExtender<Adaptor, typename Parent::
                             template NodeMap<_Value> > Parent;
    
      NodeMap(const Graph& g) 
	: Parent(g) {}
      NodeMap(const Graph& g, const _Value& v) 
	: Parent(g, v) {}
    
      NodeMap& operator=(const NodeMap& cmap) {
	return operator=<NodeMap>(cmap);
      }
    
      template <typename CMap>
      NodeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
	return *this;
      }
    };

    template <typename _Value>
    class EdgeMap 
      : public SubMapExtender<Adaptor, 
                              typename Parent::template EdgeMap<_Value> > 
    {
    public:
      typedef Adaptor Graph;
      typedef SubMapExtender<Adaptor, typename Parent::
                             template EdgeMap<_Value> > Parent;
    
      EdgeMap(const Graph& g) 
	: Parent(g) {}
      EdgeMap(const Graph& g, const _Value& v) 
	: Parent(g, v) {}
    
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
      : public SubMapExtender<Adaptor, 
                              typename Parent::template UEdgeMap<_Value> > 
    {
    public:
      typedef Adaptor Graph;
      typedef SubMapExtender<Adaptor, typename Parent::
                             template UEdgeMap<_Value> > Parent;
    
      UEdgeMap(const Graph& g) 
	: Parent(g) {}
      UEdgeMap(const Graph& g, const _Value& v) 
	: Parent(g, v) {}
    
      UEdgeMap& operator=(const UEdgeMap& cmap) {
	return operator=<UEdgeMap>(cmap);
      }
    
      template <typename CMap>
      UEdgeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
	return *this;
      }
    };
  };

  /// \ingroup graph_adaptors
  ///
  /// \brief A graph adaptor for hiding nodes and edges from an undirected 
  /// graph.
  /// 
  /// SubUGraphAdaptor shows the undirected graph with filtered node-set and 
  /// edge-set. If the \c checked parameter is true then it filters the edgeset
  /// to do not get invalid edges without source or target.
  /// 
  /// If the \c checked template parameter is false then we have to note that 
  /// the node-iterator cares only the filter on the node-set, and the 
  /// edge-iterator cares only the filter on the edge-set.
  /// This way the edge-map
  /// should filter all edges which's source or target is filtered by the 
  /// node-filter.
  ///
  template<typename _UGraph, typename NodeFilterMap, 
	   typename UEdgeFilterMap, bool checked = true>
  class SubUGraphAdaptor : 
    public UGraphAdaptorExtender<
    SubUGraphAdaptorBase<_UGraph, NodeFilterMap, UEdgeFilterMap, checked> > {
  public:
    typedef _UGraph Graph;
    typedef UGraphAdaptorExtender<
      SubUGraphAdaptorBase<_UGraph, NodeFilterMap, UEdgeFilterMap> > Parent;
  protected:
    SubUGraphAdaptor() { }
  public:
    SubUGraphAdaptor(Graph& _graph, NodeFilterMap& _node_filter_map, 
		    UEdgeFilterMap& _uedge_filter_map) { 
      setGraph(_graph);
      setNodeFilterMap(_node_filter_map);
      setUEdgeFilterMap(_uedge_filter_map);
    }
  };

  template<typename UGraph, typename NodeFilterMap, typename EdgeFilterMap>
  SubUGraphAdaptor<const UGraph, NodeFilterMap, EdgeFilterMap>
  subUGraphAdaptor(const UGraph& graph, 
                   NodeFilterMap& nfm, EdgeFilterMap& efm) {
    return SubUGraphAdaptor<const UGraph, NodeFilterMap, EdgeFilterMap>
      (graph, nfm, efm);
  }

  template<typename UGraph, typename NodeFilterMap, typename EdgeFilterMap>
  SubUGraphAdaptor<const UGraph, const NodeFilterMap, EdgeFilterMap>
  subUGraphAdaptor(const UGraph& graph, 
                   NodeFilterMap& nfm, EdgeFilterMap& efm) {
    return SubUGraphAdaptor<const UGraph, const NodeFilterMap, EdgeFilterMap>
      (graph, nfm, efm);
  }

  template<typename UGraph, typename NodeFilterMap, typename EdgeFilterMap>
  SubUGraphAdaptor<const UGraph, NodeFilterMap, const EdgeFilterMap>
  subUGraphAdaptor(const UGraph& graph, 
                   NodeFilterMap& nfm, EdgeFilterMap& efm) {
    return SubUGraphAdaptor<const UGraph, NodeFilterMap, const EdgeFilterMap>
      (graph, nfm, efm);
  }

  template<typename UGraph, typename NodeFilterMap, typename EdgeFilterMap>
  SubUGraphAdaptor<const UGraph, const NodeFilterMap, const EdgeFilterMap>
  subUGraphAdaptor(const UGraph& graph, 
                   NodeFilterMap& nfm, EdgeFilterMap& efm) {
    return SubUGraphAdaptor<const UGraph, const NodeFilterMap, 
      const EdgeFilterMap>(graph, nfm, efm);
  }

  /// \ingroup graph_adaptors
  ///
  /// \brief An adaptor for hiding nodes from an undirected graph.
  ///
  /// An adaptor for hiding nodes from an undirected graph.  This
  /// adaptor specializes SubUGraphAdaptor in the way that only the
  /// node-set can be filtered. In usual case the checked parameter is
  /// true, we get the induced subgraph. But if the checked parameter
  /// is false then we can filter only isolated nodes.
  template<typename _UGraph, typename NodeFilterMap, bool checked = true>
  class NodeSubUGraphAdaptor : 
    public SubUGraphAdaptor<_UGraph, NodeFilterMap, 
                            ConstMap<typename _UGraph::UEdge, bool>, checked> {
  public:
    typedef SubUGraphAdaptor<_UGraph, NodeFilterMap, 
                             ConstMap<typename _UGraph::UEdge, bool> > Parent;
    typedef _UGraph Graph;
  protected:
    ConstMap<typename _UGraph::UEdge, bool> const_true_map;

    NodeSubUGraphAdaptor() : const_true_map(true) {
      Parent::setUEdgeFilterMap(const_true_map);
    }

  public:
    NodeSubUGraphAdaptor(Graph& _graph, NodeFilterMap& _node_filter_map) : 
      Parent(), const_true_map(true) { 
      Parent::setGraph(_graph);
      Parent::setNodeFilterMap(_node_filter_map);
      Parent::setUEdgeFilterMap(const_true_map);
    }
  };

  template<typename UGraph, typename NodeFilterMap>
  NodeSubUGraphAdaptor<const UGraph, NodeFilterMap>
  nodeSubUGraphAdaptor(const UGraph& graph, NodeFilterMap& nfm) {
    return NodeSubUGraphAdaptor<const UGraph, NodeFilterMap>(graph, nfm);
  }

  template<typename UGraph, typename NodeFilterMap>
  NodeSubUGraphAdaptor<const UGraph, const NodeFilterMap>
  nodeSubUGraphAdaptor(const UGraph& graph, const NodeFilterMap& nfm) {
    return NodeSubUGraphAdaptor<const UGraph, const NodeFilterMap>(graph, nfm);
  }

  /// \ingroup graph_adaptors
  ///
  /// \brief An adaptor for hiding undirected edges from an undirected graph.
  ///
  /// \warning Graph adaptors are in even more experimental state
  /// than the other parts of the lib. Use them at you own risk.
  ///
  /// An adaptor for hiding undirected edges from an undirected graph.
  /// This adaptor specializes SubUGraphAdaptor in the way that
  /// only the edge-set 
  /// can be filtered.
  template<typename _UGraph, typename UEdgeFilterMap>
  class EdgeSubUGraphAdaptor : 
    public SubUGraphAdaptor<_UGraph, ConstMap<typename _UGraph::Node,bool>, 
                            UEdgeFilterMap, false> {
  public:
    typedef SubUGraphAdaptor<_UGraph, ConstMap<typename _UGraph::Node,bool>, 
                             UEdgeFilterMap, false> Parent;
    typedef _UGraph Graph;
  protected:
    ConstMap<typename Graph::Node, bool> const_true_map;

    EdgeSubUGraphAdaptor() : const_true_map(true) {
      Parent::setNodeFilterMap(const_true_map);
    }

  public:

    EdgeSubUGraphAdaptor(Graph& _graph, UEdgeFilterMap& _uedge_filter_map) : 
      Parent(), const_true_map(true) { 
      Parent::setGraph(_graph);
      Parent::setNodeFilterMap(const_true_map);
      Parent::setUEdgeFilterMap(_uedge_filter_map);
    }

  };

  template<typename UGraph, typename EdgeFilterMap>
  EdgeSubUGraphAdaptor<const UGraph, EdgeFilterMap>
  edgeSubUGraphAdaptor(const UGraph& graph, EdgeFilterMap& efm) {
    return EdgeSubUGraphAdaptor<const UGraph, EdgeFilterMap>(graph, efm);
  }

  template<typename UGraph, typename EdgeFilterMap>
  EdgeSubUGraphAdaptor<const UGraph, const EdgeFilterMap>
  edgeSubUGraphAdaptor(const UGraph& graph, const EdgeFilterMap& efm) {
    return EdgeSubUGraphAdaptor<const UGraph, const EdgeFilterMap>(graph, efm);
  }

  /// \brief Base of direct undirected graph adaptor
  ///
  /// Base class of the direct undirected graph adaptor. All public member
  /// of this class can be used with the DirUGraphAdaptor too.
  /// \sa DirUGraphAdaptor
  template <typename _UGraph, typename _DirectionMap>
  class DirUGraphAdaptorBase {
  public:
    
    typedef _UGraph Graph;
    typedef _DirectionMap DirectionMap;

    typedef typename _UGraph::Node Node;
    typedef typename _UGraph::UEdge Edge;
   
    /// \brief Reverse edge
    /// 
    /// It reverse the given edge. It simply negate the direction in the map.
    void reverseEdge(const Edge& edge) {
      direction->set(edge, !(*direction)[edge]);
    }

    void first(Node& i) const { graph->first(i); }
    void first(Edge& i) const { graph->first(i); }
    void firstIn(Edge& i, const Node& n) const {
      bool d;
      graph->firstInc(i, d, n);
      while (i != INVALID && d == (*direction)[i]) graph->nextInc(i, d);
    }
    void firstOut(Edge& i, const Node& n ) const { 
      bool d;
      graph->firstInc(i, d, n);
      while (i != INVALID && d != (*direction)[i]) graph->nextInc(i, d);
    }

    void next(Node& i) const { graph->next(i); }
    void next(Edge& i) const { graph->next(i); }
    void nextIn(Edge& i) const {
      bool d = !(*direction)[i];
      graph->nextInc(i, d);
      while (i != INVALID && d == (*direction)[i]) graph->nextInc(i, d);
    }
    void nextOut(Edge& i) const {
      bool d = (*direction)[i];
      graph->nextInc(i, d);
      while (i != INVALID && d != (*direction)[i]) graph->nextInc(i, d);
    }

    Node source(const Edge& e) const { 
      return (*direction)[e] ? graph->source(e) : graph->target(e); 
    }
    Node target(const Edge& e) const { 
      return (*direction)[e] ? graph->target(e) : graph->source(e); 
    }

    typedef NodeNumTagIndicator<Graph> NodeNumTag;
    int nodeNum() const { return graph->nodeNum(); }
    
    typedef EdgeNumTagIndicator<Graph> EdgeNumTag;
    int edgeNum() const { return graph->uEdgeNum(); }

    typedef FindEdgeTagIndicator<Graph> FindEdgeTag;
    Edge findEdge(const Node& u, const Node& v, 
		  const Edge& prev = INVALID) {
      Edge edge = prev;
      bool d = edge == INVALID ? true : (*direction)[edge];
      if (d) {
        edge = graph->findUEdge(u, v, edge);
        while (edge != INVALID && !(*direction)[edge]) {
          graph->findUEdge(u, v, edge);
        }
        if (edge != INVALID) return edge;
      }
      graph->findUEdge(v, u, edge);
      while (edge != INVALID && (*direction)[edge]) {
        graph->findUEdge(u, v, edge);
      }
      return edge;
    }
  
    Node addNode() const { 
      return Node(graph->addNode()); 
    }

    Edge addEdge(const Node& u, const Node& v) const {
      Edge edge = graph->addEdge(u, v);
      direction->set(edge, graph->source(edge) == u);
      return edge; 
    }

    void erase(const Node& i) const { graph->erase(i); }
    void erase(const Edge& i) const { graph->erase(i); }
  
    void clear() const { graph->clear(); }
    
    int id(const Node& v) const { return graph->id(v); }
    int id(const Edge& e) const { return graph->id(e); }

    int maxNodeId() const {
      return graph->maxNodeId();
    }

    int maxEdgeId() const {
      return graph->maxEdgeId();
    }

    typedef typename ItemSetTraits<Graph, Node>::ItemNotifier NodeNotifier;

    NodeNotifier& notifier(Node) const {
      return graph->notifier(Node());
    } 

    typedef typename ItemSetTraits<Graph, Edge>::ItemNotifier EdgeNotifier;

    EdgeNotifier& notifier(Edge) const {
      return graph->notifier(Edge());
    } 

    template <typename _Value>
    class NodeMap : public _UGraph::template NodeMap<_Value> {
    public:

      typedef typename _UGraph::template NodeMap<_Value> Parent;

      explicit NodeMap(const DirUGraphAdaptorBase& ga) 
	: Parent(*ga.graph) {}

      NodeMap(const DirUGraphAdaptorBase& ga, const _Value& value)
	: Parent(*ga.graph, value) {}

      NodeMap& operator=(const NodeMap& cmap) {
        return operator=<NodeMap>(cmap);
      }

      template <typename CMap>
      NodeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
        return *this;
      }

    };

    template <typename _Value>
    class EdgeMap : public _UGraph::template UEdgeMap<_Value> {
    public:

      typedef typename _UGraph::template UEdgeMap<_Value> Parent;

      explicit EdgeMap(const DirUGraphAdaptorBase& ga) 
	: Parent(*ga.graph) { }

      EdgeMap(const DirUGraphAdaptorBase& ga, const _Value& value)
	: Parent(*ga.graph, value) { }

      EdgeMap& operator=(const EdgeMap& cmap) {
        return operator=<EdgeMap>(cmap);
      }

      template <typename CMap>
      EdgeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
        return *this;
      }
    };

    

  protected:
    Graph* graph;
    DirectionMap* direction;

    void setDirectionMap(DirectionMap& _direction) {
      direction = &_direction;
    }

    void setGraph(Graph& _graph) {
      graph = &_graph;
    }

  };


  /// \ingroup graph_adaptors
  ///
  /// \brief A directed graph is made from an undirected graph by an adaptor
  ///
  /// This adaptor gives a direction for each uedge in the undirected
  /// graph. The direction of the edges stored in the
  /// DirectionMap. This map is a bool map on the undirected edges. If
  /// the uedge is mapped to true then the direction of the directed
  /// edge will be the same as the default direction of the uedge. The
  /// edges can be easily reverted by the \ref
  /// DirUGraphAdaptorBase::reverseEdge "reverseEdge()" member in the
  /// adaptor.
  ///
  /// It can be used to solve orientation problems on directed graphs.
  /// By example how can we orient an undirected graph to get the minimum
  /// number of strongly connected components. If we orient the edges with
  /// the dfs algorithm out from the source then we will get such an 
  /// orientation. 
  ///
  /// We use the \ref DfsVisitor "visitor" interface of the 
  /// \ref DfsVisit "dfs" algorithm:
  ///\code
  /// template <typename DirMap>
  /// class OrientVisitor : public DfsVisitor<UGraph> {
  /// public:
  ///
  ///   OrientVisitor(const UGraph& ugraph, DirMap& dirMap)
  ///     : _ugraph(ugraph), _dirMap(dirMap), _processed(ugraph, false) {}
  ///
  ///   void discover(const Edge& edge) {
  ///     _processed.set(edge, true);
  ///     _dirMap.set(edge, _ugraph.direction(edge));
  ///   }
  ///
  ///   void examine(const Edge& edge) {
  ///     if (_processed[edge]) return;
  ///     _processed.set(edge, true);
  ///     _dirMap.set(edge, _ugraph.direction(edge));
  ///   }  
  /// 
  /// private:
  ///   const UGraph& _ugraph;  
  ///   DirMap& _dirMap;
  ///   UGraph::UEdgeMap<bool> _processed;
  /// };
  ///\endcode
  ///
  /// And now we can use the orientation:
  ///\code
  /// UGraph::UEdgeMap<bool> dmap(ugraph);
  ///
  /// typedef OrientVisitor<UGraph::UEdgeMap<bool> > Visitor;
  /// Visitor visitor(ugraph, dmap);
  ///
  /// DfsVisit<UGraph, Visitor> dfs(ugraph, visitor);
  ///
  /// dfs.run();
  ///
  /// typedef DirUGraphAdaptor<UGraph> DGraph;
  /// DGraph dgraph(ugraph, dmap);
  ///
  /// LEMON_ASSERT(countStronglyConnectedComponents(dgraph) == 
  ///              countBiEdgeConnectedComponents(ugraph), "Wrong Orientation");
  ///\endcode
  ///
  /// The number of the bi-connected components is a lower bound for
  /// the number of the strongly connected components in the directed
  /// graph because if we contract the bi-connected components to
  /// nodes we will get a tree therefore we cannot orient edges in
  /// both direction between bi-connected components. In the other way
  /// the algorithm will orient one component to be strongly
  /// connected. The two relations proof that the assertion will
  /// be always true and the found solution is optimal.
  ///
  /// \sa DirUGraphAdaptorBase
  /// \sa dirUGraphAdaptor
  template<typename _Graph, 
           typename DirectionMap = typename _Graph::template UEdgeMap<bool> > 
  class DirUGraphAdaptor : 
    public GraphAdaptorExtender<
    DirUGraphAdaptorBase<_Graph, DirectionMap> > {
  public:
    typedef _Graph Graph;
    typedef GraphAdaptorExtender<
      DirUGraphAdaptorBase<_Graph, DirectionMap> > Parent;
  protected:
    DirUGraphAdaptor() { }
  public:
    
    /// \brief Constructor of the adaptor
    ///
    /// Constructor of the adaptor
    DirUGraphAdaptor(_Graph& _graph, DirectionMap& _direction_map) { 
      setGraph(_graph);
      setDirectionMap(_direction_map);
    }
  };

  /// \brief Just gives back a DirUGraphAdaptor
  ///
  /// Just gives back a DirUGraphAdaptor
  template<typename UGraph, typename DirectionMap>
  DirUGraphAdaptor<const UGraph, DirectionMap>
  dirUGraphAdaptor(const UGraph& graph, DirectionMap& dm) {
    return DirUGraphAdaptor<const UGraph, DirectionMap>(graph, dm);
  }

  template<typename UGraph, typename DirectionMap>
  DirUGraphAdaptor<const UGraph, const DirectionMap>
  dirUGraphAdaptor(const UGraph& graph, const DirectionMap& dm) {
    return DirUGraphAdaptor<const UGraph, const DirectionMap>(graph, dm);
  }

}

#endif
