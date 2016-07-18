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

#ifndef LEMON_BPUGRAPH_ADAPTOR_H
#define LEMON_BPUGRAPH_ADAPTOR_H

#include <lemon/bits/invalid.h>
#include <lemon/maps.h>

#include <lemon/bits/graph_adaptor_extender.h>

#include <lemon/bits/traits.h>

#include <iostream>

///\ingroup graph_adaptors
///\file
///\brief Several graph adaptors.
///
///This file contains several useful bpugraph adaptor functions.
///
///\author Balazs Dezso

namespace lemon {

  /// \ingroup graph_adaptors
  ///
  /// \brief Base type for the Bipartite Undirected Graph Adaptors
  ///
  /// This is the base type for most of LEMON bpugraph adaptors. 
  /// This class implements a trivial graph adaptor i.e. it only wraps the 
  /// functions and types of the graph. The purpose of this class is to 
  /// make easier implementing graph adaptors. E.g. if an adaptor is 
  /// considered which differs from the wrapped graph only in some of its 
  /// functions or types, then it can be derived from BpUGraphAdaptor, and 
  /// only the differences should be implemented.
  ///
  /// \author Balazs Dezso 
  template<typename _BpUGraph>
  class BpUGraphAdaptorBase {
  public:
    typedef _BpUGraph Graph;
    typedef Graph ParentGraph;

  protected:
    Graph* graph;

    BpUGraphAdaptorBase() : graph(0) {}

    void setGraph(Graph& _graph) { graph = &_graph; }

    Graph& getGraph() { return *graph; }
    const Graph& getGraph() const { return *graph; }

  public:

    BpUGraphAdaptorBase(Graph& _graph) : graph(&_graph) {}
 
    typedef typename Graph::Node Node;
    typedef typename Graph::ANode ANode;
    typedef typename Graph::BNode BNode;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::UEdge UEdge;
   
    void first(Node& i) const { graph->first(i); }
    void firstANode(Node& i) const { graph->firstANode(i); }
    void firstBNode(Node& i) const { graph->firstBNode(i); }
    void first(Edge& i) const { graph->first(i); }
    void first(UEdge& i) const { graph->first(i); }
    void firstIn(Edge& i, const Node& n) const { graph->firstIn(i, n); }
    void firstOut(Edge& i, const Node& n ) const { graph->firstOut(i, n); }
    void firstInc(UEdge &i, bool &d, const Node &n) const {
      graph->firstInc(i, d, n);
    }
    void firstFromANode(UEdge& i, const Node& n) const {
      graph->firstFromANode(i, n);
    }
    void firstFromBNode(UEdge& i, const Node& n) const {
      graph->firstFromBNode(i, n);
    }

    void next(Node& i) const { graph->next(i); }
    void nextANode(Node& i) const { graph->nextANode(i); }
    void nextBNode(Node& i) const { graph->nextBNode(i); }
    void next(Edge& i) const { graph->next(i); }
    void next(UEdge& i) const { graph->next(i); }
    void nextIn(Edge& i) const { graph->nextIn(i); }
    void nextOut(Edge& i) const { graph->nextOut(i); }
    void nextInc(UEdge &i, bool &d) const { graph->nextInc(i, d); }
    void nextFromANode(UEdge& i) const { graph->nextFromANode(i); }
    void nextFromBNode(UEdge& i) const { graph->nextFromBNode(i); }

    Node source(const UEdge& e) const { return graph->source(e); }
    Node target(const UEdge& e) const { return graph->target(e); }

    Node source(const Edge& e) const { return graph->source(e); }
    Node target(const Edge& e) const { return graph->target(e); }

    Node aNode(const UEdge& e) const { return graph->aNode(e); }
    Node bNode(const UEdge& e) const { return graph->bNode(e); }

    bool aNode(const Node& n) const { return graph->aNode(n); }
    bool bNode(const Node& n) const { return graph->bNode(n); }

    typedef NodeNumTagIndicator<Graph> NodeNumTag;
    int nodeNum() const { return graph->nodeNum(); }
    int aNodeNum() const { return graph->aNodeNum(); }
    int bNodeNum() const { return graph->bNodeNum(); }
    
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
  
    Node addANode() const { return graph->addANode(); }
    Node addBNode() const { return graph->addBNode(); }
    UEdge addEdge(const Node& u, const Node& v) const { 
      return graph->addEdge(u, v); 
    }

    void erase(const Node& i) const { graph->erase(i); }
    void erase(const UEdge& i) const { graph->erase(i); }
  
    void clear() const { graph->clear(); }

    bool direction(const Edge& e) const { return graph->direction(e); }
    Edge direct(const UEdge& e, bool d) const { return graph->direct(e, d); }
    
    int id(const Node& v) const { return graph->id(v); }
    int id(const ANode& v) const { return graph->id(v); }
    int id(const BNode& v) const { return graph->id(v); }
    int id(const Edge& e) const { return graph->id(e); }
    int id(const UEdge& e) const { return graph->id(e); }

    Node fromNodeId(int ix) const { return graph->fromNodeId(ix); }
    ANode nodeFromANodeId(int ix) const { return graph->nodeFromANodeId(ix); }
    BNode nodeFromBNodeId(int ix) const { return graph->nodeFromBNodeId(ix); }
    Edge fromEdgeId(int ix) const { return graph->fromEdgeId(ix); }
    UEdge fromUEdgeId(int ix) const { return graph->fromUEdgeId(ix); }

    int maxNodeId() const { return graph->maxNodeId(); }
    int maxANodeId() const { return graph->maxANodeId(); }
    int maxBNodeId() const { return graph->maxBNodeId(); }
    int maxEdgeId() const { return graph->maxEdgeId(); }
    int maxUEdgeId() const { return graph->maxEdgeId(); }

    typedef typename ItemSetTraits<Graph, Node>::ItemNotifier NodeNotifier;
    NodeNotifier& notifier(Node) const {
      return graph->notifier(Node()); 
    } 

    typedef typename ItemSetTraits<Graph, ANode>::ItemNotifier ANodeNotifier;
    ANodeNotifier& notifier(ANode) const {
      return graph->notifier(ANode());
    } 

    typedef typename ItemSetTraits<Graph, BNode>::ItemNotifier BNodeNotifier;
    BNodeNotifier& notifier(BNode) const {
      return graph->notifier(BNode());
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
      explicit NodeMap(const BpUGraphAdaptorBase<Graph>& ga) 
	: Parent(*ga.graph) {}
      NodeMap(const BpUGraphAdaptorBase<Graph>& ga, const _Value& value)
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
    class ANodeMap : public Graph::template ANodeMap<_Value> {
    public:
      typedef typename Graph::template ANodeMap<_Value> Parent;
      explicit ANodeMap(const BpUGraphAdaptorBase& ga) 
	: Parent(*ga.graph) {}
      ANodeMap(const BpUGraphAdaptorBase& ga, const _Value& value)
	: Parent(*ga.graph, value) {}

      ANodeMap& operator=(const ANodeMap& cmap) {
	return operator=<ANodeMap>(cmap);
      }

      template <typename CMap>
      ANodeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
	return *this;
      }
    };

    template <typename _Value>
    class BNodeMap : public Graph::template BNodeMap<_Value> {
    public:
      typedef typename Graph::template BNodeMap<_Value> Parent;
      explicit BNodeMap(const BpUGraphAdaptorBase<Graph>& ga) 
	: Parent(*ga.graph) {}
      BNodeMap(const BpUGraphAdaptorBase<Graph>& ga, const _Value& value)
	: Parent(*ga.graph, value) {}

      BNodeMap& operator=(const BNodeMap& cmap) {
	return operator=<BNodeMap>(cmap);
      }

      template <typename CMap>
      BNodeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
	return *this;
      }
    };

    template <typename _Value>
    class EdgeMap : public Graph::template EdgeMap<_Value> {
    public:
      typedef typename Graph::template EdgeMap<_Value> Parent;
      explicit EdgeMap(const BpUGraphAdaptorBase<Graph>& ga) 
	: Parent(*ga.graph) {}
      EdgeMap(const BpUGraphAdaptorBase<Graph>& ga, const _Value& value)
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
      explicit UEdgeMap(const BpUGraphAdaptorBase<Graph>& ga) 
	: Parent(*ga.graph) {}
      UEdgeMap(const BpUGraphAdaptorBase<Graph>& ga, const _Value& value)
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
  /// \brief Trivial Bipartite Undirected Graph Adaptor
  ///
  /// Trivial Bipartite Undirected Graph Adaptor. It does not change
  /// the functionality of the adapted graph.
  template <typename _BpUGraph>
  class BpUGraphAdaptor 
    : public BpUGraphAdaptorExtender< BpUGraphAdaptorBase<_BpUGraph> > { 
  public:
    typedef _BpUGraph Graph;
    typedef BpUGraphAdaptorExtender<BpUGraphAdaptorBase<_BpUGraph> > Parent;
  protected:
    BpUGraphAdaptor() : Parent() {}

  public:
    explicit BpUGraphAdaptor(Graph& _graph) { setGraph(_graph); }
  };

  template <typename _BpUGraph>
  class SwapBpUGraphAdaptorBase : public BpUGraphAdaptorBase<_BpUGraph> {
  public:
    
    typedef _BpUGraph Graph;
    typedef BpUGraphAdaptorBase<_BpUGraph> Parent;

  protected:
    
    SwapBpUGraphAdaptorBase() {}

  public:

    typedef typename Parent::Node Node;
    typedef typename Parent::BNode ANode;
    typedef typename Parent::ANode BNode;
    typedef typename Parent::Edge Edge;
    typedef typename Parent::UEdge UEdge;

    bool direction(const Edge& e) const { return !Parent::direction(e); }
    Edge direct(const UEdge& e, bool b) const { return Parent::direct(e, !b); }

    Node aNode(const UEdge& e) const { return Parent::bNode(e); }
    Node bNode(const UEdge& e) const { return Parent::aNode(e); }

    bool aNode(const Node& n) const { return Parent::bNode(n); }
    bool bNode(const Node& n) const { return Parent::aNode(n); }

    void firstANode(Node& i) const { Parent::firstBNode(i); }
    void firstBNode(Node& i) const { Parent::firstANode(i); }

    void nextANode(Node& i) const { Parent::nextBNode(i); }
    void nextBNode(Node& i) const { Parent::nextANode(i); }

    void firstFromANode(UEdge& i, const Node& n) const {
      Parent::firstFromBNode(i, n);
    }
    void firstFromBNode(UEdge& i, const Node& n) const {
      Parent::firstFromANode(i, n);
    }

    void nextFromANode(UEdge& i) const { Parent::nextFromBNode(i); }
    void nextFromBNode(UEdge& i) const { Parent::nextFromANode(i); }

    int id(const ANode& v) const { return Parent::id(v); }
    int id(const BNode& v) const { return Parent::id(v); }

    ANode nodeFromANodeId(int ix) const { return Parent::nodeFromBNodeId(ix); }
    BNode nodeFromBNodeId(int ix) const { return Parent::nodeFromANodeId(ix); }

    int maxANodeId() const { return Parent::maxBNodeId(); }
    int maxBNodeId() const { return Parent::maxANodeId(); }

    int aNodeNum() const { return Parent::bNodeNum(); }
    int bNodeNum() const { return Parent::aNodeNum(); }

    typedef typename Parent::BNodeNotifier ANodeNotifier;
    ANodeNotifier& notifier(ANode) const {
      return Parent::notifier(typename Parent::BNode());
    } 

    typedef typename Parent::ANodeNotifier BNodeNotifier;
    BNodeNotifier& notifier(BNode) const {
      return Parent::notifier(typename Parent::ANode());
    } 

    template <typename _Value>
    class ANodeMap : public Graph::template BNodeMap<_Value> {
    public:
      typedef typename Graph::template BNodeMap<_Value> Parent;
      explicit ANodeMap(const SwapBpUGraphAdaptorBase& ga) 
	: Parent(*ga.graph) {}
      ANodeMap(const SwapBpUGraphAdaptorBase& ga, const _Value& value)
	: Parent(*ga.graph, value) {}

      ANodeMap& operator=(const ANodeMap& cmap) {
	return operator=<ANodeMap>(cmap);
      }

      template <typename CMap>
      ANodeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
	return *this;
      }
    };

    template <typename _Value>
    class BNodeMap : public Graph::template ANodeMap<_Value> {
    public:
      typedef typename Graph::template ANodeMap<_Value> Parent;
      explicit BNodeMap(const SwapBpUGraphAdaptorBase<Graph>& ga) 
	: Parent(*ga.graph) {}
      BNodeMap(const SwapBpUGraphAdaptorBase<Graph>& ga, const _Value& value)
	: Parent(*ga.graph, value) {}

      BNodeMap& operator=(const BNodeMap& cmap) {
	return operator=<BNodeMap>(cmap);
      }

      template <typename CMap>
      BNodeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
	return *this;
      }
    };
    
  };

  /// \ingroup graph_adaptors
  ///
  /// \brief Bipartite graph adaptor which swaps the two nodeset.
  ///
  /// Bipartite graph adaptor which swaps the two nodeset. The adaptor's
  /// anode-set will be the original graph's bnode-set and the adaptor's
  /// bnode-set will be the original graph's anode-set.
  ///
  /// As example look at an algorithm what can be sped up with the
  /// swap bipartite undirected graph adaptor. If we want to find the
  /// maximum matching in the bipartite graph then it will be not changed
  /// if we swap the two nodesets. But the algorithm use the two nodeset
  /// different way. If we swap the nodesets it provides different
  /// running time. We run a test on random bipartite graphs with
  /// different rate of the anode-set size and bnode-set size.
  /// We always made graphs with 10000 nodes and 20000 edges and we
  /// computed the maximum cardinality matching with the Hopcroft-Karp 
  /// algorithm.
  ///
  /// The next diagram shows the running time of the tests. If the anode-set
  /// is smaller than the bnode-set the running time is better than with 
  /// the swapped graph. Other conclusion is that the running time
  /// is greater when the two nodesets size are nearly equal. 
  ///
  /// \image html swap_test.png
  /// \image latex swap_test.eps "Swapping nodesets" width=\textwidth
  /// 
  /// The next part shows how can we swap the two nodeset:
  ///
  ///\code
  /// typedef SwapBpUGraphAdaptor<BpUGraph> SBpUGraph;
  /// SBpUGraph sbpugraph(bpugraph);
  /// MaxBipartiteMatching<SBpUGraph> sbpumatch(sbpugraph);
  ///\endcode
  template <typename _BpUGraph>
  class SwapBpUGraphAdaptor 
    : public BpUGraphAdaptorExtender<SwapBpUGraphAdaptorBase<_BpUGraph> > { 
  public:

    typedef _BpUGraph Graph;
    typedef BpUGraphAdaptorExtender<SwapBpUGraphAdaptorBase<_BpUGraph> > 
    Parent;

  protected:
    SwapBpUGraphAdaptor() : Parent() {}

  public:

    /// \brief Construct a swapped graph.
    ///
    /// Construct a swapped graph.
    explicit SwapBpUGraphAdaptor(Graph& _graph) { setGraph(_graph); }

  };


  template <typename _BpUGraph, 
            typename _ANMatchingMap, typename _BNMatchingMap>
  class MatchingBpUGraphAdaptorBase 
    : public BpUGraphAdaptorBase<const _BpUGraph>
  {
  public:
    
    typedef _BpUGraph Graph;
    typedef _ANMatchingMap ANodeMatchingMap;
    typedef _BNMatchingMap BNodeMatchingMap;

    typedef BpUGraphAdaptorBase<const _BpUGraph> Parent;

  protected:
    
    MatchingBpUGraphAdaptorBase() {}

    void setANodeMatchingMap(ANodeMatchingMap& _anode_matching) {
      anode_matching = &_anode_matching;
    }

    void setBNodeMatchingMap(BNodeMatchingMap& _bnode_matching) {
      bnode_matching = &_bnode_matching;
    }

  public:

    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;

    void firstOut(Edge& edge, const Node& node) const {
      if (Parent::aNode(node)) {
        Parent::firstOut(edge, node);
        if (edge != INVALID && (*anode_matching)[node] == edge) {
          Parent::nextOut(edge);
        }
      } else {
        edge = (*bnode_matching)[node] != INVALID ?
          Parent::direct((*bnode_matching)[node], false) : INVALID;
      }
    }

    void nextOut(Edge& edge) const {
      if (Parent::aNode(Parent::source(edge))) {
        Parent::nextOut(edge);
        if (edge != INVALID && (*anode_matching)[Parent::aNode(edge)] == edge) {
          Parent::nextOut(edge);
        }
      } else {
        edge = INVALID;
      }
    }
 
    void firstIn(Edge& edge, const Node& node) const {
      if (Parent::aNode(node)) {
        edge = (*bnode_matching)[node] != INVALID ?
          Parent::direct((*anode_matching)[node], false) : INVALID;
      } else {
        Parent::firstIn(edge, node);
        if (edge != INVALID && (*bnode_matching)[node] == edge) {
          Parent::nextIn(edge);
        }
      }
    }

    void nextIn(Edge& edge) const {
      if (Parent::aNode(Parent::target(edge))) {
        edge = INVALID;
      } else {
        Parent::nextIn(edge);
        if (edge != INVALID && (*bnode_matching)[Parent::bNode(edge)] == edge) {
          Parent::nextIn(edge);
        }
      }
    } 

  private:
    ANodeMatchingMap* anode_matching;
    BNodeMatchingMap* bnode_matching;
  };


  /// \ingroup graph_adaptors
  ///
  /// \brief Bipartite graph adaptor to implement matching algorithms.
  ///
  /// Bipartite graph adaptor to implement matchings. It implements
  /// the residual graph of the matching.
  template <typename _BpUGraph, 
            typename _ANMatchingMap = 
            typename _BpUGraph::template ANodeMap<typename _BpUGraph::UEdge>, 
            typename _BNMatchingMap =
            typename _BpUGraph::template BNodeMap<typename _BpUGraph::UEdge> >
  class MatchingBpUGraphAdaptor 
    : public BpUGraphAdaptorExtender<
    MatchingBpUGraphAdaptorBase<_BpUGraph, _ANMatchingMap, _BNMatchingMap> > 
  { 
  public:

    typedef _BpUGraph BpUGraph;
    typedef _BpUGraph Graph;
    typedef _ANMatchingMap ANodeMatchingMap;
    typedef _BNMatchingMap BNodeMatchingMap;
    typedef BpUGraphAdaptorExtender<
      MatchingBpUGraphAdaptorBase<BpUGraph, 
                                  ANodeMatchingMap, BNodeMatchingMap> > 
    Parent;

  protected:
    MatchingBpUGraphAdaptor() : Parent() {}

  public:

    MatchingBpUGraphAdaptor(const Graph& _graph, 
                            ANodeMatchingMap& _anode_matching,
                            BNodeMatchingMap& _bnode_matching) {
      setGraph(_graph);
      setANodeMatchingMap(_anode_matching);
      setBNodeMatchingMap(_bnode_matching);
    }

  };

}

#endif
