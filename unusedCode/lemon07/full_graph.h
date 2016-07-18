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

#ifndef LEMON_FULL_GRAPH_H
#define LEMON_FULL_GRAPH_H

#include <lemon/math.h>

#include <lemon/bits/base_extender.h>
#include <lemon/bits/graph_extender.h>

#include <lemon/bits/invalid.h>
#include <lemon/bits/utility.h>


///\ingroup graphs
///\file
///\brief FullGraph and FullUGraph classes.


namespace lemon {

  class FullGraphBase {
    int _nodeNum;
    int _edgeNum;
  public:

    typedef FullGraphBase Graph;

    class Node;
    class Edge;

  protected:

    FullGraphBase() {}

    void construct(int n) { _nodeNum = n; _edgeNum = n * n; }

  public:
    
    typedef True NodeNumTag;
    typedef True EdgeNumTag;

    Node operator()(int ix) const { return Node(ix); }
    int index(const Node& node) const { return node.id; }

    Edge edge(const Node& u, const Node& v) const { 
      return Edge(*this, u.id, v.id); 
    }

    int nodeNum() const { return _nodeNum; }
    int edgeNum() const { return _edgeNum; }

    int maxNodeId() const { return _nodeNum-1; }
    int maxEdgeId() const { return _edgeNum-1; }

    Node source(Edge e) const { return e.id % _nodeNum; }
    Node target(Edge e) const { return e.id / _nodeNum; }


    static int id(Node v) { return v.id; }
    static int id(Edge e) { return e.id; }

    static Node nodeFromId(int id) { return Node(id);}
    
    static Edge edgeFromId(int id) { return Edge(id);}

    typedef True FindEdgeTag;

    Edge findEdge(Node u,Node v, Edge prev = INVALID) const {
      return prev.id == -1 ? Edge(*this, u.id, v.id) : INVALID;
    }
    
      
    class Node {
      friend class FullGraphBase;

    protected:
      int id;
      Node(int _id) : id(_id) {}
    public:
      Node() {}
      Node (Invalid) : id(-1) {}
      bool operator==(const Node node) const {return id == node.id;}
      bool operator!=(const Node node) const {return id != node.id;}
      bool operator<(const Node node) const {return id < node.id;}
    };
    


    class Edge {
      friend class FullGraphBase;
      
    protected:
      int id;  // _nodeNum * target + source;

      Edge(int _id) : id(_id) {}

      Edge(const FullGraphBase& _graph, int source, int target) 
	: id(_graph._nodeNum * target+source) {}
    public:
      Edge() { }
      Edge (Invalid) { id = -1; }
      bool operator==(const Edge edge) const {return id == edge.id;}
      bool operator!=(const Edge edge) const {return id != edge.id;}
      bool operator<(const Edge edge) const {return id < edge.id;}
    };

    void first(Node& node) const {
      node.id = _nodeNum-1;
    }

    static void next(Node& node) {
      --node.id;
    }

    void first(Edge& e) const {
      e.id = _edgeNum-1;
    }

    static void next(Edge& e) {
      --e.id;
    }

    void firstOut(Edge& e, const Node& n) const {
      e.id = _edgeNum + n.id - _nodeNum;
    }

    void nextOut(Edge& e) const {
      e.id -= _nodeNum;
      if (e.id < 0) e.id = -1;
    }

    void firstIn(Edge& e, const Node& n) const {
      e.id = n.id * _nodeNum;
    }
    
    void nextIn(Edge& e) const {
      ++e.id;
      if (e.id % _nodeNum == 0) e.id = -1;
    }

  };

  typedef GraphExtender<FullGraphBase> ExtendedFullGraphBase;

  /// \ingroup graphs
  ///
  /// \brief A full graph class.
  ///
  /// This is a simple and fast directed full graph implementation.
  /// It is completely static, so you can neither add nor delete either
  /// edges or nodes.
  /// Thus it conforms to
  /// the \ref concepts::Graph "Graph" concept and
  ///it also has an
  ///important extra feature that
  ///its maps are real \ref concepts::ReferenceMap "reference map"s.
  /// \sa concepts::Graph.
  ///
  /// \sa FullUGraph
  ///
  /// \author Alpar Juttner
  class FullGraph : public ExtendedFullGraphBase {
  public:

    typedef ExtendedFullGraphBase Parent;

    /// \brief Constructor
    FullGraph() { construct(0); }

    /// \brief Constructor
    ///
    FullGraph(int n) { construct(n); }

    /// \brief Resize the graph
    ///
    /// Resize the graph. The function will fully destroy and build the graph.
    /// This cause that the maps of the graph will reallocated
    /// automatically and the previous values will be lost.
    void resize(int n) {
      Parent::notifier(Edge()).clear();
      Parent::notifier(Node()).clear();
      construct(n);
      Parent::notifier(Node()).build();
      Parent::notifier(Edge()).build();
    }

    /// \brief Returns the node with the given index.
    ///
    /// Returns the node with the given index. Because it is a
    /// static size graph the node's of the graph can be indiced
    /// by the range from 0 to \e nodeNum()-1 and the index of
    /// the node can accessed by the \e index() member.
    Node operator()(int ix) const { return Parent::operator()(ix); }

    /// \brief Returns the index of the node.
    ///
    /// Returns the index of the node. Because it is a
    /// static size graph the node's of the graph can be indiced
    /// by the range from 0 to \e nodeNum()-1 and the index of
    /// the node can accessed by the \e index() member.
    int index(const Node& node) const { return Parent::index(node); }

    /// \brief Returns the edge connects the given nodes.
    ///
    /// Returns the edge connects the given nodes.
    Edge edge(const Node& u, const Node& v) const { 
      return Parent::edge(u, v); 
    }

    /// \brief Number of nodes.
    int nodeNum() const { return Parent::nodeNum(); }
    /// \brief Number of edges.
    int edgeNum() const { return Parent::edgeNum(); }
  };


  class FullUGraphBase {
    int _nodeNum;
    int _edgeNum;
  public:

    typedef FullUGraphBase Graph;

    class Node;
    class Edge;

  protected:

    FullUGraphBase() {}

    void construct(int n) { _nodeNum = n; _edgeNum = n * (n - 1) / 2; }

  public:


    Node operator()(int ix) const { return Node(ix); }
    int index(const Node& node) const { return node.id; }

    Edge edge(const Node& u, const Node& v) const { 
      return Edge(u.id * (u.id - 1) / 2 + v.id);
    }

    typedef True NodeNumTag;
    typedef True EdgeNumTag;

    int nodeNum() const { return _nodeNum; }
    int edgeNum() const { return _edgeNum; }

    int maxNodeId() const { return _nodeNum-1; }
    int maxEdgeId() const { return _edgeNum-1; }

    static Node nodeFromId(int id) { return Node(id);}
    static Edge edgeFromId(int id) { return Edge(id);}

    Node source(Edge e) const { 
      /// \todo we may do it faster
      return Node((int(sqrt(double(1 + 8 * e.id)) + 1)) / 2);
    }

    Node target(Edge e) const { 
      int s = (int(sqrt(double(1 + 8 * e.id)) + 1)) / 2;
      return Node(e.id - s * (s - 1) / 2);
    }

    static int id(Node v) { return v.id; }
    static int id(Edge e) { return e.id; }

    Edge findEdge(Node u, Node v, Edge prev = INVALID) const {
      if (prev.id != -1 || u.id <= v.id) return Edge(-1);
      return Edge(u.id * (u.id - 1) / 2 + v.id);
    }

    typedef True FindEdgeTag;
    
      
    class Node {
      friend class FullUGraphBase;

    protected:
      int id;
      Node(int _id) { id = _id;}
    public:
      Node() {}
      Node (Invalid) { id = -1; }
      bool operator==(const Node node) const {return id == node.id;}
      bool operator!=(const Node node) const {return id != node.id;}
      bool operator<(const Node node) const {return id < node.id;}
    };
    


    class Edge {
      friend class FullUGraphBase;
      
    protected:
      int id;  // _nodeNum * target + source;

      Edge(int _id) : id(_id) {}

    public:
      Edge() { }
      Edge (Invalid) { id = -1; }
      bool operator==(const Edge edge) const {return id == edge.id;}
      bool operator!=(const Edge edge) const {return id != edge.id;}
      bool operator<(const Edge edge) const {return id < edge.id;}
    };

    void first(Node& n) const {
      n.id = _nodeNum - 1;
    }

    static void next(Node& n) {
      --n.id;
    }

    void first(Edge& e) const {
      e.id = _edgeNum - 1;
    }

    static void next(Edge& e) {
      --e.id;
    }

    void firstOut(Edge& e, const Node& n) const {      
      int src = n.id;
      int trg = 0;
      e.id = (trg < src ? src * (src - 1) / 2 + trg : -1);
    }

    /// \todo with specialized iterators we can make faster iterating
    void nextOut(Edge& e) const {
      int src = source(e).id;
      int trg = target(e).id;
      ++trg;
      e.id = (trg < src ? src * (src - 1) / 2 + trg : -1);
    }

    void firstIn(Edge& e, const Node& n) const {
      int src = n.id + 1;
      int trg = n.id;
      e.id = (src < _nodeNum ? src * (src - 1) / 2 + trg : -1);
    }
    
    void nextIn(Edge& e) const {
      int src = source(e).id;
      int trg = target(e).id;
      ++src;
      e.id = (src < _nodeNum ? src * (src - 1) / 2 + trg : -1);
    }

  };

  typedef UGraphExtender<UndirGraphExtender<FullUGraphBase> > 
  ExtendedFullUGraphBase;

  /// \ingroup graphs
  ///
  /// \brief An undirected full graph class.
  ///
  /// This is a simple and fast undirected full graph implementation.
  /// It is completely static, so you can neither add nor delete either
  /// edges or nodes.
  ///
  /// The main difference beetween the \e FullGraph and \e FullUGraph class
  /// is that this class conforms to the undirected graph concept and
  /// it does not contain the loop edges.
  ///
  ///It also has an
  ///important extra feature that
  ///its maps are real \ref concepts::ReferenceMap "reference map"s.
  ///
  /// \sa FullGraph
  ///
  /// \author Balazs Dezso
  class FullUGraph : public ExtendedFullUGraphBase {
  public:

    typedef ExtendedFullUGraphBase Parent;

    /// \brief Constructor
    FullUGraph() { construct(0); }

    /// \brief Constructor
    FullUGraph(int n) { construct(n); }

    /// \brief Resize the graph
    ///
    /// Resize the graph. The function will fully destroy and build the graph.
    /// This cause that the maps of the graph will reallocated
    /// automatically and the previous values will be lost.
    void resize(int n) {
      Parent::notifier(Edge()).clear();
      Parent::notifier(UEdge()).clear();
      Parent::notifier(Node()).clear();
      construct(n);
      Parent::notifier(Node()).build();
      Parent::notifier(UEdge()).build();
      Parent::notifier(Edge()).build();
    }

    /// \brief Returns the node with the given index.
    ///
    /// Returns the node with the given index. Because it is a
    /// static size graph the node's of the graph can be indiced
    /// by the range from 0 to \e nodeNum()-1 and the index of
    /// the node can accessed by the \e index() member.
    Node operator()(int ix) const { return Parent::operator()(ix); }

    /// \brief Returns the index of the node.
    ///
    /// Returns the index of the node. Because it is a
    /// static size graph the node's of the graph can be indiced
    /// by the range from 0 to \e nodeNum()-1 and the index of
    /// the node can accessed by the \e index() member.
    int index(const Node& node) const { return Parent::index(node); }

    /// \brief Number of nodes.
    int nodeNum() const { return Parent::nodeNum(); }
    /// \brief Number of edges.
    int edgeNum() const { return Parent::edgeNum(); }
    /// \brief Number of undirected edges.
    int uEdgeNum() const { return Parent::uEdgeNum(); }

    /// \brief Returns the edge connects the given nodes.
    ///
    /// Returns the edge connects the given nodes.
    Edge edge(const Node& u, const Node& v) const { 
      if (Parent::index(u) > Parent::index(v)) {
        return Parent::direct(Parent::edge(u, v), true);
      } else if (Parent::index(u) == Parent::index(v)) {
        return INVALID;
      } else {
        return Parent::direct(Parent::edge(v, u), false); 
      }
    }

    /// \brief Returns the undirected edge connects the given nodes.
    ///
    /// Returns the undirected edge connects the given nodes.
    UEdge uEdge(const Node& u, const Node& v) const {
      if (Parent::index(u) > Parent::index(v)) {
        return Parent::edge(u, v);
      } else if (Parent::index(u) == Parent::index(v)) {
        return INVALID;
      } else {
        return Parent::edge(v, u);
      }
    }
  };


  class FullBpUGraphBase {
  protected:

    int _aNodeNum;
    int _bNodeNum;

    int _edgeNum;

  protected:

    FullBpUGraphBase() {}

    void construct(int ann, int bnn) {
      _aNodeNum = ann;
      _bNodeNum = bnn;
      _edgeNum = ann * bnn;
    }

  public:

    class NodeSetError : public LogicError {
    public:
      virtual const char* what() const throw() { 
	return "lemon::FullBpUGraph::NodeSetError";
      }
    };
  
    class Node {
      friend class FullBpUGraphBase;
    protected:
      int id;

      Node(int _id) : id(_id) {}
    public:
      Node() {}
      Node(Invalid) { id = -1; }
      bool operator==(const Node i) const {return id==i.id;}
      bool operator!=(const Node i) const {return id!=i.id;}
      bool operator<(const Node i) const {return id<i.id;}
    };

    class UEdge {
      friend class FullBpUGraphBase;
    protected:
      int id;

      UEdge(int _id) { id = _id;}
    public:
      UEdge() {}
      UEdge (Invalid) { id = -1; }
      bool operator==(const UEdge i) const {return id==i.id;}
      bool operator!=(const UEdge i) const {return id!=i.id;}
      bool operator<(const UEdge i) const {return id<i.id;}
    };

    Node aNode(int ix) const { return Node(ix << 1); }
    Node bNode(int ix) const { return Node((ix << 1) + 1); }

    int aNodeIndex(const Node& node) const { return node.id >> 1; }
    int bNodeIndex(const Node& node) const { return node.id >> 1; }

    UEdge uEdge(const Node& u, const Node& v) const { 
      if (((u.id ^ v.id) & 1) != 1) {
        return UEdge(-1);
      } else if ((u.id & 1) == 0) {
        return UEdge((u.id >> 1) * _bNodeNum + (v.id >> 1));
      } else {
        return UEdge((v.id >> 1) * _bNodeNum + (u.id >> 1));
      }
    }

    void firstANode(Node& node) const {
      node.id = 2 * _aNodeNum - 2;
      if (node.id < 0) node.id = -1; 
    }
    void nextANode(Node& node) const {
      node.id -= 2;
      if (node.id < 0) node.id = -1; 
    }

    void firstBNode(Node& node) const {
      node.id = 2 * _bNodeNum - 1;
    }
    void nextBNode(Node& node) const {
      node.id -= 2;
    }

    void first(Node& node) const {
      if (_aNodeNum > 0) {
	node.id = 2 * _aNodeNum - 2;
      } else {
	node.id = 2 * _bNodeNum - 1;
      }
    }
    void next(Node& node) const {
      node.id -= 2;
      if (node.id == -2) {
	node.id = 2 * _bNodeNum - 1;
      }
    }
  
    void first(UEdge& edge) const {
      edge.id = _edgeNum - 1;
    }
    void next(UEdge& edge) const {
      --edge.id;
    }

    void firstFromANode(UEdge& edge, const Node& node) const {
      LEMON_ASSERT((node.id & 1) == 0, NodeSetError());
      edge.id = (node.id >> 1) * _bNodeNum;
    }
    void nextFromANode(UEdge& edge) const {
      ++(edge.id);
      if (edge.id % _bNodeNum == 0) edge.id = -1;
    }

    void firstFromBNode(UEdge& edge, const Node& node) const {
      LEMON_ASSERT((node.id & 1) == 1, NodeSetError());
      edge.id = (node.id >> 1);
    }
    void nextFromBNode(UEdge& edge) const {
      edge.id += _bNodeNum;
      if (edge.id >= _edgeNum) edge.id = -1;
    }

    static int id(const Node& node) {
      return node.id;
    }
    static Node nodeFromId(int id) {
      return Node(id);
    }
    int maxNodeId() const {
      return _aNodeNum > _bNodeNum ? 
	_aNodeNum * 2 - 2 : _bNodeNum * 2 - 1;
    }
  
    static int id(const UEdge& edge) {
      return edge.id;
    }
    static UEdge uEdgeFromId(int id) {
      return UEdge(id);
    }
    int maxUEdgeId() const {
      return _edgeNum - 1;
    }
  
    static int aNodeId(const Node& node) {
      return node.id >> 1;
    }
    static Node nodeFromANodeId(int id) {
      return Node(id << 1);
    }
    int maxANodeId() const {
      return _aNodeNum;
    }

    static int bNodeId(const Node& node) {
      return node.id >> 1;
    }
    static Node nodeFromBNodeId(int id) {
      return Node((id << 1) + 1);
    }
    int maxBNodeId() const {
      return _bNodeNum;
    }

    Node aNode(const UEdge& edge) const {
      return Node((edge.id / _bNodeNum) << 1);
    }
    Node bNode(const UEdge& edge) const {
      return Node(((edge.id % _bNodeNum) << 1) + 1);
    }

    static bool aNode(const Node& node) {
      return (node.id & 1) == 0;
    }

    static bool bNode(const Node& node) {
      return (node.id & 1) == 1;
    }


    typedef True NodeNumTag;
    int nodeNum() const { return _aNodeNum + _bNodeNum; }
    int aNodeNum() const { return _aNodeNum; }
    int bNodeNum() const { return _bNodeNum; }

    typedef True EdgeNumTag;
    int uEdgeNum() const { return _edgeNum; }


    typedef True FindEdgeTag;
    UEdge findUEdge(Node u, Node v, UEdge prev = INVALID) const {
      if (prev.id != -1 || ((u.id ^ v.id) & 1) != 1) {
        return UEdge(-1);
      } else if ((u.id & 1) == 0) {
        return UEdge((u.id >> 1) * _bNodeNum + (v.id >> 1));
      } else {
        return UEdge((v.id >> 1) * _bNodeNum + (u.id >> 1));
      }
    }

  };


  typedef BpUGraphExtender<BidirBpUGraphExtender<FullBpUGraphBase> > 
  ExtendedFullBpUGraphBase;


  /// \ingroup graphs
  ///
  /// \brief An undirected full bipartite graph class.
  ///
  /// This is a simple and fast bipartite undirected full graph implementation.
  /// It is completely static, so you can neither add nor delete either
  /// edges or nodes.
  ///
  /// \author Balazs Dezso
  class FullBpUGraph : 
    public ExtendedFullBpUGraphBase {
  public:

    typedef ExtendedFullBpUGraphBase Parent;

    FullBpUGraph() {
      Parent::construct(0, 0);
    }

    FullBpUGraph(int ann, int bnn) {
      Parent::construct(ann, bnn);
    }

    /// \brief Resize the graph
    ///
    /// Resize the graph
    void resize(int n, int m) {
      Parent::notifier(Edge()).clear();
      Parent::notifier(UEdge()).clear();
      Parent::notifier(Node()).clear();
      Parent::notifier(ANode()).clear();
      Parent::notifier(BNode()).clear();
      construct(n, m);
      Parent::notifier(ANode()).build();
      Parent::notifier(BNode()).build();
      Parent::notifier(Node()).build();
      Parent::notifier(UEdge()).build();
      Parent::notifier(Edge()).build();
    }

    /// \brief Number of nodes.
    int nodeNum() const { return Parent::nodeNum(); }
    /// \brief Number of A-nodes.
    int aNodeNum() const { return Parent::aNodeNum(); }
    /// \brief Number of B-nodes.
    int bNodeNum() const { return Parent::bNodeNum(); }
    /// \brief Number of edges.
    int edgeNum() const { return Parent::edgeNum(); }
    /// \brief Number of undirected edges.
    int uEdgeNum() const { return Parent::uEdgeNum(); }


    using Parent::aNode;
    using Parent::bNode;

    /// \brief Returns the A-node with the given index.
    ///
    /// Returns the A-node with the given index. Because it is a
    /// static size graph the node's of the graph can be indiced
    /// by the range from 0 to \e aNodeNum()-1 and the index of
    /// the node can accessed by the \e aNodeIndex() member.
    Node aNode(int ix) const { return Parent::aNode(ix); }

    /// \brief Returns the B-node with the given index.
    ///
    /// Returns the B-node with the given index. Because it is a
    /// static size graph the node's of the graph can be indiced
    /// by the range from 0 to \e bNodeNum()-1 and the index of
    /// the node can accessed by the \e bNodeIndex() member.
    Node bNode(int ix) const { return Parent::bNode(ix); }

    /// \brief Returns the index of the A-node.
    ///
    /// Returns the index of the A-node. Because it is a
    /// static size graph the node's of the graph can be indiced
    /// by the range from 0 to \e aNodeNum()-1 and the index of
    /// the node can accessed by the \e aNodeIndex() member.
    int aNodeIndex(const Node& node) const { return Parent::aNodeIndex(node); }

    /// \brief Returns the index of the B-node.
    ///
    /// Returns the index of the B-node. Because it is a
    /// static size graph the node's of the graph can be indiced
    /// by the range from 0 to \e bNodeNum()-1 and the index of
    /// the node can accessed by the \e bNodeIndex() member.
    int bNodeIndex(const Node& node) const { return Parent::bNodeIndex(node); }

    /// \brief Returns the edge connects the given nodes.
    ///
    /// Returns the edge connects the given nodes.
    Edge edge(const Node& u, const Node& v) const {
      UEdge uedge = Parent::uEdge(u, v);
      if (uedge != INVALID) {
        return Parent::direct(uedge, Parent::aNode(u));
      } else {
        return INVALID;
      }
    }

    /// \brief Returns the undirected edge connects the given nodes.
    ///
    /// Returns the undirected edge connects the given nodes.
    UEdge uEdge(const Node& u, const Node& v) const {
      return Parent::uEdge(u, v);
    }
  };

} //namespace lemon


#endif //LEMON_FULL_GRAPH_H
