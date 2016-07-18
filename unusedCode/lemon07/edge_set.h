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

#ifndef LEMON_EDGE_SET_H
#define LEMON_EDGE_SET_H


#include <lemon/bits/default_map.h>
#include <lemon/bits/edge_set_extender.h>

/// \ingroup semi_adaptors
/// \file
/// \brief EdgeSet classes.
///
/// Graphs which use another graph's node-set as own. 

namespace lemon {

  template <typename _Graph>
  class ListEdgeSetBase {
  public:

    typedef _Graph Graph;
    typedef typename Graph::Node Node;
    typedef typename Graph::NodeIt NodeIt;

  protected:

    struct NodeT {
      int first_out, first_in;
      NodeT() : first_out(-1), first_in(-1) {}
    };

    typedef DefaultMap<Graph, Node, NodeT> NodesImplBase;

    NodesImplBase* nodes;

    struct EdgeT {
      Node source, target;
      int next_out, next_in;
      int prev_out, prev_in;
      EdgeT() : prev_out(-1), prev_in(-1) {}
    };

    std::vector<EdgeT> edges;

    int first_edge;
    int first_free_edge;

    const Graph* graph;

    void initalize(const Graph& _graph, NodesImplBase& _nodes) {
      graph = &_graph;
      nodes = &_nodes;
    }
    
  public:

    class Edge {
      friend class ListEdgeSetBase<Graph>;
    protected:
      Edge(int _id) : id(_id) {}
      int id;
    public:
      Edge() {}
      Edge(Invalid) : id(-1) {}
      bool operator==(const Edge& edge) const { return id == edge.id; }
      bool operator!=(const Edge& edge) const { return id != edge.id; }
      bool operator<(const Edge& edge) const { return id < edge.id; }
    };

    ListEdgeSetBase() : first_edge(-1), first_free_edge(-1) {} 

    Edge addEdge(const Node& u, const Node& v) {
      int n;
      if (first_free_edge == -1) {
	n = edges.size();
	edges.push_back(EdgeT());
      } else {
	n = first_free_edge;
	first_free_edge = edges[first_free_edge].next_in;
      }
      edges[n].next_in = (*nodes)[v].first_in;
      if ((*nodes)[v].first_in != -1) {
        edges[(*nodes)[v].first_in].prev_in = n;
      }
      (*nodes)[v].first_in = n;
      edges[n].next_out = (*nodes)[u].first_out;
      if ((*nodes)[u].first_out != -1) {
        edges[(*nodes)[u].first_out].prev_out = n;
      }
      (*nodes)[u].first_out = n;
      edges[n].source = u;
      edges[n].target = v;
      return Edge(n);
    }

    void erase(const Edge& edge) {
      int n = edge.id;
      if (edges[n].prev_in != -1) {
	edges[edges[n].prev_in].next_in = edges[n].next_in;
      } else {
	(*nodes)[edges[n].target].first_in = edges[n].next_in;
      }
      if (edges[n].next_in != -1) {
	edges[edges[n].next_in].prev_in = edges[n].prev_in;
      }

      if (edges[n].prev_out != -1) {
	edges[edges[n].prev_out].next_out = edges[n].next_out;
      } else {
	(*nodes)[edges[n].source].first_out = edges[n].next_out;
      }
      if (edges[n].next_out != -1) {
	edges[edges[n].next_out].prev_out = edges[n].prev_out;
      }
           
    }

    void clear() {
      Node node;
      for (first(node); node != INVALID; next(node)) {
        (*nodes)[node].first_in = -1;
        (*nodes)[node].first_out = -1;
      }
      edges.clear();
      first_edge = -1;
      first_free_edge = -1;
    }

    void first(Node& node) const {
      graph->first(node);
    }

    void next(Node& node) const {
      graph->next(node);
    }

    void first(Edge& edge) const {
      Node node;
      for (first(node); node != INVALID && (*nodes)[node].first_in == -1; 
	   next(node));
      edge.id = (node == INVALID) ? -1 : (*nodes)[node].first_in;
    }

    void next(Edge& edge) const {
      if (edges[edge.id].next_in != -1) {
	edge.id = edges[edge.id].next_in;
      } else {
	Node node = edges[edge.id].target;
	for (next(node); node != INVALID && (*nodes)[node].first_in == -1; 
	     next(node));
	edge.id = (node == INVALID) ? -1 : (*nodes)[node].first_in;
      }      
    }

    void firstOut(Edge& edge, const Node& node) const {
      edge.id = (*nodes)[node].first_out;    
    }
    
    void nextOut(Edge& edge) const {
      edge.id = edges[edge.id].next_out;        
    }

    void firstIn(Edge& edge, const Node& node) const {
      edge.id = (*nodes)[node].first_in;          
    }

    void nextIn(Edge& edge) const {
      edge.id = edges[edge.id].next_in;    
    }

    int id(const Node& node) const { return graph->id(node); }
    int id(const Edge& edge) const { return edge.id; }

    Node nodeFromId(int ix) const { return graph->nodeFromId(ix); }
    Edge edgeFromId(int ix) const { return Edge(ix); }

    int maxNodeId() const { return graph->maxNodeId(); };
    int maxEdgeId() const { return edges.size() - 1; }

    Node source(const Edge& edge) const { return edges[edge.id].source;}
    Node target(const Edge& edge) const { return edges[edge.id].target;}

    typedef typename ItemSetTraits<Graph, Node>::ItemNotifier NodeNotifier;

    NodeNotifier& notifier(Node) const {
      return graph->notifier(Node());
    } 

    template <typename _Value>
    class NodeMap : public Graph::template NodeMap<_Value> {
    public:

      typedef typename _Graph::template NodeMap<_Value> Parent;

      explicit NodeMap(const ListEdgeSetBase<Graph>& edgeset) 
	: Parent(*edgeset.graph) {}

      NodeMap(const ListEdgeSetBase<Graph>& edgeset, const _Value& value)
	: Parent(*edgeset.graph, value) {}

      NodeMap& operator=(const NodeMap& cmap) {
        return operator=<NodeMap>(cmap);
      }

      template <typename CMap>
      NodeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
        return *this;
      }
    };

  };

  /// \ingroup semi_adaptors
  ///
  /// \brief Graph using a node set of another graph and an
  /// own edge set.
  ///
  /// This structure can be used to establish another graph over a node set
  /// of an existing one. The node iterator will go through the nodes of the
  /// original graph.
  ///
  /// \param _Graph The type of the graph which shares its node set with 
  /// this class. Its interface must conform to the \ref concepts::Graph
  /// "Graph" concept.
  ///
  template <typename _Graph>
  class ListEdgeSet : public EdgeSetExtender<ListEdgeSetBase<_Graph> > {

  public:

    typedef EdgeSetExtender<ListEdgeSetBase<_Graph> > Parent;
    
    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;
    
    typedef _Graph Graph;


    typedef typename Parent::NodesImplBase NodesImplBase;

    void eraseNode(const Node& node) {
      Edge edge;
      Parent::firstOut(edge, node);
      while (edge != INVALID ) {
	erase(edge);
	Parent::firstOut(edge, node);
      } 

      Parent::firstIn(edge, node);
      while (edge != INVALID ) {
	erase(edge);
	Parent::firstIn(edge, node);
      }
    }
    
    void clearNodes() {
      Parent::clear();
    }

    class NodesImpl : public NodesImplBase {
    public:
      typedef NodesImplBase Parent;
      
      NodesImpl(const Graph& graph, ListEdgeSet& edgeset) 
	: Parent(graph), _edgeset(edgeset) {}

      virtual ~NodesImpl() {}
      
    protected:

      virtual void erase(const Node& node) {
	_edgeset.eraseNode(node);
	Parent::erase(node);
      }
      virtual void erase(const std::vector<Node>& nodes) {
        for (int i = 0; i < int(nodes.size()); ++i) {
          _edgeset.eraseNode(nodes[i]);
        }
	Parent::erase(nodes);
      }
      virtual void clear() {
	_edgeset.clearNodes();
	Parent::clear();
      }

    private:
      ListEdgeSet& _edgeset;
    };

    NodesImpl nodes;
    
  public:

    /// \brief Constructor of the adaptor.
    /// 
    /// Constructor of the adaptor.
    ListEdgeSet(const Graph& graph) : nodes(graph, *this) {
      Parent::initalize(graph, nodes);
    }
    
  };

  template <typename _Graph>
  class ListUEdgeSetBase {
  public:

    typedef _Graph Graph;
    typedef typename Graph::Node Node;
    typedef typename Graph::NodeIt NodeIt;

  protected:

    struct NodeT {
      int first_out;
      NodeT() : first_out(-1) {}
    };

    typedef DefaultMap<Graph, Node, NodeT> NodesImplBase;

    NodesImplBase* nodes;

    struct EdgeT {
      Node target;
      int prev_out, next_out;
      EdgeT() : prev_out(-1), next_out(-1) {}
    };

    std::vector<EdgeT> edges;

    int first_edge;
    int first_free_edge;

    const Graph* graph;

    void initalize(const Graph& _graph, NodesImplBase& _nodes) {
      graph = &_graph;
      nodes = &_nodes;
    }
    
  public:

    class UEdge {
      friend class ListUEdgeSetBase;
    protected:

      int id;
      explicit UEdge(int _id) { id = _id;}

    public:
      UEdge() {}
      UEdge (Invalid) { id = -1; }
      bool operator==(const UEdge& edge) const {return id == edge.id;}
      bool operator!=(const UEdge& edge) const {return id != edge.id;}
      bool operator<(const UEdge& edge) const {return id < edge.id;}
    };

    class Edge {
      friend class ListUEdgeSetBase;
    protected:
      Edge(int _id) : id(_id) {}
      int id;
    public:
      operator UEdge() const { return uEdgeFromId(id / 2); }

      Edge() {}
      Edge(Invalid) : id(-1) {}
      bool operator==(const Edge& edge) const { return id == edge.id; }
      bool operator!=(const Edge& edge) const { return id != edge.id; }
      bool operator<(const Edge& edge) const { return id < edge.id; }
    };

    ListUEdgeSetBase() : first_edge(-1), first_free_edge(-1) {} 

    UEdge addEdge(const Node& u, const Node& v) {
      int n;

      if (first_free_edge == -1) {
	n = edges.size();
	edges.push_back(EdgeT());
	edges.push_back(EdgeT());
      } else {
	n = first_free_edge;
	first_free_edge = edges[n].next_out;
      }
      
      edges[n].target = u;
      edges[n | 1].target = v;

      edges[n].next_out = (*nodes)[v].first_out;
      if ((*nodes)[v].first_out != -1) {
	edges[(*nodes)[v].first_out].prev_out = n;
      }
      (*nodes)[v].first_out = n;
      edges[n].prev_out = -1;
      
      if ((*nodes)[u].first_out != -1) {
	edges[(*nodes)[u].first_out].prev_out = (n | 1);
      }
      edges[n | 1].next_out = (*nodes)[u].first_out;
      (*nodes)[u].first_out = (n | 1);
      edges[n | 1].prev_out = -1;

      return UEdge(n / 2);
    }

    void erase(const UEdge& edge) {
      int n = edge.id * 2;
      
      if (edges[n].next_out != -1) {
	edges[edges[n].next_out].prev_out = edges[n].prev_out;
      } 

      if (edges[n].prev_out != -1) {
	edges[edges[n].prev_out].next_out = edges[n].next_out;
      } else {
	(*nodes)[edges[n | 1].target].first_out = edges[n].next_out;
      }

      if (edges[n | 1].next_out != -1) {
	edges[edges[n | 1].next_out].prev_out = edges[n | 1].prev_out;
      } 

      if (edges[n | 1].prev_out != -1) {
	edges[edges[n | 1].prev_out].next_out = edges[n | 1].next_out;
      } else {
	(*nodes)[edges[n].target].first_out = edges[n | 1].next_out;
      }
      
      edges[n].next_out = first_free_edge;
      first_free_edge = n;      
           
    }

    void clear() {
      Node node;
      for (first(node); node != INVALID; next(node)) {
        (*nodes)[node].first_out = -1;
      }
      edges.clear();
      first_edge = -1;
      first_free_edge = -1;
    }

    void first(Node& node) const {
      graph->first(node);
    }

    void next(Node& node) const {
      graph->next(node);
    }

    void first(Edge& edge) const {
      Node node;
      first(node);
      while (node != INVALID && (*nodes)[node].first_out == -1) {
        next(node);
      }
      edge.id = (node == INVALID) ? -1 : (*nodes)[node].first_out;
    }

    void next(Edge& edge) const {
      if (edges[edge.id].next_out != -1) {
	edge.id = edges[edge.id].next_out;
      } else {
	Node node = edges[edge.id ^ 1].target;
	next(node);
        while(node != INVALID && (*nodes)[node].first_out == -1) {
          next(node);
        }
	edge.id = (node == INVALID) ? -1 : (*nodes)[node].first_out;
      }      
    }

    void first(UEdge& uedge) const {
      Node node;
      first(node);
      while (node != INVALID) {
        uedge.id = (*nodes)[node].first_out;
        while ((uedge.id & 1) != 1) {
          uedge.id = edges[uedge.id].next_out;
        }
        if (uedge.id != -1) {
          uedge.id /= 2;
          return;
        } 
        next(node);
      }
      uedge.id = -1;
    }

    void next(UEdge& uedge) const {
      Node node = edges[uedge.id * 2].target;
      uedge.id = edges[(uedge.id * 2) | 1].next_out;
      while ((uedge.id & 1) != 1) {
        uedge.id = edges[uedge.id].next_out;
      }
      if (uedge.id != -1) {
        uedge.id /= 2;
        return;
      } 
      next(node);
      while (node != INVALID) {
        uedge.id = (*nodes)[node].first_out;
        while ((uedge.id & 1) != 1) {
          uedge.id = edges[uedge.id].next_out;
        }
        if (uedge.id != -1) {
          uedge.id /= 2;
          return;
        } 
	next(node);
      }
      uedge.id = -1;
    }

    void firstOut(Edge& edge, const Node& node) const {
      edge.id = (*nodes)[node].first_out;
    }
    
    void nextOut(Edge& edge) const {
      edge.id = edges[edge.id].next_out;        
    }

    void firstIn(Edge& edge, const Node& node) const {
      edge.id = (((*nodes)[node].first_out) ^ 1);
      if (edge.id == -2) edge.id = -1;
    }

    void nextIn(Edge& edge) const {
      edge.id = ((edges[edge.id ^ 1].next_out) ^ 1);
      if (edge.id == -2) edge.id = -1;
    }

    void firstInc(UEdge &edge, bool& dir, const Node& node) const {
      int de = (*nodes)[node].first_out;
      if (de != -1 ) {
        edge.id = de / 2;
        dir = ((de & 1) == 1);
      } else {
        edge.id = -1;
        dir = true;
      }
    }
    void nextInc(UEdge &edge, bool& dir) const {
      int de = (edges[(edge.id * 2) | (dir ? 1 : 0)].next_out);
      if (de != -1 ) {
        edge.id = de / 2;
        dir = ((de & 1) == 1);
      } else {
        edge.id = -1;
        dir = true;
      }
    }

    static bool direction(Edge edge) {
      return (edge.id & 1) == 1;
    }

    static Edge direct(UEdge uedge, bool dir) {
      return Edge(uedge.id * 2 + (dir ? 1 : 0));
    }

    int id(const Node& node) const { return graph->id(node); }
    static int id(Edge e) { return e.id; }
    static int id(UEdge e) { return e.id; }

    Node nodeFromId(int id) const { return graph->nodeFromId(id); }
    static Edge edgeFromId(int id) { return Edge(id);}
    static UEdge uEdgeFromId(int id) { return UEdge(id);}

    int maxNodeId() const { return graph->maxNodeId(); };
    int maxUEdgeId() const { return edges.size() / 2 - 1; }
    int maxEdgeId() const { return edges.size()-1; }

    Node source(Edge e) const { return edges[e.id ^ 1].target; }
    Node target(Edge e) const { return edges[e.id].target; }

    Node source(UEdge e) const { return edges[2 * e.id].target; }
    Node target(UEdge e) const { return edges[2 * e.id + 1].target; }

    typedef typename ItemSetTraits<Graph, Node>::ItemNotifier NodeNotifier;

    NodeNotifier& notifier(Node) const {
      return graph->notifier(Node());
    } 

    template <typename _Value>
    class NodeMap : public Graph::template NodeMap<_Value> {
    public:

      typedef typename _Graph::template NodeMap<_Value> Parent;

      explicit NodeMap(const ListUEdgeSetBase<Graph>& edgeset) 
	: Parent(*edgeset.graph) {}

      NodeMap(const ListUEdgeSetBase<Graph>& edgeset, const _Value& value)
	: Parent(*edgeset.graph, value) {}

      NodeMap& operator=(const NodeMap& cmap) {
        return operator=<NodeMap>(cmap);
      }

      template <typename CMap>
      NodeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
        return *this;
      }
    };

  };

  /// \ingroup semi_adaptors
  ///
  /// \brief Graph using a node set of another graph and an
  /// own uedge set.
  ///
  /// This structure can be used to establish another graph over a node set
  /// of an existing one. The node iterator will go through the nodes of the
  /// original graph.
  ///
  /// \param _Graph The type of the graph which shares its node set with 
  /// this class. Its interface must conform to the \ref concepts::Graph
  /// "Graph" concept.
  ///
  /// In the edge extension and removing it conforms to the 
  /// \ref concepts::UGraph "UGraph" concept.
  template <typename _Graph>
  class ListUEdgeSet : public UEdgeSetExtender<ListUEdgeSetBase<_Graph> > {

  public:

    typedef UEdgeSetExtender<ListUEdgeSetBase<_Graph> > Parent;
    
    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;
    
    typedef _Graph Graph;


    typedef typename Parent::NodesImplBase NodesImplBase;

    void eraseNode(const Node& node) {
      Edge edge;
      Parent::firstOut(edge, node);
      while (edge != INVALID ) {
	erase(edge);
	Parent::firstOut(edge, node);
      } 

    }
    
    void clearNodes() {
      Parent::clear();
    }

    class NodesImpl : public NodesImplBase {
    public:
      typedef NodesImplBase Parent;
      
      NodesImpl(const Graph& graph, ListUEdgeSet& edgeset) 
	: Parent(graph), _edgeset(edgeset) {}

      virtual ~NodesImpl() {}
      
    protected:

      virtual void erase(const Node& node) {
	_edgeset.eraseNode(node);
	Parent::erase(node);
      }
      virtual void erase(const std::vector<Node>& nodes) {
	for (int i = 0; i < int(nodes.size()); ++i) {
	  _edgeset.eraseNode(nodes[i]);
	}
	Parent::erase(nodes);
      }
      virtual void clear() {
	_edgeset.clearNodes();
	Parent::clear();
      }

    private:
      ListUEdgeSet& _edgeset;
    };

    NodesImpl nodes;
    
  public:

    /// \brief Constructor of the adaptor.
    /// 
    /// Constructor of the adaptor.
    ListUEdgeSet(const Graph& graph) : nodes(graph, *this) {
      Parent::initalize(graph, nodes);
    }
    
  };

  template <typename _Graph>
  class SmartEdgeSetBase {
  public:

    typedef _Graph Graph;
    typedef typename Graph::Node Node;
    typedef typename Graph::NodeIt NodeIt;

  protected:

    struct NodeT {
      int first_out, first_in;
      NodeT() : first_out(-1), first_in(-1) {}
    };

    typedef DefaultMap<Graph, Node, NodeT> NodesImplBase;

    NodesImplBase* nodes;

    struct EdgeT {
      Node source, target;
      int next_out, next_in;
      EdgeT() {}
    };

    std::vector<EdgeT> edges;

    const Graph* graph;

    void initalize(const Graph& _graph, NodesImplBase& _nodes) {
      graph = &_graph;
      nodes = &_nodes;
    }
    
  public:

    class Edge {
      friend class SmartEdgeSetBase<Graph>;
    protected:
      Edge(int _id) : id(_id) {}
      int id;
    public:
      Edge() {}
      Edge(Invalid) : id(-1) {}
      bool operator==(const Edge& edge) const { return id == edge.id; }
      bool operator!=(const Edge& edge) const { return id != edge.id; }
      bool operator<(const Edge& edge) const { return id < edge.id; }
    };

    SmartEdgeSetBase() {} 

    Edge addEdge(const Node& u, const Node& v) {
      int n = edges.size();
      edges.push_back(EdgeT());
      edges[n].next_in = (*nodes)[v].first_in;
      (*nodes)[v].first_in = n;
      edges[n].next_out = (*nodes)[u].first_out;
      (*nodes)[u].first_out = n;
      edges[n].source = u;
      edges[n].target = v;
      return Edge(n);
    }

    void clear() {
      Node node;
      for (first(node); node != INVALID; next(node)) {
        (*nodes)[node].first_in = -1;
        (*nodes)[node].first_out = -1;
      }
      edges.clear();
    }

    void first(Node& node) const {
      graph->first(node);
    }

    void next(Node& node) const {
      graph->next(node);
    }

    void first(Edge& edge) const {
      edge.id = edges.size() - 1;
    }

    void next(Edge& edge) const {
      --edge.id;
    }

    void firstOut(Edge& edge, const Node& node) const {
      edge.id = (*nodes)[node].first_out;    
    }
    
    void nextOut(Edge& edge) const {
      edge.id = edges[edge.id].next_out;        
    }

    void firstIn(Edge& edge, const Node& node) const {
      edge.id = (*nodes)[node].first_in;          
    }

    void nextIn(Edge& edge) const {
      edge.id = edges[edge.id].next_in;    
    }

    int id(const Node& node) const { return graph->id(node); }
    int id(const Edge& edge) const { return edge.id; }

    Node nodeFromId(int ix) const { return graph->nodeFromId(ix); }
    Edge edgeFromId(int ix) const { return Edge(ix); }

    int maxNodeId() const { return graph->maxNodeId(); };
    int maxEdgeId() const { return edges.size() - 1; }

    Node source(const Edge& edge) const { return edges[edge.id].source;}
    Node target(const Edge& edge) const { return edges[edge.id].target;}

    typedef typename ItemSetTraits<Graph, Node>::ItemNotifier NodeNotifier;

    NodeNotifier& notifier(Node) const {
      return graph->notifier(Node());
    } 

    template <typename _Value>
    class NodeMap : public Graph::template NodeMap<_Value> {
    public:

      typedef typename _Graph::template NodeMap<_Value> Parent;

      explicit NodeMap(const SmartEdgeSetBase<Graph>& edgeset) 
	: Parent(*edgeset.graph) { }

      NodeMap(const SmartEdgeSetBase<Graph>& edgeset, const _Value& value)
	: Parent(*edgeset.graph, value) { }

      NodeMap& operator=(const NodeMap& cmap) {
        return operator=<NodeMap>(cmap);
      }

      template <typename CMap>
      NodeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
        return *this;
      }
    };

  };


  /// \ingroup semi_adaptors
  ///
  /// \brief Graph using a node set of another graph and an
  /// own edge set.
  ///
  /// This structure can be used to establish another graph over a node set
  /// of an existing one. The node iterator will go through the nodes of the
  /// original graph.
  ///
  /// \param _Graph The type of the graph which shares its node set with 
  /// this class. Its interface must conform to the \ref concepts::Graph
  /// "Graph" concept.
  ///
  /// In the edge extension and removing it conforms to the 
  /// \ref concepts::Graph "Graph" concept.
  template <typename _Graph>
  class SmartEdgeSet : public EdgeSetExtender<SmartEdgeSetBase<_Graph> > {

  public:

    typedef EdgeSetExtender<SmartEdgeSetBase<_Graph> > Parent;
    
    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;
    
    typedef _Graph Graph;

  protected:

    typedef typename Parent::NodesImplBase NodesImplBase;

    void eraseNode(const Node& node) {
      if (Parent::InEdgeIt(*this, node) == INVALID &&
          Parent::OutEdgeIt(*this, node) == INVALID) {
        return;
      }
      throw typename NodesImplBase::Notifier::ImmediateDetach();
    }
    
    void clearNodes() {
      Parent::clear();
    }

    class NodesImpl : public NodesImplBase {
    public:
      typedef NodesImplBase Parent;
      
      NodesImpl(const Graph& graph, SmartEdgeSet& edgeset) 
	: Parent(graph), _edgeset(edgeset) {}

      virtual ~NodesImpl() {}
      
      bool attached() const {
        return Parent::attached();
      }

    protected:

      virtual void erase(const Node& node) {
        try {
          _edgeset.eraseNode(node);
          Parent::erase(node);
        } catch (const typename NodesImplBase::Notifier::ImmediateDetach&) {
          Parent::clear();
          throw;
        }
      }
      virtual void erase(const std::vector<Node>& nodes) {
        try {
          for (int i = 0; i < int(nodes.size()); ++i) {
            _edgeset.eraseNode(nodes[i]);
          }
          Parent::erase(nodes);
        } catch (const typename NodesImplBase::Notifier::ImmediateDetach&) {
          Parent::clear();
          throw;
        }
      }
      virtual void clear() {
	_edgeset.clearNodes();
	Parent::clear();
      }

    private:
      SmartEdgeSet& _edgeset;
    };

    NodesImpl nodes;
    
  public:

    /// \brief Constructor of the adaptor.
    /// 
    /// Constructor of the adaptor.
    SmartEdgeSet(const Graph& graph) : nodes(graph, *this) {
      Parent::initalize(graph, nodes);
    }

    bool valid() const {
      return nodes.attached();
    }
    
  };


  template <typename _Graph>
  class SmartUEdgeSetBase {
  public:

    typedef _Graph Graph;
    typedef typename Graph::Node Node;
    typedef typename Graph::NodeIt NodeIt;

  protected:

    struct NodeT {
      int first_out;
      NodeT() : first_out(-1) {}
    };

    typedef DefaultMap<Graph, Node, NodeT> NodesImplBase;

    NodesImplBase* nodes;

    struct EdgeT {
      Node target;
      int next_out;
      EdgeT() {}
    };

    std::vector<EdgeT> edges;

    const Graph* graph;

    void initalize(const Graph& _graph, NodesImplBase& _nodes) {
      graph = &_graph;
      nodes = &_nodes;
    }
    
  public:

    class UEdge {
      friend class SmartUEdgeSetBase;
    protected:

      int id;
      explicit UEdge(int _id) { id = _id;}

    public:
      UEdge() {}
      UEdge (Invalid) { id = -1; }
      bool operator==(const UEdge& edge) const {return id == edge.id;}
      bool operator!=(const UEdge& edge) const {return id != edge.id;}
      bool operator<(const UEdge& edge) const {return id < edge.id;}
    };

    class Edge {
      friend class SmartUEdgeSetBase;
    protected:
      Edge(int _id) : id(_id) {}
      int id;
    public:
      operator UEdge() const { return uEdgeFromId(id / 2); }

      Edge() {}
      Edge(Invalid) : id(-1) {}
      bool operator==(const Edge& edge) const { return id == edge.id; }
      bool operator!=(const Edge& edge) const { return id != edge.id; }
      bool operator<(const Edge& edge) const { return id < edge.id; }
    };

    SmartUEdgeSetBase() {} 

    UEdge addEdge(const Node& u, const Node& v) {
      int n = edges.size();
      edges.push_back(EdgeT());
      edges.push_back(EdgeT());
      
      edges[n].target = u;
      edges[n | 1].target = v;

      edges[n].next_out = (*nodes)[v].first_out;
      (*nodes)[v].first_out = n;

      edges[n | 1].next_out = (*nodes)[u].first_out;	
      (*nodes)[u].first_out = (n | 1);

      return UEdge(n / 2);
    }

    void clear() {
      Node node;
      for (first(node); node != INVALID; next(node)) {
        (*nodes)[node].first_out = -1;
      }
      edges.clear();
    }

    void first(Node& node) const {
      graph->first(node);
    }

    void next(Node& node) const {
      graph->next(node);
    }

    void first(Edge& edge) const { 
      edge.id = edges.size() - 1;
    }

    void next(Edge& edge) const {
      --edge.id;
    }

    void first(UEdge& edge) const { 
      edge.id = edges.size() / 2 - 1;
    }

    void next(UEdge& edge) const {
      --edge.id;
    }

    void firstOut(Edge& edge, const Node& node) const {
      edge.id = (*nodes)[node].first_out;    
    }
    
    void nextOut(Edge& edge) const {
      edge.id = edges[edge.id].next_out;        
    }

    void firstIn(Edge& edge, const Node& node) const {
      edge.id = (((*nodes)[node].first_out) ^ 1);
      if (edge.id == -2) edge.id = -1;
    }

    void nextIn(Edge& edge) const {
      edge.id = ((edges[edge.id ^ 1].next_out) ^ 1);
      if (edge.id == -2) edge.id = -1;
    }

    void firstInc(UEdge &edge, bool& dir, const Node& node) const {
      int de = (*nodes)[node].first_out;
      if (de != -1 ) {
        edge.id = de / 2;
        dir = ((de & 1) == 1);
      } else {
        edge.id = -1;
        dir = true;
      }
    }
    void nextInc(UEdge &edge, bool& dir) const {
      int de = (edges[(edge.id * 2) | (dir ? 1 : 0)].next_out);
      if (de != -1 ) {
        edge.id = de / 2;
        dir = ((de & 1) == 1);
      } else {
        edge.id = -1;
        dir = true;
      }
    }

    static bool direction(Edge edge) {
      return (edge.id & 1) == 1;
    }

    static Edge direct(UEdge uedge, bool dir) {
      return Edge(uedge.id * 2 + (dir ? 1 : 0));
    }

    int id(Node node) const { return graph->id(node); }
    static int id(Edge edge) { return edge.id; }
    static int id(UEdge edge) { return edge.id; }

    Node nodeFromId(int id) const { return graph->nodeFromId(id); }
    static Edge edgeFromId(int id) { return Edge(id); }
    static UEdge uEdgeFromId(int id) { return UEdge(id);}

    int maxNodeId() const { return graph->maxNodeId(); };
    int maxEdgeId() const { return edges.size() - 1; }
    int maxUEdgeId() const { return edges.size() / 2 - 1; }

    Node source(Edge e) const { return edges[e.id ^ 1].target; }
    Node target(Edge e) const { return edges[e.id].target; }

    Node source(UEdge e) const { return edges[2 * e.id].target; }
    Node target(UEdge e) const { return edges[2 * e.id + 1].target; }

    typedef typename ItemSetTraits<Graph, Node>::ItemNotifier NodeNotifier;

    NodeNotifier& notifier(Node) const {
      return graph->notifier(Node());
    } 

    template <typename _Value>
    class NodeMap : public Graph::template NodeMap<_Value> {
    public:

      typedef typename _Graph::template NodeMap<_Value> Parent;

      explicit NodeMap(const SmartUEdgeSetBase<Graph>& edgeset) 
	: Parent(*edgeset.graph) { }

      NodeMap(const SmartUEdgeSetBase<Graph>& edgeset, const _Value& value)
	: Parent(*edgeset.graph, value) { }

      NodeMap& operator=(const NodeMap& cmap) {
        return operator=<NodeMap>(cmap);
      }

      template <typename CMap>
      NodeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
        return *this;
      }
    };

  };

  /// \ingroup semi_adaptors
  ///
  /// \brief Graph using a node set of another graph and an
  /// own uedge set.
  ///
  /// This structure can be used to establish another graph over a node set
  /// of an existing one. The node iterator will go through the nodes of the
  /// original graph.
  ///
  /// \param _Graph The type of the graph which shares its node set with 
  /// this class. Its interface must conform to the \ref concepts::Graph
  /// "Graph" concept.
  ///
  /// In the edge extension and removing it conforms to the 
  /// \ref concepts::UGraph "UGraph" concept.
  template <typename _Graph>
  class SmartUEdgeSet : public UEdgeSetExtender<SmartUEdgeSetBase<_Graph> > {

  public:

    typedef UEdgeSetExtender<SmartUEdgeSetBase<_Graph> > Parent;
    
    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;
    
    typedef _Graph Graph;

  protected:

    typedef typename Parent::NodesImplBase NodesImplBase;

    void eraseNode(const Node& node) {
      if (typename Parent::IncEdgeIt(*this, node) == INVALID) {
        return;
      }
      throw typename NodesImplBase::Notifier::ImmediateDetach();
    }
    
    void clearNodes() {
      Parent::clear();
    }

    class NodesImpl : public NodesImplBase {
    public:
      typedef NodesImplBase Parent;
      
      NodesImpl(const Graph& graph, SmartUEdgeSet& edgeset) 
	: Parent(graph), _edgeset(edgeset) {}

      virtual ~NodesImpl() {}

      bool attached() const {
        return Parent::attached();
      }
      
    protected:

      virtual void erase(const Node& node) {
        try {
          _edgeset.eraseNode(node);
          Parent::erase(node);
        } catch (const typename NodesImplBase::Notifier::ImmediateDetach&) {
          Parent::clear();
          throw;
        }
      }
      virtual void erase(const std::vector<Node>& nodes) {
        try {
          for (int i = 0; i < int(nodes.size()); ++i) {
            _edgeset.eraseNode(nodes[i]);
          }
          Parent::erase(nodes);
        } catch (const typename NodesImplBase::Notifier::ImmediateDetach&) {
          Parent::clear();
          throw;
        }
      }
      virtual void clear() {
	_edgeset.clearNodes();
	Parent::clear();
      }

    private:
      SmartUEdgeSet& _edgeset;
    };

    NodesImpl nodes;
    
  public:

    /// \brief Constructor of the adaptor.
    /// 
    /// Constructor of the adaptor.
    SmartUEdgeSet(const Graph& graph) : nodes(graph, *this) {
      Parent::initalize(graph, nodes);
    }

    bool valid() const {
      return nodes.attached();
    }
    
  };

}

#endif
