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

#ifndef LEMON_STATIC_GRAPH_H
#define LEMON_STATIC_GRAPH_H

#include <lemon/bits/graph_extender.h>
#include <lemon/graph_utils.h>

namespace lemon {

  class StaticGraphBase {
  public:

    StaticGraphBase() 
      : node_num(-1), edge_num(0), 
        node_first_out(0), node_first_in(0),
        edge_source(0), edge_target(0), 
        edge_next_in(0), edge_next_out(0) {}
    
    ~StaticGraphBase() {
      if (node_num != -1) {
        delete[] node_first_out;
        delete[] node_first_in;
        delete[] edge_source;
        delete[] edge_target;
        delete[] edge_next_out;
        delete[] edge_next_in;
      }
    }

    class Node {
      friend class StaticGraphBase;
    protected:
      int id;
      Node(int _id) : id(_id) {}
    public:
      Node() {}
      Node (Invalid) : id(-1) {}
      bool operator==(const Node& node) const {return id == node.id;}
      bool operator!=(const Node& node) const {return id != node.id;}
      bool operator<(const Node& node) const {return id < node.id;}
    };

    class Edge {
      friend class StaticGraphBase;      
    protected:
      int id;
      Edge(int _id) : id(_id) {}
    public:
      Edge() { }
      Edge (Invalid) { id = -1; }
      bool operator==(const Edge& edge) const {return id == edge.id;}
      bool operator!=(const Edge& edge) const {return id != edge.id;}
      bool operator<(const Edge& edge) const {return id < edge.id;}
    };

    Node source(const Edge& e) const { return Node(edge_source[e.id]); }
    Node target(const Edge& e) const { return Node(edge_target[e.id]); }

    void first(Node& n) const { n.id = node_num - 1; }
    void next(Node& n) const { --n.id; }

    void first(Edge& n) const { n.id = edge_num - 1; }
    void next(Edge& n) const { --n.id; }

    void firstOut(Edge& e, const Node& n) const { 
      e.id = node_first_out[n.id] != node_first_out[n.id + 1] ? 
        node_first_out[n.id] : -1;
    }
    void nextOut(Edge& e) const { e.id = edge_next_out[e.id]; }

    void firstIn(Edge& e, const Node& n) const { e.id = node_first_in[n.id]; }
    void nextIn(Edge& e) const { e.id = edge_next_in[e.id]; }
    

    int id(const Node& n) const { return n.id; }
    Node nodeFromId(int id) const { return Node(id); }
    int maxNodeId() const { return node_num - 1; }

    int id(const Edge& e) const { return e.id; }
    Edge edgeFromId(int id) const { return Edge(id); }
    int maxEdgeId() const { return edge_num - 1; }

    typedef True NodeNumTag;
    int nodeNum() const { return node_num; }
    typedef True EdgeNumTag;
    int edgeNum() const { return edge_num; }

  private:

    template <typename Graph, typename NodeRefMap>
    class EdgeLess {
    public:
      typedef typename Graph::Edge Edge;

      EdgeLess(const Graph &_graph, const NodeRefMap& _nodeRef) 
        : graph(_graph), nodeRef(_nodeRef) {}
      
      bool operator()(const Edge& left, const Edge& right) const {
	return nodeRef[graph.target(left)] < nodeRef[graph.target(right)];
      }
    private:
      const Graph& graph;
      const NodeRefMap& nodeRef;
    };
    
  public:

    typedef True BuildTag;
    
    template <typename Graph, typename NodeRefMap, typename EdgeRefMap>
    void build(const Graph& graph, NodeRefMap& nodeRef, EdgeRefMap& edgeRef) {

      if (node_num != -1) {
        delete[] node_first_out;
        delete[] node_first_in;
        delete[] edge_source;
        delete[] edge_target;
        delete[] edge_next_out;
        delete[] edge_next_in;
      }

      typedef typename Graph::Node GNode;
      typedef typename Graph::Edge GEdge;

      node_num = countNodes(graph);
      edge_num = countEdges(graph);

      node_first_out = new int[node_num + 1];
      node_first_in = new int[node_num];

      edge_source = new int[edge_num];
      edge_target = new int[edge_num];
      edge_next_out = new int[edge_num];
      edge_next_in = new int[edge_num];

      int node_index = 0;
      for (typename Graph::NodeIt n(graph); n != INVALID; ++n) {
        nodeRef[n] = Node(node_index);
        node_first_in[node_index] = -1;
        ++node_index;
      }

      EdgeLess<Graph, NodeRefMap> edgeLess(graph, nodeRef);

      int edge_index = 0;
      for (typename Graph::NodeIt n(graph); n != INVALID; ++n) {
        int source = nodeRef[n].id;
        std::vector<GEdge> edges;
        for (typename Graph::OutEdgeIt e(graph, n); e != INVALID; ++e) {
          edges.push_back(e);
        }
        if (!edges.empty()) {
          node_first_out[source] = edge_index;
          std::sort(edges.begin(), edges.end(), edgeLess);
          for (typename std::vector<GEdge>::iterator it = edges.begin();
               it != edges.end(); ++it) {
            int target = nodeRef[graph.target(*it)].id;
            edgeRef[*it] = Edge(edge_index);
            edge_source[edge_index] = source; 
            edge_target[edge_index] = target;
            edge_next_in[edge_index] = node_first_in[target];
            node_first_in[target] = edge_index;
            edge_next_out[edge_index] = edge_index + 1;
            ++edge_index;
          }
          edge_next_out[edge_index - 1] = -1;
        } else {
          node_first_out[source] = edge_index;
        }
      }
      node_first_out[node_num] = edge_num;
    }

  protected:

    void fastFirstOut(Edge& e, const Node& n) const {
      e.id = node_first_out[n.id];
    }

    static void fastNextOut(Edge& e) {
      ++e.id;
    }
    void fastLastOut(Edge& e, const Node& n) const {
      e.id = node_first_out[n.id + 1];
    }
    
  private:
    int node_num;
    int edge_num;
    int *node_first_out;
    int *node_first_in;
    int *edge_source;
    int *edge_target;
    int *edge_next_in;
    int *edge_next_out;
  };

  typedef GraphExtender<StaticGraphBase> ExtendedStaticGraphBase;


  class StaticGraph : public ExtendedStaticGraphBase {
  public:

    typedef ExtendedStaticGraphBase Parent;

  protected:

    using Parent::fastFirstOut;
    using Parent::fastNextOut;
    using Parent::fastLastOut;
    
  public:
    
    class OutEdgeIt : public Edge {
    public:

      OutEdgeIt() { }

      OutEdgeIt(Invalid i) : Edge(i) { }

      OutEdgeIt(const StaticGraph& graph, const Node& node) {
	graph.fastFirstOut(*this, node);
	graph.fastLastOut(last, node);
        if (last == *this) *this = INVALID;
      }

      OutEdgeIt(const StaticGraph& graph, const Edge& edge) : Edge(edge) {
        graph.fastLastOut(last, graph.source(edge));
      }

      OutEdgeIt& operator++() { 
        StaticGraph::fastNextOut(*this);
        if (last == *this) *this = INVALID;
        return *this; 
      }

    private:
      Edge last;
    };
    
  };

}

#endif
