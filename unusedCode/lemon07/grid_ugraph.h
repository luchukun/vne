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

#ifndef GRID_UGRAPH_H
#define GRID_UGRAPH_H

#include <iostream>
#include <lemon/bits/invalid.h>
#include <lemon/bits/utility.h>

#include <lemon/bits/base_extender.h>
#include <lemon/bits/graph_extender.h>

#include <lemon/dim2.h>

///\ingroup graphs
///\file
///\brief GridUGraph class.

namespace lemon {

  class GridUGraphBase {

  public:

    typedef GridUGraphBase UGraph;

    class Node;
    class Edge;

  public:

    GridUGraphBase() {}

  protected:

    void construct(int w, int h) {
      _height = h; _width = w;
      _nodeNum = h * w; _edgeNum = 2 * _nodeNum - w - h;
      _edgeLimit = _nodeNum - w;
    }

    Edge _down(Node n) const {
      if (n.id < _nodeNum - _width) {
	return Edge(n.id);
      } else {
	return INVALID;
      }
    }

    Edge _up(Node n) const {
      if (n.id >= _width) {
	return Edge(n.id - _width);
      } else {
	return INVALID;
      }
    }

    Edge _right(Node n) const {
      if (n.id % _width < _width - 1) {
	return _edgeLimit + n.id % _width + (n.id / _width) * (_width - 1);
      } else {
	return INVALID;
      }
    }

    Edge _left(Node n) const {
      if (n.id % _width > 0) {
	return _edgeLimit + n.id % _width + (n.id / _width) * (_width - 1) - 1;
      } else {
	return INVALID;
      }
    }

  public:

    class IndexError : public RuntimeError {
    public:
      virtual const char* what() const throw() {
        return "lemon::GridUGraph::IndexError";
      }  
    };

    
    Node operator()(int i, int j) const {
      LEMON_ASSERT(0 <= i && i < width() && 0 <= j  && 
                   j < height(), IndexError());
      return Node(i + j * _width);
    }

    int row(Node n) const {
      return n.id / _width;
    }
    
    int col(Node n) const {
      return n.id % _width;    
    }

    int width() const {
      return _width;
    }

    int height() const {
      return _height;
    }

    typedef True NodeNumTag;
    typedef True EdgeNumTag;

    int nodeNum() const { return _nodeNum; }
    int edgeNum() const { return _edgeNum; }

    int maxNodeId() const { return nodeNum() - 1; }
    int maxEdgeId() const { return edgeNum() - 1; }

    Node source(Edge e) const {
      if (e.id < _edgeLimit) {
	return e.id;
      } else {
	return (e.id - _edgeLimit) % (_width - 1) +
	  (e.id - _edgeLimit) / (_width - 1) * _width;
      }
    }

    Node target(Edge e) const {
      if (e.id < _edgeLimit) {
	return e.id + _width;
      } else {
	return (e.id - _edgeLimit) % (_width - 1) +
	  (e.id - _edgeLimit) / (_width - 1) * _width + 1;
      }
    }

    static int id(Node v) { return v.id; }
    static int id(Edge e) { return e.id; }

    static Node nodeFromId(int id) { return Node(id);}
    
    static Edge edgeFromId(int id) { return Edge(id);}

    typedef True FindEdgeTag;

    Edge findEdge(Node u, Node v, Edge prev = INVALID) const {
      if (prev != INVALID) return INVALID;
      if (v.id - u.id == _width) return Edge(u.id);
      if (v.id - u.id == 1 && u.id % _width < _width - 1) {
	return Edge(u.id / _width * (_width - 1) +
		    u.id % _width + _edgeLimit);
      }
      return INVALID;
    }
    
      
    class Node {
      friend class GridUGraphBase;

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
      friend class GridUGraphBase;
      
    protected:
      int id; 

      Edge(int _id) : id(_id) {}

    public:
      Edge() { }
      Edge (Invalid) { id = -1; }
      bool operator==(const Edge edge) const {return id == edge.id;}
      bool operator!=(const Edge edge) const {return id != edge.id;}
      bool operator<(const Edge edge) const {return id < edge.id;}
    };

    void first(Node& node) const {
      node.id = nodeNum() - 1;
    }

    static void next(Node& node) {
      --node.id;
    }

    void first(Edge& edge) const {
      edge.id = edgeNum() - 1;
    }

    static void next(Edge& edge) {
      --edge.id;
    }

    void firstOut(Edge& edge, const Node& node) const {
      if (node.id < _nodeNum - _width) {
	edge.id = node.id;
      } else if (node.id % _width < _width - 1) {
	edge.id = _edgeLimit + node.id % _width +
	  (node.id / _width) * (_width - 1);
      } else {
	edge.id = -1;
      }
    }

    void nextOut(Edge& edge) const {
      if (edge.id >= _edgeLimit) {
	edge.id = -1;
      } else if (edge.id % _width < _width - 1) {
	edge.id = _edgeLimit + edge.id % _width +
	  (edge.id / _width) * (_width - 1);
      } else {
	edge.id = -1;
      }
    }

    void firstIn(Edge& edge, const Node& node) const {
      if (node.id >= _width) {
	edge.id = node.id - _width;
      } else if (node.id % _width > 0) {
	edge.id = _edgeLimit + node.id % _width +
	  (node.id / _width) * (_width - 1) - 1;
      } else {
	edge.id = -1;
      }
    }
    
    void nextIn(Edge& edge) const {
      if (edge.id >= _edgeLimit) {
	edge.id = -1;
      } else if (edge.id % _width > 0) {
	edge.id = _edgeLimit + edge.id % _width +
	  (edge.id / _width + 1) * (_width - 1) - 1;
      } else {
	edge.id = -1;
      }
    }

  private:
    int _width, _height;
    int _nodeNum, _edgeNum;
    int _edgeLimit;
  };


  typedef UGraphExtender<UndirGraphExtender<GridUGraphBase> > 
  ExtendedGridUGraphBase;

  /// \ingroup graphs
  ///
  /// \brief Grid graph class
  ///
  /// This class implements a special graph type. The nodes of the
  /// graph can be indiced by two integer \c (i,j) value where \c i
  /// is in the \c [0,width) range and j is in the [0, height) range.
  /// Two nodes are connected in the graph if the indices differ only
  /// on one position and only one is the difference. 
  ///
  /// \image html grid_ugraph.png
  /// \image latex grid_ugraph.eps "Grid graph" width=\textwidth
  ///
  /// The graph can be indiced in the following way:
  ///\code
  /// GridUGraph graph(w, h);
  /// GridUGraph::NodeMap<int> val(graph); 
  /// for (int i = 0; i < graph.width(); ++i) {
  ///   for (int j = 0; j < graph.height(); ++j) {
  ///     val[graph(i, j)] = i + j;
  ///   }
  /// }
  ///\endcode
  ///
  /// The graph type is fully conform to the \ref concepts::UGraph
  /// "Undirected Graph" concept,  and it also has an
  ///important extra feature that
  ///its maps are real \ref concepts::ReferenceMap "reference map"s.
  ///
  ///
  /// \author Balazs Dezso
  class GridUGraph : public ExtendedGridUGraphBase {
  public:

    typedef ExtendedGridUGraphBase Parent;

    /// \brief Map to get the indices of the nodes as dim2::Point<int>.
    ///
    /// Map to get the indices of the nodes as dim2::Point<int>.
    class IndexMap {
    public:
      /// \brief The key type of the map
      typedef GridUGraph::Node Key;
      /// \brief The value type of the map
      typedef dim2::Point<int> Value;

      /// \brief Constructor
      ///
      /// Constructor
      IndexMap(const GridUGraph& _graph) : graph(_graph) {}

      /// \brief The subscript operator
      ///
      /// The subscript operator.
      Value operator[](Key key) const {
	return dim2::Point<int>(graph.row(key), graph.col(key));
      }

    private:
      const GridUGraph& graph;
    };

    /// \brief Map to get the row of the nodes.
    ///
    /// Map to get the row of the nodes.
    class RowMap {
    public:
      /// \brief The key type of the map
      typedef GridUGraph::Node Key;
      /// \brief The value type of the map
      typedef int Value;

      /// \brief Constructor
      ///
      /// Constructor
      RowMap(const GridUGraph& _graph) : graph(_graph) {}

      /// \brief The subscript operator
      ///
      /// The subscript operator.
      Value operator[](Key key) const {
	return graph.row(key);
      }

    private:
      const GridUGraph& graph;
    };

    /// \brief Map to get the column of the nodes.
    ///
    /// Map to get the column of the nodes.
    class ColMap {
    public:
      /// \brief The key type of the map
      typedef GridUGraph::Node Key;
      /// \brief The value type of the map
      typedef int Value;

      /// \brief Constructor
      ///
      /// Constructor
      ColMap(const GridUGraph& _graph) : graph(_graph) {}

      /// \brief The subscript operator
      ///
      /// The subscript operator.
      Value operator[](Key key) const {
	return graph.col(key);
      }

    private:
      const GridUGraph& graph;
    };

    /// \brief Constructor
    ///
    /// 
    GridUGraph(int n, int m) { construct(n, m); }

    /// \brief Resize the graph
    ///
    void resize(int n, int m) {
      Parent::notifier(Edge()).clear();
      Parent::notifier(UEdge()).clear();
      Parent::notifier(Node()).clear();
      construct(n, m);
      Parent::notifier(Node()).build();
      Parent::notifier(UEdge()).build();
      Parent::notifier(Edge()).build();
    }
    
    /// \brief The node on the given position.
    /// 
    /// Gives back the node on the given position.
    Node operator()(int i, int j) const {
      return Parent::operator()(i, j);
    }

    /// \brief Gives back the row index of the node.
    ///
    /// Gives back the row index of the node.
    int row(Node n) const {
      return Parent::row(n);
    }
    
    /// \brief Gives back the coloumn index of the node.
    ///
    /// Gives back the coloumn index of the node.
    int col(Node n) const {
      return Parent::col(n);
    }

    /// \brief Gives back the width of the graph.
    ///
    /// Gives back the width of the graph.
    int width() const {
      return Parent::width();
    }

    /// \brief Gives back the height of the graph.
    ///
    /// Gives back the height of the graph.
    int height() const {
      return Parent::height();
    }

    /// \brief Gives back the edge goes down from the node.
    ///
    /// Gives back the edge goes down from the node. If there is not
    /// outgoing edge then it gives back INVALID.
    Edge down(Node n) const {
      UEdge ue = _down(n);
      return ue != INVALID ? direct(ue, true) : INVALID;
    }
    
    /// \brief Gives back the edge goes up from the node.
    ///
    /// Gives back the edge goes up from the node. If there is not
    /// outgoing edge then it gives back INVALID.
    Edge up(Node n) const {
      UEdge ue = _up(n);
      return ue != INVALID ? direct(ue, false) : INVALID;
    }

    /// \brief Gives back the edge goes right from the node.
    ///
    /// Gives back the edge goes right from the node. If there is not
    /// outgoing edge then it gives back INVALID.
    Edge right(Node n) const {
      UEdge ue = _right(n);
      return ue != INVALID ? direct(ue, true) : INVALID;
    }

    /// \brief Gives back the edge goes left from the node.
    ///
    /// Gives back the edge goes left from the node. If there is not
    /// outgoing edge then it gives back INVALID.
    Edge left(Node n) const {
      UEdge ue = _left(n);
      return ue != INVALID ? direct(ue, false) : INVALID;
    }
    
  };

  /// \brief Index map of the grid graph
  ///
  /// Just returns an IndexMap for the grid graph.
  inline GridUGraph::IndexMap indexMap(const GridUGraph& graph) {
    return GridUGraph::IndexMap(graph);
  }

  /// \brief Row map of the grid graph
  ///
  /// Just returns a RowMap for the grid graph.
  inline GridUGraph::RowMap rowMap(const GridUGraph& graph) {
    return GridUGraph::RowMap(graph);
  }

  /// \brief Column map of the grid graph
  ///
  /// Just returns a ColMap for the grid graph.
  inline GridUGraph::ColMap colMap(const GridUGraph& graph) {
    return GridUGraph::ColMap(graph);
  }
}
#endif
