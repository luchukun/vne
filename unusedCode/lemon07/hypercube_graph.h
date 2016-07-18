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

#ifndef HYPERCUBE_GRAPH_H
#define HYPERCUBE_GRAPH_H

#include <iostream>
#include <vector>
#include <lemon/bits/invalid.h>
#include <lemon/bits/utility.h>
#include <lemon/error.h>

#include <lemon/bits/base_extender.h>
#include <lemon/bits/graph_extender.h>

///\ingroup graphs
///\file
///\brief HyperCubeGraph class.

namespace lemon {

  class HyperCubeGraphBase {

  public:

    typedef HyperCubeGraphBase Graph;

    class Node;
    class Edge;

  public:

    HyperCubeGraphBase() {}

  protected:

    void construct(int dim) {
      _dim = dim;
      _nodeNum = 1 << dim;
    }

  public:
    

    typedef True NodeNumTag;
    typedef True EdgeNumTag;

    int nodeNum() const { return _nodeNum; }
    int edgeNum() const { return _nodeNum * _dim; }

    int maxNodeId() const { return nodeNum() - 1; }
    int maxEdgeId() const { return edgeNum() - 1; }

    Node source(Edge e) const {
      return e.id / _dim;
    }

    Node target(Edge e) const {
      return (e.id / _dim) ^ ( 1 << (e.id % _dim));
    }

    static int id(Node v) { return v.id; }
    static int id(Edge e) { return e.id; }

    static Node nodeFromId(int id) { return Node(id);}
    
    static Edge edgeFromId(int id) { return Edge(id);}

    class Node {
      friend class HyperCubeGraphBase;

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
      friend class HyperCubeGraphBase;
      
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
      edge.id = node.id * _dim;
    }

    void nextOut(Edge& edge) const {
      ++edge.id;
      if (edge.id % _dim == 0) edge.id = -1;
    }

    void firstIn(Edge& edge, const Node& node) const {
      edge.id = (node.id ^ 1) * _dim;
    }
    
    void nextIn(Edge& edge) const {
      int cnt = edge.id % _dim;
      if ((cnt + 1) % _dim == 0) {
	edge.id = -1;
      } else {
	edge.id = ((edge.id / _dim) ^ ((1 << cnt) * 3)) * _dim + cnt + 1; 
      }
    }

    int dimension() const {
      return _dim;
    }

    bool projection(Node node, int n) const {
      return static_cast<bool>(node.id & (1 << n));
    }

    int dimension(Edge edge) const {
      return edge.id % _dim;
    }

    int index(Node node) const {
      return node.id;
    }

    Node operator()(int ix) const {
      return Node(ix);
    }
    
  private:
    int _dim, _nodeNum;
  };


  typedef GraphExtender<HyperCubeGraphBase> ExtendedHyperCubeGraphBase;

  /// \ingroup graphs
  ///
  /// \brief HyperCube graph class
  ///
  /// This class implements a special graph type. The nodes of the
  /// graph can be indiced with integers with at most \c dim binary length.
  /// Two nodes are connected in the graph if the indices differ only
  /// on one position in the binary form. 
  ///
  /// \note The type of the \c ids is chosen to \c int because efficiency
  /// reasons. This way the maximal dimension of this implementation
  /// is 26. 
  ///
  /// The graph type is fully conform to the \ref concepts::Graph
  /// concept but it does not conform to the \ref concepts::UGraph.
  ///
  /// \author Balazs Dezso
  class HyperCubeGraph : public ExtendedHyperCubeGraphBase {
  public:

    typedef ExtendedHyperCubeGraphBase Parent;

    /// \brief Construct a graph with \c dim dimension.
    ///
    /// Construct a graph with \c dim dimension.
    HyperCubeGraph(int dim) { construct(dim); }

    /// \brief Gives back the number of the dimensions.
    ///
    /// Gives back the number of the dimensions.
    int dimension() const {
      return Parent::dimension();
    }

    /// \brief Returns true if the n'th bit of the node is one.
    ///
    /// Returns true if the n'th bit of the node is one. 
    bool projection(Node node, int n) const {
      return Parent::projection(node, n);
    }

    /// \brief The dimension id of the edge.
    ///
    /// It returns the dimension id of the edge. It can
    /// be in the \f$ \{0, 1, \dots, dim-1\} \f$ intervall.
    int dimension(Edge edge) const {
      return Parent::dimension(edge);
    }

    /// \brief Gives back the index of the node.
    ///
    /// Gives back the index of the node. The lower bits of the
    /// integer describes the node.
    int index(Node node) const {
      return Parent::index(node);
    }

    /// \brief Gives back the node by its index.
    ///
    /// Gives back the node by its index.
    Node operator()(int ix) const {
      return Parent::operator()(ix);
    }

    /// \brief Number of nodes.
    int nodeNum() const { return Parent::nodeNum(); }
    /// \brief Number of edges.
    int edgeNum() const { return Parent::edgeNum(); }

    /// \brief Linear combination map.
    ///
    /// It makes possible to give back a linear combination
    /// for each node. This function works like the \c std::accumulate
    /// so it accumulates the \c bf binary function with the \c fv
    /// first value. The map accumulates only on that dimensions where
    /// the node's index is one. The accumulated values should be
    /// given by the \c begin and \c end iterators and this range's length
    /// should be the dimension number of the graph.
    /// 
    ///\code
    /// const int DIM = 3;
    /// HyperCubeGraph graph(DIM);
    /// dim2::Point<double> base[DIM];
    /// for (int k = 0; k < DIM; ++k) {
    ///   base[k].x = rnd();
    ///   base[k].y = rnd();
    /// } 
    /// HyperCubeGraph::HyperMap<dim2::Point<double> > 
    ///   pos(graph, base, base + DIM, dim2::Point<double>(0.0, 0.0));
    ///\endcode
    ///
    /// \see HyperCubeGraph
    template <typename T, typename BF = std::plus<T> >
    class HyperMap {
    public:

      typedef Node Key;
      typedef T Value;
    
      
      /// \brief Constructor for HyperMap. 
      ///
      /// Construct a HyperMap for the given graph. The accumulated values 
      /// should be given by the \c begin and \c end iterators and this 
      /// range's length should be the dimension number of the graph.
      ///
      /// This function accumulates the \c bf binary function with 
      /// the \c fv first value. The map accumulates only on that dimensions 
      /// where the node's index is one.           
      template <typename It>
      HyperMap(const Graph& graph, It begin, It end, 
		   T fv = 0.0, const BF& bf = BF()) 
	: _graph(graph), _values(begin, end), _first_value(fv), _bin_func(bf) {
	LEMON_ASSERT(_values.size() == graph.dimension(), 
		     "Wrong size of dimension");
      }

      /// \brief Gives back the partial accumulated value.
      ///
      /// Gives back the partial accumulated value.
      Value operator[](Key k) const {
	Value val = _first_value;
	int id = _graph.index(k); 
	int n = 0;
	while (id != 0) {
	  if (id & 1) {
	    val = _bin_func(val, _values[n]);
	  }
	  id >>= 1;
	  ++n;
	}
	return val;
      }
      
    private:
      const Graph& _graph;
      std::vector<T> _values;
      T _first_value;
      BF _bin_func;
    };    
  };
}
#endif

