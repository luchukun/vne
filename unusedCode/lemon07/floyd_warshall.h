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

#ifndef LEMON_FLOYD_WARSHALL_H
#define LEMON_FLOYD_WARSHALL_H

///\ingroup shortest_path
/// \file
/// \brief FloydWarshall algorithm.
///

#include <lemon/list_graph.h>
#include <lemon/graph_utils.h>
#include <lemon/bits/path_dump.h>
#include <lemon/bits/invalid.h>
#include <lemon/error.h>
#include <lemon/matrix_maps.h>
#include <lemon/maps.h>

#include <limits>

namespace lemon {

  /// \brief Default OperationTraits for the FloydWarshall algorithm class.
  ///  
  /// It defines all computational operations and constants which are
  /// used in the Floyd-Warshall algorithm. The default implementation
  /// is based on the numeric_limits class. If the numeric type does not
  /// have infinity value then the maximum value is used as extremal
  /// infinity value.
  template <
    typename Value, 
    bool has_infinity = std::numeric_limits<Value>::has_infinity>
  struct FloydWarshallDefaultOperationTraits {
    /// \brief Gives back the zero value of the type.
    static Value zero() {
      return static_cast<Value>(0);
    }
    /// \brief Gives back the positive infinity value of the type.
    static Value infinity() {
      return std::numeric_limits<Value>::infinity();
    }
    /// \brief Gives back the sum of the given two elements.
    static Value plus(const Value& left, const Value& right) {
      return left + right;
    }
    /// \brief Gives back true only if the first value less than the second.
    static bool less(const Value& left, const Value& right) {
      return left < right;
    }
  };

  template <typename Value>
  struct FloydWarshallDefaultOperationTraits<Value, false> {
    static Value zero() {
      return static_cast<Value>(0);
    }
    static Value infinity() {
      return std::numeric_limits<Value>::max();
    }
    static Value plus(const Value& left, const Value& right) {
      if (left == infinity() || right == infinity()) return infinity();
      return left + right;
    }
    static bool less(const Value& left, const Value& right) {
      return left < right;
    }
  };
  
  /// \brief Default traits class of FloydWarshall class.
  ///
  /// Default traits class of FloydWarshall class.
  /// \param _Graph Graph type.
  /// \param _LegthMap Type of length map.
  template<class _Graph, class _LengthMap>
  struct FloydWarshallDefaultTraits {
    /// The graph type the algorithm runs on. 
    typedef _Graph Graph;

    /// \brief The type of the map that stores the edge lengths.
    ///
    /// The type of the map that stores the edge lengths.
    /// It must meet the \ref concepts::ReadMap "ReadMap" concept.
    typedef _LengthMap LengthMap;

    // The type of the length of the edges.
    typedef typename _LengthMap::Value Value;

    /// \brief Operation traits for floyd-warshall algorithm.
    ///
    /// It defines the infinity type on the given Value type
    /// and the used operation.
    /// \see FloydWarshallDefaultOperationTraits
    typedef FloydWarshallDefaultOperationTraits<Value> OperationTraits;
 
    /// \brief The type of the matrix map that stores the last edges of the 
    /// shortest paths.
    /// 
    /// The type of the map that stores the last edges of the shortest paths.
    /// It must be a matrix map with \c Graph::Edge value type.
    ///
    typedef DynamicMatrixMap<Graph, typename Graph::Node, 
			     typename Graph::Edge> PredMap;

    /// \brief Instantiates a PredMap.
    /// 
    /// This function instantiates a \ref PredMap. 
    /// \param graph is the graph,
    /// to which we would like to define the PredMap.
    /// \todo The graph alone may be insufficient for the initialization
    static PredMap *createPredMap(const _Graph& graph) {
      return new PredMap(graph);
    }

    /// \brief The type of the map that stores the dists of the nodes.
    ///
    /// The type of the map that stores the dists of the nodes.
    /// It must meet the \ref concepts::WriteMatrixMap "WriteMatrixMap" concept.
    ///
    typedef DynamicMatrixMap<Graph, typename Graph::Node, Value> DistMap;

    /// \brief Instantiates a DistMap.
    ///
    /// This function instantiates a \ref DistMap. 
    /// \param graph is the graph, to which we would like to define the 
    /// \ref DistMap
    static DistMap *createDistMap(const _Graph& graph) {
      return new DistMap(graph);
    }

  };
  
  /// \brief %FloydWarshall algorithm class.
  ///
  /// \ingroup shortest_path
  /// This class provides an efficient implementation of \c Floyd-Warshall 
  /// algorithm. The edge lengths are passed to the algorithm using a
  /// \ref concepts::ReadMap "ReadMap", so it is easy to change it to any 
  /// kind of length.
  ///
  /// The algorithm solves the shortest path problem for each pair
  /// of node when the edges can have negative length but the graph should
  /// not contain cycles with negative sum of length. If we can assume
  /// that all edge is non-negative in the graph then the dijkstra algorithm
  /// should be used from each node rather and if the graph is sparse and
  /// there are negative circles then the johnson algorithm.
  ///
  /// The complexity of this algorithm is \f$ O(n^3+e) \f$.
  ///
  /// The type of the length is determined by the
  /// \ref concepts::ReadMap::Value "Value" of the length map.
  ///
  /// \param _Graph The graph type the algorithm runs on. The default value
  /// is \ref ListGraph. The value of _Graph is not used directly by
  /// FloydWarshall, it is only passed to \ref FloydWarshallDefaultTraits.
  /// \param _LengthMap This read-only EdgeMap determines the lengths of the
  /// edges. It is read once for each edge, so the map may involve in
  /// relatively time consuming process to compute the edge length if
  /// it is necessary. The default map type is \ref
  /// concepts::Graph::EdgeMap "Graph::EdgeMap<int>".  The value
  /// of _LengthMap is not used directly by FloydWarshall, it is only passed 
  /// to \ref FloydWarshallDefaultTraits.  \param _Traits Traits class to set
  /// various data types used by the algorithm.  The default traits
  /// class is \ref FloydWarshallDefaultTraits
  /// "FloydWarshallDefaultTraits<_Graph,_LengthMap>".  See \ref
  /// FloydWarshallDefaultTraits for the documentation of a FloydWarshall 
  /// traits class.
  ///
  /// \author Balazs Dezso
  /// \todo A function type interface would be nice.
  /// \todo Implement \c nextNode() and \c nextEdge()
#ifdef DOXYGEN
  template <typename _Graph, typename _LengthMap, typename _Traits >
#else
  template <typename _Graph=ListGraph,
	    typename _LengthMap=typename _Graph::template EdgeMap<int>,
	    typename _Traits=FloydWarshallDefaultTraits<_Graph,_LengthMap> >
#endif
  class FloydWarshall {
  public:
    
    /// \brief \ref Exception for uninitialized parameters.
    ///
    /// This error represents problems in the initialization
    /// of the parameters of the algorithms.

    class UninitializedParameter : public lemon::UninitializedParameter {
    public:
      virtual const char* what() const throw() {
	return "lemon::FloydWarshall::UninitializedParameter";
      }
    };

    typedef _Traits Traits;
    ///The type of the underlying graph.
    typedef typename _Traits::Graph Graph;

    typedef typename Graph::Node Node;
    typedef typename Graph::NodeIt NodeIt;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::EdgeIt EdgeIt;
    
    /// \brief The type of the length of the edges.
    typedef typename _Traits::LengthMap::Value Value;
    /// \brief The type of the map that stores the edge lengths.
    typedef typename _Traits::LengthMap LengthMap;
    /// \brief The type of the map that stores the last
    /// edges of the shortest paths. The type of the PredMap
    /// is a matrix map for Edges
    typedef typename _Traits::PredMap PredMap;
    /// \brief The type of the map that stores the dists of the nodes.
    /// The type of the DistMap is a matrix map for Values
    ///
    /// \todo It should rather be
    /// called \c DistMatrix
    typedef typename _Traits::DistMap DistMap;
    /// \brief The operation traits.
    typedef typename _Traits::OperationTraits OperationTraits;
  private:
    /// Pointer to the underlying graph.
    const Graph *graph;
    /// Pointer to the length map
    const LengthMap *length;
    ///Pointer to the map of predecessors edges.
    PredMap *_pred;
    ///Indicates if \ref _pred is locally allocated (\c true) or not.
    bool local_pred;
    ///Pointer to the map of distances.
    DistMap *_dist;
    ///Indicates if \ref _dist is locally allocated (\c true) or not.
    bool local_dist;

    /// Creates the maps if necessary.
    void create_maps() {
      if(!_pred) {
	local_pred = true;
	_pred = Traits::createPredMap(*graph);
      }
      if(!_dist) {
	local_dist = true;
	_dist = Traits::createDistMap(*graph);
      }
    }
    
  public :
 
    /// \name Named template parameters

    ///@{

    template <class T>
    struct DefPredMapTraits : public Traits {
      typedef T PredMap;
      static PredMap *createPredMap(const Graph& graph) {
	throw UninitializedParameter();
      }
    };

    /// \brief \ref named-templ-param "Named parameter" for setting PredMap 
    /// type
    /// \ref named-templ-param "Named parameter" for setting PredMap type
    ///
    template <class T>
    struct DefPredMap 
      : public FloydWarshall< Graph, LengthMap, DefPredMapTraits<T> > {
      typedef FloydWarshall< Graph, LengthMap, DefPredMapTraits<T> > Create;
    };
    
    template <class T>
    struct DefDistMapTraits : public Traits {
      typedef T DistMap;
      static DistMap *createDistMap(const Graph& graph) {
	throw UninitializedParameter();
      }
    };
    /// \brief \ref named-templ-param "Named parameter" for setting DistMap 
    /// type
    ///
    /// \ref named-templ-param "Named parameter" for setting DistMap type
    ///
    template <class T>
    struct DefDistMap 
      : public FloydWarshall< Graph, LengthMap, DefDistMapTraits<T> > {
      typedef FloydWarshall< Graph, LengthMap, DefDistMapTraits<T> > Create;
    };
    
    template <class T>
    struct DefOperationTraitsTraits : public Traits {
      typedef T OperationTraits;
    };
    
    /// \brief \ref named-templ-param "Named parameter" for setting 
    /// OperationTraits type
    ///
    /// \ref named-templ-param "Named parameter" for setting PredMap type
    template <class T>
    struct DefOperationTraits
      : public FloydWarshall< Graph, LengthMap, DefOperationTraitsTraits<T> > {
      typedef FloydWarshall< Graph, LengthMap, DefOperationTraitsTraits<T> >
      Create;
    };
    
    ///@}

  protected:

    FloydWarshall() {}

  public:      

    typedef FloydWarshall Create;
    
    /// \brief Constructor.
    ///
    /// \param _graph the graph the algorithm will run on.
    /// \param _length the length map used by the algorithm.
    FloydWarshall(const Graph& _graph, const LengthMap& _length) :
      graph(&_graph), length(&_length),
      _pred(0), local_pred(false),
      _dist(0), local_dist(false) {}
    
    ///Destructor.
    ~FloydWarshall() {
      if(local_pred) delete _pred;
      if(local_dist) delete _dist;
    }

    /// \brief Sets the length map.
    ///
    /// Sets the length map.
    /// \return \c (*this)
    FloydWarshall &lengthMap(const LengthMap &m) {
      length = &m;
      return *this;
    }

    /// \brief Sets the map storing the predecessor edges.
    ///
    /// Sets the map storing the predecessor edges.
    /// If you don't use this function before calling \ref run(),
    /// it will allocate one. The destuctor deallocates this
    /// automatically allocated map, of course.
    /// \return \c (*this)
    FloydWarshall &predMap(PredMap &m) {
      if(local_pred) {
	delete _pred;
	local_pred=false;
      }
      _pred = &m;
      return *this;
    }

    /// \brief Sets the map storing the distances calculated by the algorithm.
    ///
    /// Sets the map storing the distances calculated by the algorithm.
    /// If you don't use this function before calling \ref run(),
    /// it will allocate one. The destuctor deallocates this
    /// automatically allocated map, of course.
    /// \return \c (*this)
    FloydWarshall &distMap(DistMap &m) {
      if(local_dist) {
	delete _dist;
	local_dist=false;
      }
      _dist = &m;
      return *this;
    }

    ///\name Execution control
    /// The simplest way to execute the algorithm is to use
    /// one of the member functions called \c run(...).
    /// \n
    /// If you need more control on the execution,
    /// Finally \ref start() will perform the actual path
    /// computation.

    ///@{

    /// \brief Initializes the internal data structures.
    /// 
    /// Initializes the internal data structures.
    void init() {
      create_maps();
      for (NodeIt it(*graph); it != INVALID; ++it) {
	for (NodeIt jt(*graph); jt != INVALID; ++jt) {
	  _pred->set(it, jt, INVALID);
	  _dist->set(it, jt, OperationTraits::infinity());
	}
	_dist->set(it, it, OperationTraits::zero());
      }
      for (EdgeIt it(*graph); it != INVALID; ++it) {
	Node source = graph->source(it);
	Node target = graph->target(it);	
	if (OperationTraits::less((*length)[it], (*_dist)(source, target))) {
	  _dist->set(source, target, (*length)[it]);
	  _pred->set(source, target, it);
	}
      }
    }
    
    /// \brief Executes the algorithm.
    ///
    /// This method runs the %FloydWarshall algorithm in order to compute 
    /// the shortest path to each node pairs. The algorithm 
    /// computes 
    /// - The shortest path tree for each node.
    /// - The distance between each node pairs.
    void start() {
      for (NodeIt kt(*graph); kt != INVALID; ++kt) {
	for (NodeIt it(*graph); it != INVALID; ++it) {
	  for (NodeIt jt(*graph); jt != INVALID; ++jt) {
	    Value relaxed = OperationTraits::plus((*_dist)(it, kt),
						  (*_dist)(kt, jt));
	    if (OperationTraits::less(relaxed, (*_dist)(it, jt))) {
	      _dist->set(it, jt, relaxed);
	      _pred->set(it, jt, (*_pred)(kt, jt));
	    }
	  }
	}
      }
    }

    /// \brief Executes the algorithm and checks the negative cycles.
    ///
    /// This method runs the %FloydWarshall algorithm in order to compute 
    /// the shortest path to each node pairs. If there is a negative cycle 
    /// in the graph it gives back false. 
    /// The algorithm computes 
    /// - The shortest path tree for each node.
    /// - The distance between each node pairs.
    bool checkedStart() {
      start();
      for (NodeIt it(*graph); it != INVALID; ++it) {
	if (OperationTraits::less((*dist)(it, it), OperationTraits::zero())) {
	  return false;
	}
      }
      return true;
    }
    
    /// \brief Runs %FloydWarshall algorithm.
    ///    
    /// This method runs the %FloydWarshall algorithm from a each node
    /// in order to compute the shortest path to each node pairs. 
    /// The algorithm computes
    /// - The shortest path tree for each node.
    /// - The distance between each node pairs.
    ///
    /// \note d.run(s) is just a shortcut of the following code.
    ///\code
    ///  d.init();
    ///  d.start();
    ///\endcode
    void run() {
      init();
      start();
    }
    
    ///@}

    /// \name Query Functions
    /// The result of the %FloydWarshall algorithm can be obtained using these
    /// functions.\n
    /// Before the use of these functions,
    /// either run() or start() must be called.
    
    ///@{

    typedef PredMatrixMapPath<Graph, PredMap> Path;

    ///Gives back the shortest path.
    
    ///Gives back the shortest path.
    ///\pre The \c t should be reachable from the \c t.
    Path path(Node s, Node t) 
    {
      return Path(*graph, *_pred, s, t);
    }
	  
    /// \brief The distance between two nodes.
    ///
    /// Returns the distance between two nodes.
    /// \pre \ref run() must be called before using this function.
    /// \warning If node \c v in unreachable from the root the return value
    /// of this funcion is undefined.
    Value dist(Node source, Node target) const { 
      return (*_dist)(source, target); 
    }

    /// \brief Returns the 'previous edge' of the shortest path tree.
    ///
    /// For the node \c node it returns the 'previous edge' of the shortest 
    /// path tree to direction of the node \c root 
    /// i.e. it returns the last edge of a shortest path from the node \c root 
    /// to \c node. It is \ref INVALID if \c node is unreachable from the root
    /// or if \c node=root. The shortest path tree used here is equal to the 
    /// shortest path tree used in \ref predNode(). 
    /// \pre \ref run() must be called before using this function.
    Edge predEdge(Node root, Node node) const { 
      return (*_pred)(root, node); 
    }

    /// \brief Returns the 'previous node' of the shortest path tree.
    ///
    /// For a node \c node it returns the 'previous node' of the shortest path 
    /// tree to direction of the node \c root, i.e. it returns the last but 
    /// one node from a shortest path from the \c root to \c node. It is 
    /// INVALID if \c node is unreachable from the root or if \c node=root. 
    /// The shortest path tree used here is equal to the 
    /// shortest path tree used in \ref predEdge().  
    /// \pre \ref run() must be called before using this function.
    Node predNode(Node root, Node node) const { 
      return (*_pred)(root, node) == INVALID ? 
      INVALID : graph->source((*_pred)(root, node)); 
    }
    
    /// \brief Returns a reference to the matrix node map of distances.
    ///
    /// Returns a reference to the matrix node map of distances. 
    ///
    /// \pre \ref run() must be called before using this function.
    const DistMap &distMap() const { return *_dist;}
 
    /// \brief Returns a reference to the shortest path tree map.
    ///
    /// Returns a reference to the matrix node map of the edges of the
    /// shortest path tree.
    /// \pre \ref run() must be called before using this function.
    const PredMap &predMap() const { return *_pred;}
 
    /// \brief Checks if a node is reachable from the root.
    ///
    /// Returns \c true if \c v is reachable from the root.
    /// \pre \ref run() must be called before using this function.
    ///
    bool connected(Node source, Node target) { 
      return (*_dist)(source, target) != OperationTraits::infinity(); 
    }
    
    ///@}
  };
 
} //END OF NAMESPACE LEMON

#endif

