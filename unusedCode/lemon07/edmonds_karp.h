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

#ifndef LEMON_EDMONDS_KARP_H
#define LEMON_EDMONDS_KARP_H

/// \file
/// \ingroup max_flow
/// \brief Implementation of the Edmonds-Karp algorithm.

#include <lemon/tolerance.h>
#include <vector>

namespace lemon {

  /// \brief Default traits class of EdmondsKarp class.
  ///
  /// Default traits class of EdmondsKarp class.
  /// \param _Graph Graph type.
  /// \param _CapacityMap Type of capacity map.
  template <typename _Graph, typename _CapacityMap>
  struct EdmondsKarpDefaultTraits {

    /// \brief The graph type the algorithm runs on. 
    typedef _Graph Graph;

    /// \brief The type of the map that stores the edge capacities.
    ///
    /// The type of the map that stores the edge capacities.
    /// It must meet the \ref concepts::ReadMap "ReadMap" concept.
    typedef _CapacityMap CapacityMap;

    /// \brief The type of the length of the edges.
    typedef typename CapacityMap::Value Value;

    /// \brief The map type that stores the flow values.
    ///
    /// The map type that stores the flow values. 
    /// It must meet the \ref concepts::ReadWriteMap "ReadWriteMap" concept.
    typedef typename Graph::template EdgeMap<Value> FlowMap;

    /// \brief Instantiates a FlowMap.
    ///
    /// This function instantiates a \ref FlowMap. 
    /// \param graph The graph, to which we would like to define the flow map.
    static FlowMap* createFlowMap(const Graph& graph) {
      return new FlowMap(graph);
    }

    /// \brief The tolerance used by the algorithm
    ///
    /// The tolerance used by the algorithm to handle inexact computation.
    typedef Tolerance<Value> Tolerance;

  };

  /// \ingroup max_flow
  ///
  /// \brief Edmonds-Karp algorithms class.
  ///
  /// This class provides an implementation of the \e Edmonds-Karp \e
  /// algorithm producing a flow of maximum value in a directed
  /// graphs. The Edmonds-Karp algorithm is slower than the Preflow
  /// algorithm but it has an advantage of the step-by-step execution
  /// control with feasible flow solutions. The \e source node, the \e
  /// target node, the \e capacity of the edges and the \e starting \e
  /// flow value of the edges should be passed to the algorithm
  /// through the constructor.
  ///
  /// The time complexity of the algorithm is \f$ O(nm^2) \f$ in
  /// worst case.  Always try the preflow algorithm instead of this if
  /// you just want to compute the optimal flow.
  ///
  /// \param _Graph The directed graph type the algorithm runs on.
  /// \param _CapacityMap The capacity map type.
  /// \param _Traits Traits class to set various data types used by
  /// the algorithm.  The default traits class is \ref
  /// EdmondsKarpDefaultTraits.  See \ref EdmondsKarpDefaultTraits for the
  /// documentation of a Edmonds-Karp traits class. 
  ///
  /// \author Balazs Dezso 
#ifdef DOXYGEN
  template <typename _Graph, typename _CapacityMap, typename _Traits>
#else 
  template <typename _Graph,
	    typename _CapacityMap = typename _Graph::template EdgeMap<int>,
            typename _Traits = EdmondsKarpDefaultTraits<_Graph, _CapacityMap> >
#endif
  class EdmondsKarp {
  public:

    typedef _Traits Traits;
    typedef typename Traits::Graph Graph;
    typedef typename Traits::CapacityMap CapacityMap;
    typedef typename Traits::Value Value; 

    typedef typename Traits::FlowMap FlowMap;
    typedef typename Traits::Tolerance Tolerance;

    /// \brief \ref Exception for the case when the source equals the target.
    ///
    /// \ref Exception for the case when the source equals the target.
    ///
    class InvalidArgument : public lemon::LogicError {
    public:
      virtual const char* what() const throw() {
	return "lemon::EdmondsKarp::InvalidArgument";
      }
    };


  private:

    GRAPH_TYPEDEFS(typename Graph);
    typedef typename Graph::template NodeMap<Edge> PredMap;
    
    const Graph& _graph;
    const CapacityMap* _capacity;

    Node _source, _target;

    FlowMap* _flow;
    bool _local_flow;

    PredMap* _pred;
    std::vector<Node> _queue;
    
    Tolerance _tolerance;
    Value _flow_value;

    void createStructures() {
      if (!_flow) {
	_flow = Traits::createFlowMap(_graph);
	_local_flow = true;
      }
      if (!_pred) {
	_pred = new PredMap(_graph);
      }
      _queue.resize(countNodes(_graph));
    }

    void destroyStructures() {
      if (_local_flow) {
	delete _flow;
      }
      if (_pred) {
	delete _pred;
      }
    }
    
  public:

    ///\name Named template parameters

    ///@{

    template <typename _FlowMap>
    struct DefFlowMapTraits : public Traits {
      typedef _FlowMap FlowMap;
      static FlowMap *createFlowMap(const Graph&) {
	throw UninitializedParameter();
      }
    };

    /// \brief \ref named-templ-param "Named parameter" for setting
    /// FlowMap type
    ///
    /// \ref named-templ-param "Named parameter" for setting FlowMap
    /// type
    template <typename _FlowMap>
    struct DefFlowMap 
      : public EdmondsKarp<Graph, CapacityMap, DefFlowMapTraits<_FlowMap> > {
      typedef EdmondsKarp<Graph, CapacityMap, DefFlowMapTraits<_FlowMap> > 
      Create;
    };


    /// @}

  protected:
    
    EdmondsKarp() {}

  public:

    /// \brief The constructor of the class.
    ///
    /// The constructor of the class. 
    /// \param graph The directed graph the algorithm runs on. 
    /// \param capacity The capacity of the edges. 
    /// \param source The source node.
    /// \param target The target node.
    EdmondsKarp(const Graph& graph, const CapacityMap& capacity,
		Node source, Node target)
      : _graph(graph), _capacity(&capacity), _source(source), _target(target),
	_flow(0), _local_flow(false), _pred(0), _tolerance(), _flow_value()
    {
      if (_source == _target) {
        throw InvalidArgument();
      }
    }

    /// \brief Destrcutor.
    ///
    /// Destructor.
    ~EdmondsKarp() {
      destroyStructures();
    }

    /// \brief Sets the capacity map.
    ///
    /// Sets the capacity map.
    /// \return \c (*this)
    EdmondsKarp& capacityMap(const CapacityMap& map) {
      _capacity = &map;
      return *this;
    }

    /// \brief Sets the flow map.
    ///
    /// Sets the flow map.
    /// \return \c (*this)
    EdmondsKarp& flowMap(FlowMap& map) {
      if (_local_flow) {
	delete _flow;
	_local_flow = false;
      }
      _flow = &map;
      return *this;
    }

    /// \brief Returns the flow map.
    ///
    /// \return The flow map.
    const FlowMap& flowMap() {
      return *_flow;
    }

    /// \brief Sets the source node.
    ///
    /// Sets the source node.
    /// \return \c (*this)
    EdmondsKarp& source(const Node& node) {
      _source = node;
      return *this;
    }

    /// \brief Sets the target node.
    ///
    /// Sets the target node.
    /// \return \c (*this)
    EdmondsKarp& target(const Node& node) {
      _target = node;
      return *this;
    }

    /// \brief Sets the tolerance used by algorithm.
    ///
    /// Sets the tolerance used by algorithm.
    EdmondsKarp& tolerance(const Tolerance& tolerance) const {
      _tolerance = tolerance;
      return *this;
    } 

    /// \brief Returns the tolerance used by algorithm.
    ///
    /// Returns the tolerance used by algorithm.
    const Tolerance& tolerance() const {
      return tolerance;
    } 

    /// \name Execution control The simplest way to execute the
    /// algorithm is to use the \c run() member functions.
    /// \n
    /// If you need more control on initial solution or
    /// execution then you have to call one \ref init() function and then
    /// the start() or multiple times the \c augment() member function.  
    
    ///@{

    /// \brief Initializes the algorithm
    /// 
    /// It sets the flow to empty flow.
    void init() {
      createStructures();
      for (EdgeIt it(_graph); it != INVALID; ++it) {
        _flow->set(it, 0);
      }
      _flow_value = 0;
    }
    
    /// \brief Initializes the algorithm
    /// 
    /// Initializes the flow to the \c flowMap. The \c flowMap should
    /// contain a feasible flow, ie. in each node excluding the source
    /// and the target the incoming flow should be equal to the
    /// outgoing flow.
    template <typename FlowMap>
    void flowInit(const FlowMap& flowMap) {
      createStructures();
      for (EdgeIt e(_graph); e != INVALID; ++e) {
	_flow->set(e, flowMap[e]);
      }
      _flow_value = 0;
      for (OutEdgeIt jt(_graph, _source); jt != INVALID; ++jt) {
        _flow_value += (*_flow)[jt];
      }
      for (InEdgeIt jt(_graph, _source); jt != INVALID; ++jt) {
        _flow_value -= (*_flow)[jt];
      }
    }

    /// \brief Initializes the algorithm
    /// 
    /// Initializes the flow to the \c flowMap. The \c flowMap should
    /// contain a feasible flow, ie. in each node excluding the source
    /// and the target the incoming flow should be equal to the
    /// outgoing flow.  
    /// \return %False when the given flowMap does not contain
    /// feasible flow.
    template <typename FlowMap>
    bool checkedFlowInit(const FlowMap& flowMap) {
      createStructures();
      for (EdgeIt e(_graph); e != INVALID; ++e) {
	_flow->set(e, flowMap[e]);
      }
      for (NodeIt it(_graph); it != INVALID; ++it) {
        if (it == _source || it == _target) continue;
        Value outFlow = 0;
        for (OutEdgeIt jt(_graph, it); jt != INVALID; ++jt) {
          outFlow += (*_flow)[jt];
        }
        Value inFlow = 0;
        for (InEdgeIt jt(_graph, it); jt != INVALID; ++jt) {
          inFlow += (*_flow)[jt];
        }
        if (_tolerance.different(outFlow, inFlow)) {
          return false;
        }
      }
      for (EdgeIt it(_graph); it != INVALID; ++it) {
        if (_tolerance.less((*_flow)[it], 0)) return false;
        if (_tolerance.less((*_capacity)[it], (*_flow)[it])) return false;
      }
      _flow_value = 0;
      for (OutEdgeIt jt(_graph, _source); jt != INVALID; ++jt) {
        _flow_value += (*_flow)[jt];
      }
      for (InEdgeIt jt(_graph, _source); jt != INVALID; ++jt) {
        _flow_value -= (*_flow)[jt];
      }
      return true;
    }

    /// \brief Augment the solution on an edge shortest path.
    /// 
    /// Augment the solution on an edge shortest path. It search an
    /// edge shortest path between the source and the target
    /// in the residual graph with the bfs algoritm.
    /// Then it increase the flow on this path with the minimal residual
    /// capacity on the path. If there is not such path it gives back
    /// false.
    /// \return %False when the augmenting is not success so the
    /// current flow is a feasible and optimal solution.
    bool augment() {
      for (NodeIt n(_graph); n != INVALID; ++n) {
	_pred->set(n, INVALID);
      }
      
      int first = 0, last = 1;
      
      _queue[0] = _source;
      _pred->set(_source, OutEdgeIt(_graph, _source));

      while (first != last && (*_pred)[_target] == INVALID) {
	Node n = _queue[first++];
	
	for (OutEdgeIt e(_graph, n); e != INVALID; ++e) {
	  Value rem = (*_capacity)[e] - (*_flow)[e];
	  Node t = _graph.target(e);
	  if (_tolerance.positive(rem) && (*_pred)[t] == INVALID) {
	    _pred->set(t, e);
	    _queue[last++] = t;
	  }
	}
	for (InEdgeIt e(_graph, n); e != INVALID; ++e) {
	  Value rem = (*_flow)[e];
	  Node t = _graph.source(e);
	  if (_tolerance.positive(rem) && (*_pred)[t] == INVALID) {
	    _pred->set(t, e);
	    _queue[last++] = t;
	  }
	}
      }

      if ((*_pred)[_target] != INVALID) {
	Node n = _target;
	Edge e = (*_pred)[n];

	Value prem = (*_capacity)[e] - (*_flow)[e];
	n = _graph.source(e);
	while (n != _source) {
	  e = (*_pred)[n];
	  if (_graph.target(e) == n) {
	    Value rem = (*_capacity)[e] - (*_flow)[e];
	    if (rem < prem) prem = rem;
	    n = _graph.source(e);
	  } else {
	    Value rem = (*_flow)[e];
	    if (rem < prem) prem = rem;
	    n = _graph.target(e);   
	  } 
	}

	n = _target;
	e = (*_pred)[n];

	_flow->set(e, (*_flow)[e] + prem);
	n = _graph.source(e);
	while (n != _source) {
	  e = (*_pred)[n];
	  if (_graph.target(e) == n) {
	    _flow->set(e, (*_flow)[e] + prem);
	    n = _graph.source(e);
	  } else {
	    _flow->set(e, (*_flow)[e] - prem);
	    n = _graph.target(e);   
	  } 
	}

	_flow_value += prem;	
	return true;
      } else {
	return false;
      }
    }

    /// \brief Executes the algorithm
    ///
    /// It runs augmenting phases until the optimal solution is reached. 
    void start() {
      while (augment()) {}
    }

    /// \brief runs the algorithm.
    /// 
    /// It is just a shorthand for:
    ///
    ///\code 
    /// ek.init();
    /// ek.start();
    ///\endcode
    void run() {
      init();
      start();
    }

    /// @}

    /// \name Query Functions
    /// The result of the Edmonds-Karp algorithm can be obtained using these
    /// functions.\n
    /// Before the use of these functions,
    /// either run() or start() must be called.
    
    ///@{

    /// \brief Returns the value of the maximum flow.
    ///
    /// Returns the value of the maximum flow by returning the excess
    /// of the target node \c t. This value equals to the value of
    /// the maximum flow already after the first phase.
    Value flowValue() const {
      return _flow_value;
    }


    /// \brief Returns the flow on the edge.
    ///
    /// Sets the \c flowMap to the flow on the edges. This method can
    /// be called after the second phase of algorithm.
    Value flow(const Edge& edge) const {
      return (*_flow)[edge];
    }

    /// \brief Returns true when the node is on the source side of minimum cut.
    ///

    /// Returns true when the node is on the source side of minimum
    /// cut. This method can be called both after running \ref
    /// startFirstPhase() and \ref startSecondPhase().
    bool minCut(const Node& node) const {
      return (*_pred)[node] != INVALID;
    }

    /// \brief Returns a minimum value cut.
    ///
    /// Sets \c cut to the characteristic vector of a minimum value cut
    /// It simply calls the minMinCut member.
    /// \retval cut Write node bool map. 
    template <typename CutMap>
    void minCutMap(CutMap& cutMap) const {
      for (NodeIt n(_graph); n != INVALID; ++n) {
	cutMap.set(n, (*_pred)[n] != INVALID);
      }
      cutMap.set(_source, true);
    }    

    /// @}

  };

}

#endif
