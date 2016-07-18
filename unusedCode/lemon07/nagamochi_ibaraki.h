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

#ifndef LEMON_NAGAMOCHI_IBARAKI_H
#define LEMON_NAGAMOCHI_IBARAKI_H


/// \ingroup min_cut
/// \file 
/// \brief Maximum cardinality search and minimum cut in undirected
/// graphs.

#include <lemon/list_graph.h>
#include <lemon/bin_heap.h>
#include <lemon/bucket_heap.h>

#include <lemon/unionfind.h>
#include <lemon/topology.h>

#include <lemon/bits/invalid.h>
#include <lemon/error.h>
#include <lemon/maps.h>

#include <functional>

#include <lemon/graph_writer.h>
#include <lemon/time_measure.h>

namespace lemon {

  namespace _min_cut_bits {

    template <typename CapacityMap>
    struct HeapSelector {
      template <typename Value, typename Ref>
      struct Selector {
        typedef BinHeap<Value, Ref, std::greater<Value> > Heap;
      };
    };

    template <typename CapacityKey>
    struct HeapSelector<ConstMap<CapacityKey, Const<int, 1> > > {
      template <typename Value, typename Ref>
      struct Selector {
        typedef BucketHeap<Ref, false > Heap;
      };
    };

  }

  /// \brief Default traits class of MaxCardinalitySearch class.
  ///
  /// Default traits class of MaxCardinalitySearch class.
  /// \param Graph Graph type.
  /// \param CapacityMap Type of length map.
  template <typename _Graph, typename _CapacityMap>
  struct MaxCardinalitySearchDefaultTraits {
    /// The graph type the algorithm runs on. 
    typedef _Graph Graph;

    /// \brief The type of the map that stores the edge capacities.
    ///
    /// The type of the map that stores the edge capacities.
    /// It must meet the \ref concepts::ReadMap "ReadMap" concept.
    typedef _CapacityMap CapacityMap;

    /// \brief The type of the capacity of the edges.
    typedef typename CapacityMap::Value Value;

    /// \brief The cross reference type used by heap.
    ///
    /// The cross reference type used by heap.
    /// Usually it is \c Graph::NodeMap<int>.
    typedef typename Graph::template NodeMap<int> HeapCrossRef;

    /// \brief Instantiates a HeapCrossRef.
    ///
    /// This function instantiates a \ref HeapCrossRef. 
    /// \param graph is the graph, to which we would like to define the 
    /// HeapCrossRef.
    static HeapCrossRef *createHeapCrossRef(const Graph &graph) {
      return new HeapCrossRef(graph);
    }
    
    /// \brief The heap type used by MaxCardinalitySearch algorithm.
    ///
    /// The heap type used by MaxCardinalitySearch algorithm. It should
    /// maximalize the priorities. The default heap type is
    /// the \ref BinHeap, but it is specialized when the
    /// CapacityMap is ConstMap<Graph::Node, Const<int, 1> >
    /// to BucketHeap.
    ///
    /// \sa MaxCardinalitySearch
    typedef typename _min_cut_bits
    ::HeapSelector<CapacityMap>
    ::template Selector<Value, HeapCrossRef>
    ::Heap Heap;

    /// \brief Instantiates a Heap.
    ///
    /// This function instantiates a \ref Heap. 
    /// \param crossref The cross reference of the heap.
    static Heap *createHeap(HeapCrossRef& crossref) {
      return new Heap(crossref);
    }

    /// \brief The type of the map that stores whether a nodes is processed.
    ///
    /// The type of the map that stores whether a nodes is processed.
    /// It must meet the \ref concepts::WriteMap "WriteMap" concept.
    /// By default it is a NullMap.
    typedef NullMap<typename Graph::Node, bool> ProcessedMap;

    /// \brief Instantiates a ProcessedMap.
    ///
    /// This function instantiates a \ref ProcessedMap. 
    /// \param graph is the graph, to which
    /// we would like to define the \ref ProcessedMap
#ifdef DOXYGEN
    static ProcessedMap *createProcessedMap(const Graph &graph)
#else
    static ProcessedMap *createProcessedMap(const Graph &)
#endif
    {
      return new ProcessedMap();
    }

    /// \brief The type of the map that stores the cardinalties of the nodes.
    /// 
    /// The type of the map that stores the cardinalities of the nodes.
    /// It must meet the \ref concepts::WriteMap "WriteMap" concept.
    typedef typename Graph::template NodeMap<Value> CardinalityMap;

    /// \brief Instantiates a CardinalityMap.
    ///
    /// This function instantiates a \ref CardinalityMap. 
    /// \param graph is the graph, to which we would like to define the \ref 
    /// CardinalityMap
    static CardinalityMap *createCardinalityMap(const Graph &graph) {
      return new CardinalityMap(graph);
    }


  };
  
  /// \ingroup search
  ///
  /// \brief Maximum Cardinality Search algorithm class.
  ///
  /// This class provides an efficient implementation of Maximum Cardinality 
  /// Search algorithm. The maximum cardinality search chooses first time any 
  /// node of the graph. Then every time it chooses the node which is connected
  /// to the processed nodes at most in the sum of capacities on the out 
  /// edges. If there is a cut in the graph the algorithm should choose
  /// again any unprocessed node of the graph. Each node cardinality is
  /// the sum of capacities on the out edges to the nodes which are processed
  /// before the given node.
  ///
  /// The edge capacities are passed to the algorithm using a
  /// \ref concepts::ReadMap "ReadMap", so it is easy to change it to any 
  /// kind of capacity.
  ///
  /// The type of the capacity is determined by the \ref 
  /// concepts::ReadMap::Value "Value" of the capacity map.
  ///
  /// It is also possible to change the underlying priority heap.
  ///
  ///
  /// \param _Graph The graph type the algorithm runs on. The default value
  /// is \ref ListGraph. The value of Graph is not used directly by
  /// the search algorithm, it is only passed to 
  /// \ref MaxCardinalitySearchDefaultTraits.
  /// \param _CapacityMap This read-only EdgeMap determines the capacities of 
  /// the edges. It is read once for each edge, so the map may involve in
  /// relatively time consuming process to compute the edge capacity if
  /// it is necessary. The default map type is \ref
  /// concepts::Graph::EdgeMap "Graph::EdgeMap<int>".  The value
  /// of CapacityMap is not used directly by search algorithm, it is only 
  /// passed to \ref MaxCardinalitySearchDefaultTraits.  
  /// \param _Traits Traits class to set various data types used by the 
  /// algorithm.  The default traits class is 
  /// \ref MaxCardinalitySearchDefaultTraits 
  /// "MaxCardinalitySearchDefaultTraits<_Graph, _CapacityMap>".  
  /// See \ref MaxCardinalitySearchDefaultTraits 
  /// for the documentation of a MaxCardinalitySearch traits class.
  ///
  /// \author Balazs Dezso

#ifdef DOXYGEN
  template <typename _Graph, typename _CapacityMap, typename _Traits>
#else
  template <typename _Graph = ListUGraph,
	    typename _CapacityMap = typename _Graph::template EdgeMap<int>,
	    typename _Traits = 
            MaxCardinalitySearchDefaultTraits<_Graph, _CapacityMap> >
#endif
  class MaxCardinalitySearch {
  public:
    /// \brief \ref Exception for uninitialized parameters.
    ///
    /// This error represents problems in the initialization
    /// of the parameters of the algorithms.
    class UninitializedParameter : public lemon::UninitializedParameter {
    public:
      virtual const char* what() const throw() {
	return "lemon::MaxCardinalitySearch::UninitializedParameter";
      }
    };

    typedef _Traits Traits;
    ///The type of the underlying graph.
    typedef typename Traits::Graph Graph;
    
    ///The type of the capacity of the edges.
    typedef typename Traits::CapacityMap::Value Value;
    ///The type of the map that stores the edge capacities.
    typedef typename Traits::CapacityMap CapacityMap;
    ///The type of the map indicating if a node is processed.
    typedef typename Traits::ProcessedMap ProcessedMap;
    ///The type of the map that stores the cardinalities of the nodes.
    typedef typename Traits::CardinalityMap CardinalityMap;
    ///The cross reference type used for the current heap.
    typedef typename Traits::HeapCrossRef HeapCrossRef;
    ///The heap type used by the algorithm. It maximize the priorities.
    typedef typename Traits::Heap Heap;
  private:
    /// Pointer to the underlying graph.
    const Graph *_graph;
    /// Pointer to the capacity map
    const CapacityMap *_capacity;
    ///Pointer to the map of cardinality.
    CardinalityMap *_cardinality;
    ///Indicates if \ref _cardinality is locally allocated (\c true) or not.
    bool local_cardinality;
    ///Pointer to the map of processed status of the nodes.
    ProcessedMap *_processed;
    ///Indicates if \ref _processed is locally allocated (\c true) or not.
    bool local_processed;
    ///Pointer to the heap cross references.
    HeapCrossRef *_heap_cross_ref;
    ///Indicates if \ref _heap_cross_ref is locally allocated (\c true) or not.
    bool local_heap_cross_ref;
    ///Pointer to the heap.
    Heap *_heap;
    ///Indicates if \ref _heap is locally allocated (\c true) or not.
    bool local_heap;

  public :

    typedef MaxCardinalitySearch Create;
 
    ///\name Named template parameters

    ///@{

    template <class T>
    struct DefCardinalityMapTraits : public Traits {
      typedef T CardinalityMap;
      static CardinalityMap *createCardinalityMap(const Graph &) 
      {
	throw UninitializedParameter();
      }
    };
    /// \brief \ref named-templ-param "Named parameter" for setting 
    /// CardinalityMap type
    ///
    /// \ref named-templ-param "Named parameter" for setting CardinalityMap 
    /// type
    template <class T>
    struct DefCardinalityMap 
      : public MaxCardinalitySearch<Graph, CapacityMap, 
                                    DefCardinalityMapTraits<T> > { 
      typedef MaxCardinalitySearch<Graph, CapacityMap, 
                                   DefCardinalityMapTraits<T> > Create;
    };
    
    template <class T>
    struct DefProcessedMapTraits : public Traits {
      typedef T ProcessedMap;
      static ProcessedMap *createProcessedMap(const Graph &) {
	throw UninitializedParameter();
      }
    };
    /// \brief \ref named-templ-param "Named parameter" for setting 
    /// ProcessedMap type
    ///
    /// \ref named-templ-param "Named parameter" for setting ProcessedMap type
    ///
    template <class T>
    struct DefProcessedMap 
      : public MaxCardinalitySearch<Graph, CapacityMap, 
                                    DefProcessedMapTraits<T> > { 
      typedef MaxCardinalitySearch<Graph, CapacityMap, 
                                   DefProcessedMapTraits<T> > Create;
    };
    
    template <class H, class CR>
    struct DefHeapTraits : public Traits {
      typedef CR HeapCrossRef;
      typedef H Heap;
      static HeapCrossRef *createHeapCrossRef(const Graph &) {
	throw UninitializedParameter();
      }
      static Heap *createHeap(HeapCrossRef &) {
	throw UninitializedParameter();
      }
    };
    /// \brief \ref named-templ-param "Named parameter" for setting heap 
    /// and cross reference type
    ///
    /// \ref named-templ-param "Named parameter" for setting heap and cross 
    /// reference type
    template <class H, class CR = typename Graph::template NodeMap<int> >
    struct DefHeap
      : public MaxCardinalitySearch<Graph, CapacityMap, 
                                    DefHeapTraits<H, CR> > { 
      typedef MaxCardinalitySearch< Graph, CapacityMap, 
                                    DefHeapTraits<H, CR> > Create;
    };

    template <class H, class CR>
    struct DefStandardHeapTraits : public Traits {
      typedef CR HeapCrossRef;
      typedef H Heap;
      static HeapCrossRef *createHeapCrossRef(const Graph &graph) {
	return new HeapCrossRef(graph);
      }
      static Heap *createHeap(HeapCrossRef &crossref) {
	return new Heap(crossref);
      }
    };

    /// \brief \ref named-templ-param "Named parameter" for setting heap and 
    /// cross reference type with automatic allocation
    ///
    /// \ref named-templ-param "Named parameter" for setting heap and cross 
    /// reference type. It can allocate the heap and the cross reference 
    /// object if the cross reference's constructor waits for the graph as 
    /// parameter and the heap's constructor waits for the cross reference.
    template <class H, class CR = typename Graph::template NodeMap<int> >
    struct DefStandardHeap
      : public MaxCardinalitySearch<Graph, CapacityMap, 
                                    DefStandardHeapTraits<H, CR> > { 
      typedef MaxCardinalitySearch<Graph, CapacityMap, 
                                   DefStandardHeapTraits<H, CR> > 
      Create;
    };
    
    ///@}


  protected:

    MaxCardinalitySearch() {}

  public:      
    
    /// \brief Constructor.
    ///
    ///\param graph the graph the algorithm will run on.
    ///\param capacity the capacity map used by the algorithm.
    MaxCardinalitySearch(const Graph& graph, const CapacityMap& capacity) :
      _graph(&graph), _capacity(&capacity),
      _cardinality(0), local_cardinality(false),
      _processed(0), local_processed(false),
      _heap_cross_ref(0), local_heap_cross_ref(false),
      _heap(0), local_heap(false)
    { }
    
    /// \brief Destructor.
    ~MaxCardinalitySearch() {
      if(local_cardinality) delete _cardinality;
      if(local_processed) delete _processed;
      if(local_heap_cross_ref) delete _heap_cross_ref;
      if(local_heap) delete _heap;
    }

    /// \brief Sets the capacity map.
    ///
    /// Sets the capacity map.
    /// \return <tt> (*this) </tt>
    MaxCardinalitySearch &capacityMap(const CapacityMap &m) {
      _capacity = &m;
      return *this;
    }

    /// \brief Sets the map storing the cardinalities calculated by the 
    /// algorithm.
    ///
    /// Sets the map storing the cardinalities calculated by the algorithm.
    /// If you don't use this function before calling \ref run(),
    /// it will allocate one. The destuctor deallocates this
    /// automatically allocated map, of course.
    /// \return <tt> (*this) </tt>
    MaxCardinalitySearch &cardinalityMap(CardinalityMap &m) {
      if(local_cardinality) {
	delete _cardinality;
	local_cardinality=false;
      }
      _cardinality = &m;
      return *this;
    }

    /// \brief Sets the map storing the processed nodes.
    ///
    /// Sets the map storing the processed nodes.
    /// If you don't use this function before calling \ref run(),
    /// it will allocate one. The destuctor deallocates this
    /// automatically allocated map, of course.
    /// \return <tt> (*this) </tt>
    MaxCardinalitySearch &processedMap(ProcessedMap &m) 
    {
      if(local_processed) {
	delete _processed;
	local_processed=false;
      }
      _processed = &m;
      return *this;
    }

    /// \brief Sets the heap and the cross reference used by algorithm.
    ///
    /// Sets the heap and the cross reference used by algorithm.
    /// If you don't use this function before calling \ref run(),
    /// it will allocate one. The destuctor deallocates this
    /// automatically allocated map, of course.
    /// \return <tt> (*this) </tt>
    MaxCardinalitySearch &heap(Heap& hp, HeapCrossRef &cr) {
      if(local_heap_cross_ref) {
	delete _heap_cross_ref;
	local_heap_cross_ref = false;
      }
      _heap_cross_ref = &cr;
      if(local_heap) {
	delete _heap;
	local_heap = false;
      }
      _heap = &hp;
      return *this;
    }

  private:

    typedef typename Graph::Node Node;
    typedef typename Graph::NodeIt NodeIt;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::InEdgeIt InEdgeIt;

    void create_maps() {
      if(!_cardinality) {
	local_cardinality = true;
	_cardinality = Traits::createCardinalityMap(*_graph);
      }
      if(!_processed) {
	local_processed = true;
	_processed = Traits::createProcessedMap(*_graph);
      }
      if (!_heap_cross_ref) {
	local_heap_cross_ref = true;
	_heap_cross_ref = Traits::createHeapCrossRef(*_graph);
      }
      if (!_heap) {
	local_heap = true;
	_heap = Traits::createHeap(*_heap_cross_ref);
      }
    }
    
    void finalizeNodeData(Node node, Value capacity) {
      _processed->set(node, true);
      _cardinality->set(node, capacity);
    }

  public:
    /// \name Execution control
    /// The simplest way to execute the algorithm is to use
    /// one of the member functions called \c run(...).
    /// \n
    /// If you need more control on the execution,
    /// first you must call \ref init(), then you can add several source nodes
    /// with \ref addSource().
    /// Finally \ref start() will perform the actual path
    /// computation.

    ///@{

    /// \brief Initializes the internal data structures.
    ///
    /// Initializes the internal data structures.
    void init() {
      create_maps();
      _heap->clear();
      for (NodeIt it(*_graph) ; it != INVALID ; ++it) {
	_processed->set(it, false);
	_heap_cross_ref->set(it, Heap::PRE_HEAP);
      }
    }
    
    /// \brief Adds a new source node.
    /// 
    /// Adds a new source node to the priority heap.
    ///
    /// It checks if the node has not yet been added to the heap.
    void addSource(Node source, Value capacity = 0) {
      if(_heap->state(source) == Heap::PRE_HEAP) {
	_heap->push(source, capacity);
      } 
    }
    
    /// \brief Processes the next node in the priority heap
    ///
    /// Processes the next node in the priority heap.
    ///
    /// \return The processed node.
    ///
    /// \warning The priority heap must not be empty!
    Node processNextNode() {
      Node node = _heap->top(); 
      finalizeNodeData(node, _heap->prio());
      _heap->pop();
      
      for (InEdgeIt it(*_graph, node); it != INVALID; ++it) {
	Node source = _graph->source(it); 
	switch (_heap->state(source)) {
	case Heap::PRE_HEAP:
	  _heap->push(source, (*_capacity)[it]); 
	  break;
	case Heap::IN_HEAP:
	  _heap->decrease(source, (*_heap)[source] + (*_capacity)[it]); 
	  break;
	case Heap::POST_HEAP:
	  break;
	}
      }
      return node;
    }

    /// \brief Next node to be processed.
    ///
    /// Next node to be processed.
    ///
    /// \return The next node to be processed or INVALID if the 
    /// priority heap is empty.
    Node nextNode() { 
      return _heap->empty() ? _heap->top() : INVALID;
    }
 
    /// \brief Returns \c false if there are nodes
    /// to be processed in the priority heap
    ///
    /// Returns \c false if there are nodes
    /// to be processed in the priority heap
    bool emptyQueue() { return _heap->empty(); }
    /// \brief Returns the number of the nodes to be processed 
    /// in the priority heap
    ///
    /// Returns the number of the nodes to be processed in the priority heap
    int queueSize() { return _heap->size(); }
    
    /// \brief Executes the algorithm.
    ///
    /// Executes the algorithm.
    ///
    ///\pre init() must be called and at least one node should be added
    /// with addSource() before using this function.
    ///
    /// This method runs the Maximum Cardinality Search algorithm from the 
    /// source node(s).
    void start() {
      while ( !_heap->empty() ) processNextNode();
    }
    
    /// \brief Executes the algorithm until \c dest is reached.
    ///
    /// Executes the algorithm until \c dest is reached.
    ///
    /// \pre init() must be called and at least one node should be added
    /// with addSource() before using this function.
    ///
    /// This method runs the %MaxCardinalitySearch algorithm from the source 
    /// nodes.
    void start(Node dest) {
      while ( !_heap->empty() && _heap->top()!=dest ) processNextNode();
      if ( !_heap->empty() ) finalizeNodeData(_heap->top(), _heap->prio());
    }
    
    /// \brief Executes the algorithm until a condition is met.
    ///
    /// Executes the algorithm until a condition is met.
    ///
    /// \pre init() must be called and at least one node should be added
    /// with addSource() before using this function.
    ///
    /// \param nm must be a bool (or convertible) node map. The algorithm
    /// will stop when it reaches a node \c v with <tt>nm[v]==true</tt>.
    template <typename NodeBoolMap>
    void start(const NodeBoolMap &nm) {
      while ( !_heap->empty() && !nm[_heap->top()] ) processNextNode();
      if ( !_heap->empty() ) finalizeNodeData(_heap->top(),_heap->prio());
    }
    
    /// \brief Runs the maximal cardinality search algorithm from node \c s.
    ///
    /// This method runs the %MaxCardinalitySearch algorithm from a root 
    /// node \c s.
    ///
    ///\note d.run(s) is just a shortcut of the following code.
    ///\code
    ///  d.init();
    ///  d.addSource(s);
    ///  d.start();
    ///\endcode
    void run(Node s) {
      init();
      addSource(s);
      start();
    }

    /// \brief Runs the maximal cardinality search algorithm for the 
    /// whole graph.
    ///
    /// This method runs the %MaxCardinalitySearch algorithm from all 
    /// unprocessed node of the graph.
    ///
    ///\note d.run(s) is just a shortcut of the following code.
    ///\code
    ///  d.init();
    ///  for (NodeIt it(graph); it != INVALID; ++it) {
    ///    if (!d.reached(it)) {
    ///      d.addSource(s);
    ///      d.start();
    ///    }
    ///  }
    ///\endcode
    void run() {
      init();
      for (NodeIt it(*_graph); it != INVALID; ++it) {
        if (!reached(it)) {
          addSource(it);
          start();
        }
      }
    }
    
    ///@}

    /// \name Query Functions
    /// The result of the maximum cardinality search algorithm can be 
    /// obtained using these functions.
    /// \n
    /// Before the use of these functions, either run() or start() must be 
    /// called.
    
    ///@{

    /// \brief The cardinality of a node.
    ///
    /// Returns the cardinality of a node.
    /// \pre \ref run() must be called before using this function.
    /// \warning If node \c v in unreachable from the root the return value
    /// of this funcion is undefined.
    Value cardinality(Node node) const { return (*_cardinality)[node]; }

    /// \brief The current cardinality of a node.
    ///
    /// Returns the current cardinality of a node.
    /// \pre the given node should be reached but not processed
    Value currentCardinality(Node node) const { return (*_heap)[node]; }

    /// \brief Returns a reference to the NodeMap of cardinalities.
    ///
    /// Returns a reference to the NodeMap of cardinalities. \pre \ref run() 
    /// must be called before using this function.
    const CardinalityMap &cardinalityMap() const { return *_cardinality;}
 
    /// \brief Checks if a node is reachable from the root.
    ///
    /// Returns \c true if \c v is reachable from the root.
    /// \warning The source nodes are inditated as unreached.
    /// \pre \ref run() must be called before using this function.
    bool reached(Node v) { return (*_heap_cross_ref)[v] != Heap::PRE_HEAP; }

    /// \brief Checks if a node is processed.
    ///
    /// Returns \c true if \c v is processed, i.e. the shortest
    /// path to \c v has already found.
    /// \pre \ref run() must be called before using this function.
    bool processed(Node v) { return (*_heap_cross_ref)[v] == Heap::POST_HEAP; }
    
    ///@}
  };

  /// \brief Default traits class of NagamochiIbaraki class.
  ///
  /// Default traits class of NagamochiIbaraki class.
  /// \param Graph Graph type.
  /// \param CapacityMap Type of length map.
  template <typename _Graph, typename _CapacityMap>
  struct NagamochiIbarakiDefaultTraits {
    /// \brief The type of the capacity of the edges.
    typedef typename _CapacityMap::Value Value;

    /// The graph type the algorithm runs on. 
    typedef _Graph Graph;

    /// The AuxGraph type which is an Graph
    typedef ListUGraph AuxGraph;

    /// \brief Instantiates a AuxGraph.
    ///
    /// This function instantiates a \ref AuxGraph. 
    static AuxGraph *createAuxGraph() {
      return new AuxGraph();
    }

    /// \brief The type of the map that stores the edge capacities.
    ///
    /// The type of the map that stores the edge capacities.
    /// It must meet the \ref concepts::ReadMap "ReadMap" concept.
    typedef _CapacityMap CapacityMap;

    /// \brief Instantiates a CapacityMap.
    ///
    /// This function instantiates a \ref CapacityMap.
#ifdef DOXYGEN
    static CapacityMap *createCapacityMap(const Graph& graph) 
#else
    static CapacityMap *createCapacityMap(const Graph&)
#endif
    {
      throw UninitializedParameter();
    }

    /// \brief The CutValueMap type
    ///
    /// The type of the map that stores the cut value of a node.
    typedef AuxGraph::NodeMap<Value> AuxCutValueMap;

    /// \brief Instantiates a AuxCutValueMap.
    ///
    /// This function instantiates a \ref AuxCutValueMap. 
    static AuxCutValueMap *createAuxCutValueMap(const AuxGraph& graph) {
      return new AuxCutValueMap(graph);
    }

    /// \brief The AuxCapacityMap type
    ///
    /// The type of the map that stores the auxiliary edge capacities.
    typedef AuxGraph::UEdgeMap<Value> AuxCapacityMap;

    /// \brief Instantiates a AuxCapacityMap.
    ///
    /// This function instantiates a \ref AuxCapacityMap. 
    static AuxCapacityMap *createAuxCapacityMap(const AuxGraph& graph) {
      return new AuxCapacityMap(graph);
    }

    /// \brief The cross reference type used by heap.
    ///
    /// The cross reference type used by heap.
    /// Usually it is \c Graph::NodeMap<int>.
    typedef AuxGraph::NodeMap<int> HeapCrossRef;

    /// \brief Instantiates a HeapCrossRef.
    ///
    /// This function instantiates a \ref HeapCrossRef. 
    /// \param graph is the graph, to which we would like to define the 
    /// HeapCrossRef.
    static HeapCrossRef *createHeapCrossRef(const AuxGraph &graph) {
      return new HeapCrossRef(graph);
    }
    
    /// \brief The heap type used by NagamochiIbaraki algorithm.
    ///
    /// The heap type used by NagamochiIbaraki algorithm. It should
    /// maximalize the priorities and the heap's key type is
    /// the aux graph's node.
    ///
    /// \sa BinHeap
    /// \sa NagamochiIbaraki
    typedef typename _min_cut_bits
    ::HeapSelector<CapacityMap>
    ::template Selector<Value, HeapCrossRef>
    ::Heap Heap;
    
    /// \brief Instantiates a Heap.
    ///
    /// This function instantiates a \ref Heap. 
    /// \param crossref The cross reference of the heap.
    static Heap *createHeap(HeapCrossRef& crossref) {
      return new Heap(crossref);
    }

    /// \brief Map from the AuxGraph's node type to the Graph's node type.
    ///
    /// Map from the AuxGraph's node type to the Graph's node type.
    typedef typename AuxGraph
    ::template NodeMap<typename Graph::Node> NodeRefMap;

    /// \brief Instantiates a NodeRefMap.
    ///
    /// This function instantiates a \ref NodeRefMap. 
    static NodeRefMap *createNodeRefMap(const AuxGraph& graph) {
      return new NodeRefMap(graph);
    }

    /// \brief Map from the Graph's node type to the Graph's node type.
    ///
    /// Map from the Graph's node type to the Graph's node type.
    typedef typename Graph
    ::template NodeMap<typename Graph::Node> ListRefMap;

    /// \brief Instantiates a ListRefMap.
    ///
    /// This function instantiates a \ref ListRefMap. 
    static ListRefMap *createListRefMap(const Graph& graph) {
      return new ListRefMap(graph);
    }
    

  };

  /// \ingroup min_cut
  ///
  /// \brief Calculates the minimum cut in an undirected graph.
  ///
  /// Calculates the minimum cut in an undirected graph with the
  /// Nagamochi-Ibaraki algorithm. The algorithm separates the graph's
  /// nodes into two partitions with the minimum sum of edge capacities
  /// between the two partitions. The algorithm can be used to test
  /// the network reliability specifically to test how many links have
  /// to be destroyed in the network to split it at least two
  /// distinict subnetwork.
  ///
  /// The complexity of the algorithm is \f$ O(ne\log(n)) \f$ but with
  /// Fibonacci heap it can be decreased to \f$ O(ne+n^2\log(n)) \f$.
  /// When unit capacity minimum cut is computed then it uses
  /// BucketHeap which results \f$ O(ne) \f$ time complexity.
  ///
  /// \warning The value type of the capacity map should be able to hold
  /// any cut value of the graph, otherwise the result can overflow.
#ifdef DOXYGEN
  template <typename _Graph, typename _CapacityMap, typename _Traits>
#else
  template <typename _Graph = ListUGraph, 
	    typename _CapacityMap = typename _Graph::template UEdgeMap<int>, 
	    typename _Traits 
            = NagamochiIbarakiDefaultTraits<_Graph, _CapacityMap> >
#endif
  class NagamochiIbaraki {
  public:
    /// \brief \ref Exception for uninitialized parameters.
    ///
    /// This error represents problems in the initialization
    /// of the parameters of the algorithms.
    class UninitializedParameter : public lemon::UninitializedParameter {
    public:
      virtual const char* what() const throw() {
	return "lemon::NagamochiIbaraki::UninitializedParameter";
      }
    };


  private:

    typedef _Traits Traits;
    /// The type of the underlying graph.
    typedef typename Traits::Graph Graph;
    
    /// The type of the capacity of the edges.
    typedef typename Traits::CapacityMap::Value Value;
    /// The type of the map that stores the edge capacities.
    typedef typename Traits::CapacityMap CapacityMap;
    /// The type of the aux graph
    typedef typename Traits::AuxGraph AuxGraph;
    /// The type of the aux capacity map
    typedef typename Traits::AuxCapacityMap AuxCapacityMap;
    /// The type of the aux cut value map
    typedef typename Traits::AuxCutValueMap AuxCutValueMap;
    /// The cross reference type used for the current heap.
    typedef typename Traits::HeapCrossRef HeapCrossRef;
    /// The heap type used by the max cardinality algorithm.
    typedef typename Traits::Heap Heap;
    /// The node refrefernces between the original and aux graph type.
    typedef typename Traits::NodeRefMap NodeRefMap;
    /// The list node refrefernces in the original graph type.
    typedef typename Traits::ListRefMap ListRefMap;

  public:

    ///\name Named template parameters

    ///@{

    struct DefUnitCapacityTraits : public Traits {
      typedef ConstMap<typename Graph::UEdge, Const<int, 1> > CapacityMap;
      static CapacityMap *createCapacityMap(const Graph&) {
	return new CapacityMap();
      }
    };
    /// \brief \ref named-templ-param "Named parameter" for setting  
    /// the capacity type to constMap<UEdge, int, 1>()
    ///
    /// \ref named-templ-param "Named parameter" for setting 
    /// the capacity type to constMap<UEdge, int, 1>()
    struct DefUnitCapacity
      : public NagamochiIbaraki<Graph, CapacityMap, 
                                DefUnitCapacityTraits> { 
      typedef NagamochiIbaraki<Graph, CapacityMap, 
                               DefUnitCapacityTraits> Create;
    };


    template <class H, class CR>
    struct DefHeapTraits : public Traits {
      typedef CR HeapCrossRef;
      typedef H Heap;
      static HeapCrossRef *createHeapCrossRef(const AuxGraph &) {
	throw UninitializedParameter();
      }
      static Heap *createHeap(HeapCrossRef &) {
	throw UninitializedParameter();
      }
    };
    /// \brief \ref named-templ-param "Named parameter" for setting heap 
    /// and cross reference type
    ///
    /// \ref named-templ-param "Named parameter" for setting heap and cross 
    /// reference type
    template <class H, class CR = typename Graph::template NodeMap<int> >
    struct DefHeap
      : public NagamochiIbaraki<Graph, CapacityMap, 
                                DefHeapTraits<H, CR> > { 
      typedef NagamochiIbaraki< Graph, CapacityMap, 
                                DefHeapTraits<H, CR> > Create;
    };

    template <class H, class CR>
    struct DefStandardHeapTraits : public Traits {
      typedef CR HeapCrossRef;
      typedef H Heap;
      static HeapCrossRef *createHeapCrossRef(const AuxGraph &graph) {
	return new HeapCrossRef(graph);
      }
      static Heap *createHeap(HeapCrossRef &crossref) {
	return new Heap(crossref);
      }
    };

    /// \brief \ref named-templ-param "Named parameter" for setting heap and 
    /// cross reference type with automatic allocation
    ///
    /// \ref named-templ-param "Named parameter" for setting heap and cross 
    /// reference type. It can allocate the heap and the cross reference 
    /// object if the cross reference's constructor waits for the graph as 
    /// parameter and the heap's constructor waits for the cross reference.
    template <class H, class CR = typename Graph::template NodeMap<int> >
    struct DefStandardHeap
      : public NagamochiIbaraki<Graph, CapacityMap, 
                                DefStandardHeapTraits<H, CR> > { 
      typedef NagamochiIbaraki<Graph, CapacityMap, 
                               DefStandardHeapTraits<H, CR> > 
      Create;
    };

    ///@}


  private:
    /// Pointer to the underlying graph.
    const Graph *_graph;
    /// Pointer to the capacity map
    const CapacityMap *_capacity;
    /// \brief Indicates if \ref _capacity is locally allocated 
    /// (\c true) or not.
    bool local_capacity;

    /// Pointer to the aux graph.
    AuxGraph *_aux_graph;
    /// \brief Indicates if \ref _aux_graph is locally allocated 
    /// (\c true) or not.
    bool local_aux_graph;
    /// Pointer to the aux capacity map
    AuxCapacityMap *_aux_capacity;
    /// \brief Indicates if \ref _aux_capacity is locally allocated 
    /// (\c true) or not.
    bool local_aux_capacity;
    /// Pointer to the aux cut value map
    AuxCutValueMap *_aux_cut_value;
    /// \brief Indicates if \ref _aux_cut_value is locally allocated 
    /// (\c true) or not.
    bool local_aux_cut_value;
    /// Pointer to the heap cross references.
    HeapCrossRef *_heap_cross_ref;
    /// \brief Indicates if \ref _heap_cross_ref is locally allocated 
    /// (\c true) or not.
    bool local_heap_cross_ref;
    /// Pointer to the heap.
    Heap *_heap;
    /// Indicates if \ref _heap is locally allocated (\c true) or not.
    bool local_heap;

    /// The min cut value.
    Value _min_cut;
    /// The number of the nodes of the aux graph.
    int _node_num;
    /// The first and last node of the min cut in the next list.
    std::vector<typename Graph::Node> _cut;

    /// \brief The first and last element in the list associated
    /// to the aux graph node.
    NodeRefMap *_first, *_last;
    /// \brief The next node in the node lists.
    ListRefMap *_next;
    
    void createStructures() {
      if (!_capacity) {
        local_capacity = true;
        _capacity = Traits::createCapacityMap(*_graph);
      }
      if(!_aux_graph) {
	local_aux_graph = true;
	_aux_graph = Traits::createAuxGraph();
      }
      if(!_aux_capacity) {
	local_aux_capacity = true;
	_aux_capacity = Traits::createAuxCapacityMap(*_aux_graph);
      }
      if(!_aux_cut_value) {
	local_aux_cut_value = true;
	_aux_cut_value = Traits::createAuxCutValueMap(*_aux_graph);
      }

      _first = Traits::createNodeRefMap(*_aux_graph);
      _last = Traits::createNodeRefMap(*_aux_graph);

      _next = Traits::createListRefMap(*_graph);

      if (!_heap_cross_ref) {
	local_heap_cross_ref = true;
	_heap_cross_ref = Traits::createHeapCrossRef(*_aux_graph);
      }
      if (!_heap) {
	local_heap = true;
	_heap = Traits::createHeap(*_heap_cross_ref);
      }
    }

    void createAuxGraph() {
      typename Graph::template NodeMap<typename AuxGraph::Node> ref(*_graph);

      for (typename Graph::NodeIt n(*_graph); n != INVALID; ++n) {
        _next->set(n, INVALID);
        typename AuxGraph::Node node = _aux_graph->addNode();
        _first->set(node, n);
        _last->set(node, n);
        ref.set(n, node);
      }

      typename AuxGraph::template NodeMap<typename AuxGraph::UEdge> 
      edges(*_aux_graph, INVALID);

      for (typename Graph::NodeIt n(*_graph); n != INVALID; ++n) {
        for (typename Graph::IncEdgeIt e(*_graph, n); e != INVALID; ++e) {
          typename Graph::Node tn = _graph->runningNode(e);
          if (n < tn || n == tn) continue;
          if (edges[ref[tn]] != INVALID) {
            Value value = 
              (*_aux_capacity)[edges[ref[tn]]] + (*_capacity)[e];
            _aux_capacity->set(edges[ref[tn]], value);
          } else {
            edges.set(ref[tn], _aux_graph->addEdge(ref[n], ref[tn]));
            Value value = (*_capacity)[e];
            _aux_capacity->set(edges[ref[tn]], value);            
          }
        }
        for (typename Graph::IncEdgeIt e(*_graph, n); e != INVALID; ++e) {
          typename Graph::Node tn = _graph->runningNode(e);
          edges.set(ref[tn], INVALID);
        }
      }

      _cut.resize(1, INVALID);
      _min_cut = std::numeric_limits<Value>::max();
      for (typename AuxGraph::NodeIt n(*_aux_graph); n != INVALID; ++n) {
        Value value = 0;
        for (typename AuxGraph::IncEdgeIt e(*_aux_graph, n); 
             e != INVALID; ++e) {
          value += (*_aux_capacity)[e];
        }
        if (_min_cut > value) {
          _min_cut = value;
          _cut[0] = (*_first)[n];
        } 
        (*_aux_cut_value)[n] = value;
      }
    }
    

  public :

    typedef NagamochiIbaraki Create;


    /// \brief Constructor.
    ///
    ///\param graph the graph the algorithm will run on.
    ///\param capacity the capacity map used by the algorithm.
    NagamochiIbaraki(const Graph& graph, const CapacityMap& capacity) 
      : _graph(&graph), 
        _capacity(&capacity), local_capacity(false),
        _aux_graph(0), local_aux_graph(false),
        _aux_capacity(0), local_aux_capacity(false),
        _aux_cut_value(0), local_aux_cut_value(false),
        _heap_cross_ref(0), local_heap_cross_ref(false),
        _heap(0), local_heap(false),
        _first(0), _last(0), _next(0) {}

    /// \brief Constructor.
    ///
    /// This constructor can be used only when the Traits class
    /// defines how can we instantiate a local capacity map.
    /// If the DefUnitCapacity used the algorithm automatically
    /// construct the capacity map.
    ///
    ///\param graph the graph the algorithm will run on.
    NagamochiIbaraki(const Graph& graph) 
      : _graph(&graph), 
        _capacity(0), local_capacity(false),
        _aux_graph(0), local_aux_graph(false),
        _aux_capacity(0), local_aux_capacity(false),
        _aux_cut_value(0), local_aux_cut_value(false),
        _heap_cross_ref(0), local_heap_cross_ref(false),
        _heap(0), local_heap(false),
        _first(0), _last(0), _next(0) {}

    /// \brief Destructor.
    ///
    /// Destructor.
    ~NagamochiIbaraki() {
      if (local_heap) delete _heap;
      if (local_heap_cross_ref) delete _heap_cross_ref;
      if (_first) delete _first;
      if (_last) delete _last;
      if (_next) delete _next;
      if (local_aux_capacity) delete _aux_capacity;
      if (local_aux_cut_value) delete _aux_cut_value;
      if (local_aux_graph) delete _aux_graph;
      if (local_capacity) delete _capacity;
    }

    /// \brief Sets the heap and the cross reference used by algorithm.
    ///
    /// Sets the heap and the cross reference used by algorithm.
    /// If you don't use this function before calling \ref run(),
    /// it will allocate one. The destuctor deallocates this
    /// automatically allocated heap and cross reference, of course.
    /// \return <tt> (*this) </tt>
    NagamochiIbaraki &heap(Heap& hp, HeapCrossRef &cr)
    {
      if (local_heap_cross_ref) {
	delete _heap_cross_ref;
	local_heap_cross_ref=false;
      }
      _heap_cross_ref = &cr;
      if (local_heap) {
	delete _heap;
	local_heap=false;
      }
      _heap = &hp;
      return *this;
    }

    /// \brief Sets the aux graph.
    ///
    /// Sets the aux graph used by algorithm.
    /// If you don't use this function before calling \ref run(),
    /// it will allocate one. The destuctor deallocates this
    /// automatically allocated graph, of course.
    /// \return <tt> (*this) </tt>
    NagamochiIbaraki &auxGraph(AuxGraph& aux_graph)
    {
      if(local_aux_graph) {
	delete _aux_graph;
	local_aux_graph=false;
      }
      _aux_graph = &aux_graph;
      return *this;
    }

    /// \brief Sets the aux capacity map.
    ///
    /// Sets the aux capacity map used by algorithm.
    /// If you don't use this function before calling \ref run(),
    /// it will allocate one. The destuctor deallocates this
    /// automatically allocated graph, of course.
    /// \return <tt> (*this) </tt>
    NagamochiIbaraki &auxCapacityMap(AuxCapacityMap& aux_capacity_map)
    {
      if(local_aux_capacity) {
	delete _aux_capacity;
	local_aux_capacity=false;
      }
      _aux_capacity = &aux_capacity_map;
      return *this;
    }

    /// \name Execution control
    /// The simplest way to execute the algorithm is to use
    /// one of the member functions called \c run().
    /// \n
    /// If you need more control on the execution,
    /// first you must call \ref init() and then call the start()
    /// or proper times the processNextPhase() member functions.

    ///@{

    /// \brief Initializes the internal data structures.
    ///
    /// Initializes the internal data structures.
    void init() {
      _node_num = countNodes(*_graph);
      createStructures();
      createAuxGraph();
    }

  private:

    struct EdgeInfo {
      typename AuxGraph::Node source, target;
      Value capacity;
    };
    
    struct EdgeInfoLess {
      bool operator()(const EdgeInfo& left, const EdgeInfo& right) const {
        return (left.source < right.source) || 
          (left.source == right.source && left.target < right.target);
      }
    };

  public:


    /// \brief Processes the next phase
    ///
    /// Processes the next phase in the algorithm. The function should
    /// be called at most countNodes(graph) - 1 times to get surely
    /// the min cut in the graph.
    ///
    ///\return %True when the algorithm finished.
    bool processNextPhase() {
      if (_node_num <= 1) return true;

      typedef typename AuxGraph::Node Node;
      typedef typename AuxGraph::NodeIt NodeIt;
      typedef typename AuxGraph::UEdge UEdge;
      typedef typename AuxGraph::UEdgeIt UEdgeIt;
      typedef typename AuxGraph::IncEdgeIt IncEdgeIt;
      
      typename AuxGraph::template UEdgeMap<Value> _edge_cut(*_aux_graph);


      for (NodeIt n(*_aux_graph); n != INVALID; ++n) {
        _heap_cross_ref->set(n, Heap::PRE_HEAP);
      }

      std::vector<Node> nodes;
      nodes.reserve(_node_num);
      int sep = 0;

      Value alpha = 0;
      Value emc = std::numeric_limits<Value>::max();

      _heap->push(typename AuxGraph::NodeIt(*_aux_graph), 0);
      while (!_heap->empty()) {
        Node node = _heap->top();
        Value value = _heap->prio();
        
        _heap->pop();
        for (typename AuxGraph::IncEdgeIt e(*_aux_graph, node); 
             e != INVALID; ++e) {
          Node tn = _aux_graph->runningNode(e);
          switch (_heap->state(tn)) {
          case Heap::PRE_HEAP:
            _heap->push(tn, (*_aux_capacity)[e]);
            _edge_cut[e] = (*_heap)[tn];
            break;
          case Heap::IN_HEAP:
            _heap->decrease(tn, (*_aux_capacity)[e] + (*_heap)[tn]);
            _edge_cut[e] = (*_heap)[tn];
            break;
          case Heap::POST_HEAP:
            break;
          }
        }

        alpha += (*_aux_cut_value)[node];
        alpha -= 2 * value;

        nodes.push_back(node);
        if (!_heap->empty()) {
          if (alpha < emc) {
            emc = alpha;
            sep = nodes.size();
          }
        }
      }

      if (int(nodes.size()) < _node_num) {
        _cut.clear();
        for (int i = 0; i < int(nodes.size()); ++i) {
          typename Graph::Node n = (*_first)[nodes[i]];
          while (n != INVALID) {
            _cut.push_back(n);
            n = (*_next)[n];
          }
        }
	_aux_graph->clear();
        _node_num = 1;
        _min_cut = 0;
        return true;
      }

      if (emc < _min_cut) {
        _cut.clear();
        for (int i = 0; i < sep; ++i) {
          typename Graph::Node n = (*_first)[nodes[i]];
          while (n != INVALID) {
            _cut.push_back(n);
            n = (*_next)[n];
          }
        }
        _min_cut = emc;
      }

      typedef typename AuxGraph::template NodeMap<int> UfeCr;
      UfeCr ufecr(*_aux_graph);
      typedef UnionFindEnum<UfeCr> Ufe; 
      Ufe ufe(ufecr);

      for (typename AuxGraph::NodeIt n(*_aux_graph); n != INVALID; ++n) {
        ufe.insert(n);
      }

      for (typename AuxGraph::UEdgeIt e(*_aux_graph); e != INVALID; ++e) {
        if (_edge_cut[e] >= emc) {
          ufe.join(_aux_graph->source(e), _aux_graph->target(e));
        }
      }

      typedef typename Ufe::ClassIt UfeCIt;
      if (ufe.size(UfeCIt(ufe)) == _node_num) {
        _aux_graph->clear();
        _node_num = 1;
        return true;
      }
      
      std::vector<typename AuxGraph::Node> remnodes;

      typename AuxGraph::template NodeMap<UEdge> edges(*_aux_graph, INVALID);
      for (typename Ufe::ClassIt c(ufe); c != INVALID; ++c) {
        if (ufe.size(c) == 1) continue;
	Node cn = ufe.item(c);
        for (typename Ufe::ItemIt r(ufe, c); r != INVALID; ++r) {
          if (static_cast<Node>(r) == static_cast<Node>(cn)) continue;
          _next->set((*_last)[cn], (*_first)[r]);
          _last->set(cn, (*_last)[r]);
          remnodes.push_back(r);
          --_node_num;
        }
      }

      std::vector<EdgeInfo> addedges;
      std::vector<UEdge> remedges;

      for (typename AuxGraph::UEdgeIt e(*_aux_graph);
           e != INVALID; ++e) {
	int sc = ufe.find(_aux_graph->source(e));
	int tc = ufe.find(_aux_graph->target(e));
        if ((ufe.size(sc) == 1 && ufe.size(tc) == 1)) {
          continue;
        }
        if (sc == tc) {
          remedges.push_back(e);
          continue;
        }
        Node sn = ufe.item(sc);
        Node tn = ufe.item(tc);

        EdgeInfo info;
        if (sn < tn) {
          info.source = sn;
          info.target = tn;
        } else {
          info.source = tn;
          info.target = sn;
        }
        info.capacity = (*_aux_capacity)[e];
        addedges.push_back(info);
        remedges.push_back(e);
      }

      for (int i = 0; i < int(remedges.size()); ++i) {
        _aux_graph->erase(remedges[i]);
      }

      sort(addedges.begin(), addedges.end(), EdgeInfoLess());

      {
        int i = 0;
        while (i < int(addedges.size())) {
          Node sn = addedges[i].source;
          Node tn = addedges[i].target;
          Value ec = addedges[i].capacity;
          ++i;
          while (i < int(addedges.size()) && 
                 sn == addedges[i].source && tn == addedges[i].target) {
            ec += addedges[i].capacity;
            ++i;
          }
          typename AuxGraph::UEdge ne = _aux_graph->addEdge(sn, tn);
          (*_aux_capacity)[ne] = ec;
        }
      }

      for (typename Ufe::ClassIt c(ufe); c != INVALID; ++c) {
        if (ufe.size(c) == 1) continue;
	Node cn = ufe.item(c);
        Value cutvalue = 0;
        for (typename AuxGraph::IncEdgeIt e(*_aux_graph, cn);
             e != INVALID; ++e) {
          cutvalue += (*_aux_capacity)[e];
        }
        
        (*_aux_cut_value)[cn] = cutvalue;
        
      }

      for (int i = 0; i < int(remnodes.size()); ++i) {
        _aux_graph->erase(remnodes[i]);
      }

      return _node_num == 1;
    }

    /// \brief Executes the algorithm.
    ///
    /// Executes the algorithm.
    ///
    /// \pre init() must be called
    void start() {
      while (!processNextPhase());
    }


    /// \brief Runs %NagamochiIbaraki algorithm.
    ///
    /// This method runs the %Min cut algorithm
    ///
    /// \note mc.run(s) is just a shortcut of the following code.
    ///\code
    ///  mc.init();
    ///  mc.start();
    ///\endcode
    void run() {
      init();
      start();
    }

    ///@}

    /// \name Query Functions 
    ///
    /// The result of the %NagamochiIbaraki
    /// algorithm can be obtained using these functions.\n 
    /// Before the use of these functions, either run() or start()
    /// must be called.
    
    ///@{

    /// \brief Returns the min cut value.
    ///
    /// Returns the min cut value if the algorithm finished.
    /// After the first processNextPhase() it is a value of a
    /// valid cut in the graph.
    Value minCut() const {
      return _min_cut;
    }

    /// \brief Returns a min cut in a NodeMap.
    ///
    /// It sets the nodes of one of the two partitions to true in
    /// the given BoolNodeMap. The map contains a valid cut if the
    /// map have been set false previously. 
    template <typename NodeMap>
    Value quickMinCut(NodeMap& nodeMap) const { 
      for (int i = 0; i < int(_cut.size()); ++i) {
        nodeMap.set(_cut[i], true);
      }
      return minCut();
    }

    /// \brief Returns a min cut in a NodeMap.
    ///
    /// It sets the nodes of one of the two partitions to true and
    /// the other partition to false. The function first set all of the
    /// nodes to false and after it call the quickMinCut() member.
    template <typename NodeMap>
    Value minCut(NodeMap& nodeMap) const { 
      for (typename Graph::NodeIt it(*_graph); it != INVALID; ++it) {
        nodeMap.set(it, false);      
      }
      quickMinCut(nodeMap);
      return minCut();
    }

    /// \brief Returns a min cut in an EdgeMap.
    ///
    /// If an undirected edge is in a min cut then it will be
    /// set to true and the others will be set to false in the given map.
    template <typename EdgeMap>
    Value cutEdges(EdgeMap& edgeMap) const {
      typename Graph::template NodeMap<bool> cut(*_graph, false);
      quickMinCut(cut);
      for (typename Graph::EdgeIt it(*_graph); it != INVALID; ++it) {
        edgeMap.set(it, cut[_graph->source(it)] ^ cut[_graph->target(it)]);
      }
      return minCut();
    }

    ///@}

  private:


  };
}

#endif
