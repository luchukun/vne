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

#ifndef LEMON_PRIM_H
#define LEMON_PRIM_H

///\ingroup spantree
///\file
///\brief Prim algorithm to compute minimum spanning tree.

#include <lemon/list_graph.h>
#include <lemon/bin_heap.h>
#include <lemon/bits/invalid.h>
#include <lemon/error.h>
#include <lemon/maps.h>
#include <lemon/bits/traits.h>

#include <lemon/concepts/ugraph.h>

namespace lemon {

  ///Default traits class of Prim class.

  ///Default traits class of Prim class.
  ///\param GR Graph type.
  ///\param CM Type of cost map.
  template<class GR, class CM>
  struct PrimDefaultTraits{
    ///The graph type the algorithm runs on. 
    typedef GR UGraph;
    ///The type of the map that stores the edge costs.

    ///The type of the map that stores the edge costs.
    ///It must meet the \ref concepts::ReadMap "ReadMap" concept.
    typedef CM CostMap;
    //The type of the cost of the edges.
    typedef typename CM::Value Value;
    /// The cross reference type used by heap.

    /// The cross reference type used by heap.
    /// Usually it is \c UGraph::NodeMap<int>.
    typedef typename UGraph::template NodeMap<int> HeapCrossRef;
    ///Instantiates a HeapCrossRef.

    ///This function instantiates a \ref HeapCrossRef. 
    /// \param _graph is the graph, to which we would like to define the 
    /// HeapCrossRef.
    static HeapCrossRef *createHeapCrossRef(const GR &_graph){
      return new HeapCrossRef(_graph);
    }
    
    ///The heap type used by Prim algorithm.

    ///The heap type used by Prim algorithm.
    ///
    ///\sa BinHeap
    ///\sa Prim
    typedef BinHeap<typename CM::Value,
		    HeapCrossRef, std::less<Value> > Heap;

    static Heap *createHeap(HeapCrossRef& _ref){
      return new Heap(_ref);
    }

    ///\brief The type of the map that stores the last
    ///edges of the minimum spanning tree.
    /// 
    ///The type of the map that stores the last
    ///edges of the minimum spanning tree.
    ///It must meet the \ref concepts::WriteMap "WriteMap" concept.
    ///
    typedef typename UGraph::template NodeMap<typename GR::UEdge> PredMap;
    ///Instantiates a PredMap.
 
    ///This function instantiates a \ref PredMap. 
    ///\param _graph is the graph, to which we would like to define the PredMap.
    static PredMap *createPredMap(const GR &_graph){
      return new PredMap(_graph);
    }

    ///The type of the map that stores whether an edge is in the
    ///spanning tree or not.

    ///The type of the map that stores whether an edge is in the
    ///spanning tree or not.
    ///By default it is a NullMap.
    typedef NullMap<typename UGraph::UEdge,bool> TreeMap;
    ///Instantiates a TreeMap.

    ///This function instantiates a \ref TreeMap.
    ///
    ///The first parameter is the graph, to which
    ///we would like to define the \ref TreeMap
    static TreeMap *createTreeMap(const GR &){
      return new TreeMap();
    }

    ///The type of the map that stores whether a nodes is processed.
 
    ///The type of the map that stores whether a nodes is processed.
    ///It must meet the \ref concepts::WriteMap "WriteMap" concept.
    ///By default it is a NodeMap<bool>.
    typedef NullMap<typename UGraph::Node,bool> ProcessedMap;
    ///Instantiates a ProcessedMap.
 
    ///This function instantiates a \ref ProcessedMap. 
    ///\param _graph is the graph, to which
    ///we would like to define the \ref ProcessedMap
#ifdef DOXYGEN
    static ProcessedMap *createProcessedMap(const GR &_graph)
#else
    static ProcessedMap *createProcessedMap(const GR &)
#endif
    {
      return new ProcessedMap();
    }
  };
  
  ///%Prim algorithm class to find a minimum spanning tree.
  
  /// \ingroup spantree
  ///This class provides an efficient implementation of %Prim algorithm.
  ///
  ///The running time is \f$ O(e\log(n)) \f$ where e is the number of edges and
  ///n is the number of nodes in the graph.
  ///
  ///The edge costs are passed to the algorithm using a
  ///\ref concepts::ReadMap "ReadMap",
  ///so it is easy to change it to any kind of cost.
  ///
  ///The type of the cost is determined by the
  ///\ref concepts::ReadMap::Value "Value" of the cost map.
  ///
  ///It is also possible to change the underlying priority heap.
  ///
  ///\param GR The graph type the algorithm runs on. The default value
  ///is \ref ListUGraph. The value of GR is not used directly by
  ///Prim, it is only passed to \ref PrimDefaultTraits.
  ///
  ///\param CM This read-only UEdgeMap determines the costs of the
  ///edges. It is read once for each edge, so the map may involve in
  ///relatively time consuming process to compute the edge cost if
  ///it is necessary. The default map type is \ref
  ///concepts::UGraph::UEdgeMap "UGraph::UEdgeMap<int>".  The value
  ///of CM is not used directly by Prim, it is only passed to \ref
  ///PrimDefaultTraits.
  ///
  ///\param TR Traits class to set
  ///various data types used by the algorithm.  The default traits
  ///class is \ref PrimDefaultTraits
  ///"PrimDefaultTraits<GR,CM>".  See \ref
  ///PrimDefaultTraits for the documentation of a Prim traits
  ///class.
  ///
  ///\author Balazs Attila Mihaly

#ifdef DOXYGEN
  template <typename GR,
	    typename CM,
	    typename TR>
#else
  template <typename GR=ListUGraph,
	    typename CM=typename GR::template UEdgeMap<int>,
	    typename TR=PrimDefaultTraits<GR,CM> >
#endif
  class Prim {
  public:
    
    /// \brief \ref Exception for uninitialized parameters.
    ///
    /// This error represents problems in the initialization
    /// of the parameters of the algorithms.
    class UninitializedParameter : public lemon::UninitializedParameter {
    public:
      virtual const char* what() const throw() {
	return "lemon::Prim::UninitializedParameter";
      }
    };

    typedef TR Traits;
    ///The type of the underlying graph.
    typedef typename TR::UGraph UGraph;
    ///\e
    typedef typename UGraph::Node Node;
    ///\e
    typedef typename UGraph::NodeIt NodeIt;
    ///\e
    typedef typename UGraph::UEdge UEdge;
    ///\e
    typedef typename UGraph::IncEdgeIt IncEdgeIt;
    
    ///The type of the cost of the edges.
    typedef typename TR::CostMap::Value Value;
    ///The type of the map that stores the edge costs.
    typedef typename TR::CostMap CostMap;
    ///\brief The type of the map that stores the last
    ///predecessor edges of the spanning tree.
    typedef typename TR::PredMap PredMap;
    ///Edges of the spanning tree.
    typedef typename TR::TreeMap TreeMap;
    ///The type of the map indicating if a node is processed.
    typedef typename TR::ProcessedMap ProcessedMap;
    ///The cross reference type used for the current heap.
    typedef typename TR::HeapCrossRef HeapCrossRef;
    ///The heap type used by the prim algorithm.
    typedef typename TR::Heap Heap;
  private:
    /// Pointer to the underlying graph.
    const UGraph *graph;
    /// Pointer to the cost map
    const CostMap *cost;
    ///Pointer to the map of predecessors edges.
    PredMap *_pred;
    ///Indicates if \ref _pred is locally allocated (\c true) or not.
    bool local_pred;
    ///Pointer to the map of tree edges.
    TreeMap *_tree;
    ///Indicates if \ref _tree is locally allocated (\c true) or not.
    bool local_tree;
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

    ///Creates the maps if necessary.
    void create_maps(){
      if(!_pred) {
	local_pred = true;
	_pred = Traits::createPredMap(*graph);
      }
      if(!_tree) {
	local_tree = true;
	_tree = Traits::createTreeMap(*graph);
      }
      if(!_processed) {
	local_processed = true;
	_processed = Traits::createProcessedMap(*graph);
      }
      if (!_heap_cross_ref) {
	local_heap_cross_ref = true;
	_heap_cross_ref = Traits::createHeapCrossRef(*graph);
      }
      if (!_heap) {
	local_heap = true;
	_heap = Traits::createHeap(*_heap_cross_ref);
      }
    }
    
  public :

    typedef Prim Create;
 
    ///\name Named template parameters

    ///@{

    template <class T>
    struct DefPredMapTraits : public Traits {
      typedef T PredMap;
      static PredMap *createPredMap(const UGraph &_graph){
	throw UninitializedParameter();
      }
    };
    ///\brief \ref named-templ-param "Named parameter" for setting PredMap type
    ///
    ///\ref named-templ-param "Named parameter" for setting PredMap type
    ///
    template <class T>
    struct DefPredMap 
      : public Prim< UGraph, CostMap, DefPredMapTraits<T> > {
      typedef Prim< UGraph, CostMap, DefPredMapTraits<T> > Create;
    };
    
    template <class T>
    struct DefProcessedMapTraits : public Traits {
      typedef T ProcessedMap;
      static ProcessedMap *createProcessedMap(const UGraph &_graph){
	throw UninitializedParameter();
      }
    };
    ///\brief \ref named-templ-param "Named parameter" for setting
    ///ProcessedMap type
    ///
    ///\ref named-templ-param "Named parameter" for setting ProcessedMap type
    ///
    template <class T>
    struct DefProcessedMap 
      : public Prim< UGraph, CostMap, DefProcessedMapTraits<T> > { 
      typedef Prim< UGraph, CostMap, DefProcessedMapTraits<T> > Create;
    };
    
    struct DefGraphProcessedMapTraits : public Traits {
      typedef typename UGraph::template NodeMap<bool> ProcessedMap;
      static ProcessedMap *createProcessedMap(const UGraph &_graph){
	return new ProcessedMap(_graph);
      }
    };


    template <class H, class CR>
    struct DefHeapTraits : public Traits {
      typedef CR HeapCrossRef;
      typedef H Heap;
      static HeapCrossRef *createHeapCrossRef(const UGraph &) {
	throw UninitializedParameter();
      }
      static Heap *createHeap(HeapCrossRef &){
	return UninitializedParameter();
      }
    };
    ///\brief \ref named-templ-param "Named parameter" for setting
    ///heap and cross reference type
    ///
    ///\ref named-templ-param "Named parameter" for setting heap and cross 
    ///reference type
    ///
    template <class H, class CR = typename UGraph::template NodeMap<int> >
    struct DefHeap
      : public Prim< UGraph, CostMap, DefHeapTraits<H, CR> > {
      typedef Prim< UGraph, CostMap, DefHeapTraits<H, CR> > Create;
    };

    template <class H, class CR>
    struct DefStandardHeapTraits : public Traits {
      typedef CR HeapCrossRef;
      typedef H Heap;
      static HeapCrossRef *createHeapCrossRef(const UGraph &_graph) {
	return new HeapCrossRef(_graph);
      }
      static Heap *createHeap(HeapCrossRef &ref){
	return new Heap(ref);
      }
    };
    ///\brief \ref named-templ-param "Named parameter" for setting
    ///heap and cross reference type with automatic allocation
    ///
    ///\ref named-templ-param "Named parameter" for setting heap and cross 
    ///reference type. It can allocate the heap and the cross reference 
    ///object if the cross reference's constructor waits for the graph as 
    ///parameter and the heap's constructor waits for the cross reference.
    template <class H, class CR = typename UGraph::template NodeMap<int> >
    struct DefStandardHeap
      : public Prim< UGraph, CostMap, DefStandardHeapTraits<H, CR> > { 
      typedef Prim< UGraph, CostMap, DefStandardHeapTraits<H, CR> > 
      Create;
    };

    template <class TM>
    struct DefTreeMapTraits : public Traits {
      typedef TM TreeMap;
      static TreeMap *createTreeMap(const UGraph &) {
        throw UninitializedParameter();
      }
    };
    ///\brief \ref named-templ-param "Named parameter" for setting
    ///TreeMap
    ///
    ///\ref named-templ-param "Named parameter" for setting TreeMap
    ///
    template <class TM>
    struct DefTreeMap
      : public Prim< UGraph, CostMap, DefTreeMapTraits<TM> > {
      typedef Prim< UGraph, CostMap, DefTreeMapTraits<TM> > Create;
    };    

    struct DefGraphTreeMapTraits : public Traits {
      typedef typename UGraph::template NodeMap<bool> TreeMap;
      static TreeMap *createTreeMap(const UGraph &_graph){
	return new TreeMap(_graph);
      }
    };

    ///@}


  protected:

    Prim() {}

  public:      
    
    ///Constructor.
    ///
    ///\param _graph the graph the algorithm will run on.
    ///\param _cost the cost map used by the algorithm.
    Prim(const UGraph& _graph, const CostMap& _cost) :
      graph(&_graph), cost(&_cost),
      _pred(0), local_pred(false),
      _tree(0), local_tree(false),
      _processed(0), local_processed(false),
      _heap_cross_ref(0), local_heap_cross_ref(false),
      _heap(0), local_heap(false)
    {
      checkConcept<concepts::UGraph, UGraph>();
    }
    
    ///Destructor.
    ~Prim(){
      if(local_pred) delete _pred;
      if(local_tree) delete _tree;
      if(local_processed) delete _processed;
      if(local_heap_cross_ref) delete _heap_cross_ref;
      if(local_heap) delete _heap;
    }

    ///\brief Sets the cost map.
    ///
    ///Sets the cost map.
    ///\return <tt> (*this) </tt>
    Prim &costMap(const CostMap &m){
      cost = &m;
      return *this;
    }

    ///\brief Sets the map storing the predecessor edges.
    ///
    ///Sets the map storing the predecessor edges.
    ///If you don't use this function before calling \ref run(),
    ///it will allocate one. The destuctor deallocates this
    ///automatically allocated map, of course.
    ///\return <tt> (*this) </tt>
    Prim &predMap(PredMap &m){
      if(local_pred) {
	delete _pred;
	local_pred=false;
      }
      _pred = &m;
      return *this;
    }

    ///\brief Sets the map storing the tree edges.
    ///
    ///Sets the map storing the tree edges.
    ///If you don't use this function before calling \ref run(),
    ///it will allocate one. The destuctor deallocates this
    ///automatically allocated map, of course.
    ///By default this is a NullMap.
    ///\return <tt> (*this) </tt>
    Prim &treeMap(TreeMap &m){
      if(local_tree) {
	delete _tree;
	local_tree=false;
      }
      _tree = &m;
      return *this;
    }

    ///\brief Sets the heap and the cross reference used by algorithm.
    ///
    ///Sets the heap and the cross reference used by algorithm.
    ///If you don't use this function before calling \ref run(),
    ///it will allocate one. The destuctor deallocates this
    ///automatically allocated map, of course.
    ///\return <tt> (*this) </tt>
    Prim &heap(Heap& heap, HeapCrossRef &crossRef){
      if(local_heap_cross_ref) {
	delete _heap_cross_ref;
	local_heap_cross_ref=false;
      }
      _heap_cross_ref = &crossRef;
      if(local_heap) {
	delete _heap;
	local_heap=false;
      }
      _heap = &heap;
      return *this;
    }

  public:
    ///\name Execution control
    ///The simplest way to execute the algorithm is to use
    ///one of the member functions called \c run(...).
    ///\n
    ///If you need more control on the execution,
    ///first you must call \ref init(), then you can add several source nodes
    ///with \ref addSource().
    ///Finally \ref start() will perform the actual path
    ///computation.

    ///@{

    ///\brief Initializes the internal data structures.
    ///
    ///Initializes the internal data structures.
    ///
    void init(){
      create_maps();
      _heap->clear();
      for ( NodeIt u(*graph) ; u!=INVALID ; ++u ) {
	_pred->set(u,INVALID);
	_processed->set(u,false);
	_heap_cross_ref->set(u,Heap::PRE_HEAP);
      }
    }
    
    ///\brief Adds a new source node.
    ///
    ///Adds a new source node to the priority heap.
    ///
    ///It checks if the node has already been added to the heap and
    ///it is pushed to the heap only if it was not in the heap.
    void addSource(Node s){
      if(_heap->state(s) != Heap::IN_HEAP) {
	_heap->push(s,Value());
      }
    }
    ///\brief Processes the next node in the priority heap
    ///
    ///Processes the next node in the priority heap.
    ///
    ///\return The processed node.
    ///
    ///\warning The priority heap must not be empty!
    Node processNextNode(){
      Node v=_heap->top(); 
      _heap->pop();
      _processed->set(v,true);
      
      for(IncEdgeIt e(*graph,v); e!=INVALID; ++e) {
	Node w=graph->oppositeNode(v,e);
	switch(_heap->state(w)) {
	case Heap::PRE_HEAP:
	  _heap->push(w,(*cost)[e]);
	  _pred->set(w,e);
	  break;
	case Heap::IN_HEAP:
	  if ( (*cost)[e] < (*_heap)[w] ) {
	    _heap->decrease(w,(*cost)[e]); 
	    _pred->set(w,e);
	  }
	  break;
	case Heap::POST_HEAP:
	  break;
	}
      }
      if ((*_pred)[v]!=INVALID)_tree->set((*_pred)[v],true);
      return v;
    }

    ///\brief Next node to be processed.
    ///
    ///Next node to be processed.
    ///
    ///\return The next node to be processed or INVALID if the priority heap
    /// is empty.
    Node nextNode(){ 
      return _heap->empty()?_heap->top():INVALID;
    }
 
    ///\brief Returns \c false if there are nodes to be processed in
    ///the priority heap
    ///
    ///Returns \c false if there are nodes
    ///to be processed in the priority heap
    bool emptyQueue() { return _heap->empty(); }

    ///\brief Returns the number of the nodes to be processed in the
    ///priority heap
    ///
    ///Returns the number of the nodes to be processed in the priority heap
    ///
    int queueSize() { return _heap->size(); }
    
    ///\brief Executes the algorithm.
    ///
    ///Executes the algorithm.
    ///
    ///\pre init() must be called and at least one node should be added
    ///with addSource() before using this function.
    ///
    ///This method runs the %Prim algorithm from the node(s)
    ///in order to compute the
    ///minimum spanning tree.
    ///
    void start(){
      while ( !_heap->empty() ) processNextNode();
    }
    
    ///\brief Executes the algorithm until a condition is met.
    ///
    ///Executes the algorithm until a condition is met.
    ///
    ///\pre init() must be called and at least one node should be added
    ///with addSource() before using this function.
    ///
    ///\param nm must be a bool (or convertible) node map. The algorithm
    ///will stop when it reaches a node \c v with <tt>nm[v]==true</tt>.
    template<class NodeBoolMap>
    void start(const NodeBoolMap &nm){
      while ( !_heap->empty() && !nm[_heap->top()] ) processNextNode();
      if ( !_heap->empty() ) _processed->set(_heap->top(),true);
    }
    
    ///\brief Runs %Prim algorithm.
    ///
    ///This method runs the %Prim algorithm
    ///in order to compute the
    ///minimum spanning tree (or minimum spanning forest).
    ///The method also works on graphs that has more than one components.
    ///In this case it computes the minimum spanning forest.
    void run() {
      init();
      for(NodeIt it(*graph);it!=INVALID;++it){
	if(!processed(it)){
	  addSource(it);
	  start();
	}
      }
    }

    ///\brief Runs %Prim algorithm from node \c s.
    ///
    ///This method runs the %Prim algorithm from node \c s
    ///in order to
    ///compute the
    ///minimun spanning tree
    ///
    ///\note p.run(s) is just a shortcut of the following code.
    ///\code
    ///  p.init();
    ///  p.addSource(s);
    ///  p.start();
    ///\endcode
    ///\note If the graph has more than one components, the method
    ///will compute the minimun spanning tree for only one component.
    ///
    ///See \ref run() if you want to compute the minimal spanning forest.
    void run(Node s){
      init();
      addSource(s);
      start();
    }
    
    ///@}

    ///\name Query Functions
    ///The result of the %Prim algorithm can be obtained using these
    ///functions.\n
    ///Before the use of these functions,
    ///either run() or start() must be called.
    
    ///@{

    ///\brief Returns the 'previous edge' of the minimum spanning tree.

    ///For a node \c v it returns the 'previous edge' of the minimum
    ///spanning tree, i.e. it returns the edge from where \c v was
    ///reached. For a source node or an unreachable node it is \ref
    ///INVALID.  The minimum spanning tree used here is equal to the
    ///minimum spanning tree used in \ref predNode().  
    ///\pre \ref run() or \ref start() must be called before using
    ///this function.
    UEdge predEdge(Node v) const { return (*_pred)[v]; }

    ///\brief Returns the 'previous node' of the minimum spanning
    ///tree.
    ///
    ///For a node \c v it returns the 'previous node' of the minimum
    ///spanning tree, i.e. it returns the node from where \c v was
    ///reached. For a source node or an unreachable node it is \ref
    ///INVALID.  //The minimum spanning tree used here is equal to the
    ///minimum spanning tree used in \ref predEdge().  
    ///\pre \ref run() or \ref start() must be called before using
    ///this function.
    Node predNode(Node v) const { return (*_pred)[v]==INVALID ? INVALID:
				  graph->source((*_pred)[v]); }
    
    ///\brief Returns a reference to the NodeMap of the edges of the
    ///minimum spanning tree.
    ///
    ///Returns a reference to the NodeMap of the edges of the minimum
    ///spanning tree.
    ///\pre \ref run() or \ref start() must be called before using
    ///this function.
    const PredMap &predMap() const { return *_pred;}
 
    ///\brief Returns a reference to the tree edges map.

    ///Returns a reference to the TreeEdgeMap of the edges of the
    ///minimum spanning tree. The value of the map is \c true only if
    ///the edge is in the minimum spanning tree.
    ///
    ///\warning By default, the TreeEdgeMap is a NullMap.
    ///
    ///If it is not set before the execution of the algorithm, use the
    ///\ref treeMap(TreeMap&) function (after the execution) to set an
    ///UEdgeMap with the edges of the minimum spanning tree in O(n)
    ///time where n is the number of nodes in the graph.
    ///\pre \ref run() or \ref start() must be called before using
    ///this function.
    const TreeMap &treeMap() const { return *_tree;}
 
    ///\brief Sets the tree edges map.
    ///
    ///Sets the TreeMap of the edges of the minimum spanning tree.
    ///The map values belonging to the edges of the minimum
    ///spanning tree are set to \c tree_edge_value or \c true by default,
    ///the other map values remain untouched.
    ///
    ///\pre \ref run() or \ref start() must be called before using
    ///this function.

    template<class TreeMap>
    void quickTreeEdges(TreeMap& tree) const {
      for(NodeIt i(*graph);i!=INVALID;++i){
        if((*_pred)[i]!=INVALID) tree.set((*_pred)[i],true);
      }
    }

    ///\brief Sets the tree edges map.
    ///
    ///Sets the TreeMap of the edges of the minimum spanning tree.
    ///The map values belonging to the edges of the minimum
    ///spanning tree are set to \c tree_edge_value or \c true by default while
    ///the edge values not belonging to the minimum spanning tree are set to
    ///\c tree_default_value or \c false by default.
    ///
    ///\pre \ref run() or \ref start() must be called before using
    ///this function.
    template <class TreeMap>
    void treeEdges(TreeMap& tree) const {
      typedef typename ItemSetTraits<UGraph,UEdge>::ItemIt TreeMapIt;
      for(TreeMapIt i(*graph); i != INVALID; ++i) {
	tree.set(i,false);
      }
      for(NodeIt i(*graph); i != INVALID; ++i){
        if((*_pred)[i] != INVALID) tree.set((*_pred)[i],true);
      }
    }

    ///\brief Checks if a node is reachable from the starting node.
    ///
    ///Returns \c true if \c v is reachable from the starting node.
    ///\warning The source nodes are inditated as unreached.
    ///\pre \ref run() or \ref start() must be called before using
    ///this function.
    ///
    bool reached(Node v) { return (*_heap_cross_ref)[v] != Heap::PRE_HEAP; }

    ///\brief Checks if a node is processed.
    ///
    ///Returns \c true if \c v is processed, i.e. \c v is already
    ///connencted to the minimum spanning tree.  \pre \ref run() or
    ///\ref start() must be called before using this function.
    bool processed(Node v) { return (*_heap_cross_ref)[v] == Heap::POST_HEAP; }
    

    ///\brief Checks if an edge is in the spanning tree or not.
    ///
    ///Checks if an edge is in the spanning tree or not.
    ///\param e is the edge that will be checked
    ///\return \c true if e is in the spanning tree, \c false otherwise
    bool tree(UEdge e){
      return (*_pred)[*graph.source(e)]==e || (*_pred)[*graph.target(e)]==e;
    }

    ///\brief Returns the value of the total cost of the spanning tree.
    ///
    /// Returns the value of the total cost of the spanning tree.
    Value treeValue() const {
      Value value = 0;
      for(NodeIt i(*graph); i!= INVALID; ++i){
        if ((*_pred)[i] != INVALID) value += (*cost)[(*_pred)[i]];
      }
      return value;
    }
    ///@}
  };


  /// \ingroup spantree
  ///
  /// \brief Function type interface for Prim algorithm.
  ///
  /// Function type interface for Prim algorithm.
  /// \param graph the UGraph that the algorithm runs on
  /// \param cost the CostMap of the edges
  /// \retval tree the EdgeMap that contains whether an edge is in 
  /// the spanning tree or not
  /// \return The total cost of the spanning tree
  ///
  ///\sa Prim
  template<class Graph,class CostMap,class TreeMap>
  typename CostMap::Value prim(const Graph& graph, 
                               const CostMap& cost,
                               TreeMap& tree){
    typename Prim<Graph,CostMap>::
      template DefTreeMap<TreeMap>::
      Create prm(graph,cost);
    prm.treeMap(tree);
    prm.run();
    return prm.treeValue();
  }

  /// \ingroup spantree
  ///
  /// \brief Function type interface for Prim algorithm.
  ///
  /// Function type interface for Prim algorithm.
  /// \param graph the UGraph that the algorithm runs on
  /// \param cost the CostMap of the edges
  /// the spanning tree or not
  /// \return The total cost of the spanning tree
  ///
  ///\sa Prim
  template<class Graph,class CostMap,class TreeMap>
  typename CostMap::Value prim(const Graph& graph, 
                               const CostMap& cost){
    typename Prim<Graph,CostMap>::
      Create prm(graph,cost);
    prm.run();
    return prm.treeValue();
  }

} //END OF NAMESPACE LEMON

#endif
