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

#ifndef LEMON_FREDMAN_TARJAN_H
#define LEMON_FREDMAN_TARJAN_H

///\ingroup spantree
///\file
///\brief FredmanTarjan algorithm to compute minimum spanning forest.

#include <limits>
#include <vector>

#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/fib_heap.h>
#include <lemon/radix_sort.h>
#include <lemon/bits/invalid.h>
#include <lemon/error.h>
#include <lemon/maps.h>
#include <lemon/bits/traits.h>
#include <lemon/graph_utils.h>

#include <lemon/concepts/ugraph.h>

namespace lemon {

  ///Default traits class of FredmanTarjan class.

  ///Default traits class of FredmanTarjan class.
  ///\param GR Graph type.
  ///\param CM Type of cost map.
  template<class GR, class CM>
  struct FredmanTarjanDefaultTraits{
    ///The graph type the algorithm runs on. 
    typedef GR UGraph;
    ///The type of the map that stores the edge costs.

    ///The type of the map that stores the edge costs.
    ///It must meet the \ref concepts::ReadMap "ReadMap" concept.
    typedef CM CostMap;
    //The type of the cost of the edges.
    typedef typename CM::Value Value;
    ///The type of the map that stores whether an edge is in the
    ///spanning tree or not.

    ///The type of the map that stores whether an edge is in the
    ///spanning tree or not.
    ///It must meet the \ref concepts::ReadWriteMap "ReadWriteMap" concept.
    ///By default it is a BoolEdgeMap.
    typedef typename UGraph::template UEdgeMap<bool> TreeMap;
    ///Instantiates a TreeMap.

    ///This function instantiates a \ref TreeMap.
    ///\param _graph is the graph, to which
    ///we would like to define the \ref TreeMap
    static TreeMap *createTreeMap(const GR &_graph){
      return new TreeMap(_graph);
    }
  };
  
  ///%FredmanTarjan algorithm class to find a minimum spanning tree.
  
  /// \ingroup spantree 
  ///
  ///This class provides an efficient implementation of %FredmanTarjan
  ///algorithm whitch is sometimes a bit quicker than the Prim
  ///algorithm on larger graphs.  Due to the structure of the
  ///algorithm, it has less controll functions than Prim.
  ///
  /// The running time is \f$ O(e\beta(e,n)) \f$ where 
  /// \f$ \beta(e,n) = \min\{ i | \log^{i}(n) \le e/n\} \f$ and 
  /// \f$ \log^{i+1}(n)=\log(\log^{i}(n)) \f$
  ///
  ///The edge costs are passed to the algorithm using a \ref
  ///concepts::ReadMap "ReadMap", so it is easy to change it to any
  ///kind of cost.
  ///
  ///The type of the cost is determined by the \ref
  ///concepts::ReadMap::Value "Value" of the cost map.
  ///
  ///\param GR The graph type the algorithm runs on. The default value
  ///is \ref ListUGraph. The value of GR is not used directly by
  ///FredmanTarjan, it is only passed to \ref
  ///FredmanTarjanDefaultTraits.
  ///
  ///\param CM This read-only UEdgeMap determines the costs of the
  ///edges. It is read once for each edge, so the map may involve in
  ///relatively time consuming process to compute the edge cost if it
  ///is necessary. The default map type is \ref
  ///concepts::UGraph::UEdgeMap "UGraph::UEdgeMap<int>". The value of
  ///CM is not used directly by FredmanTarjan, it is only passed to
  ///\ref FredmanTarjanDefaultTraits.
  ///
  ///\param TR Traits class to set various data types used by the
  ///algorithm.  The default traits class is \ref
  ///FredmanTarjanDefaultTraits "FredmanTarjanDefaultTraits<GR,CM>".
  ///See \ref FredmanTarjanDefaultTraits for the documentation of a
  ///FredmanTarjan traits class.
  ///
  ///\author Balazs Attila Mihaly

#ifdef DOXYGEN
  template <typename GR,
	    typename CM,
	    typename TR>
#else
  template <typename GR=ListUGraph,
	    typename CM=typename GR::template UEdgeMap<int>,
	    typename TR=FredmanTarjanDefaultTraits<GR,CM> >
#endif
  class FredmanTarjan {
  public:
    ///\brief \ref Exception for uninitialized parameters.
    ///
    ///This error represents problems in the initialization
    ///of the parameters of the algorithms.
    class UninitializedParameter : public lemon::UninitializedParameter {
    public:
      virtual const char* what() const throw() {
	return "lemon::FredmanTarjan::UninitializedParameter";
      }
    };

    typedef GR Graph;
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
    typedef typename UGraph::UEdgeIt UEdgeIt;
    ///\e
    typedef typename UGraph::IncEdgeIt IncEdgeIt;
    
    ///The type of the cost of the edges.
    typedef typename TR::CostMap::Value Value;
    ///The type of the map that stores the edge costs.
    typedef typename TR::CostMap CostMap;
    ///Edges of the spanning tree.
    typedef typename TR::TreeMap TreeMap;
  private:
    ///Pointer to the underlying graph.
    const UGraph *graph;
    ///Pointer to the cost map
    const CostMap *cost;
    ///Pointer to the map of tree edges.
    TreeMap *_tree;
    ///Indicates if \ref _tree is locally allocated (\c true) or not.
    bool local_tree;

    ///Creates the maps if necessary.
    
    void create_maps(){
      if(!_tree){
	local_tree=true;
	_tree=Traits::createTreeMap(*graph);
      }
    }
    
  public :

    typedef FredmanTarjan Create;
 
    ///\name Named template parameters

    ///@{

    template <class TM>
    struct DefTreeMapTraits : public Traits {
      typedef TM TreeMap;
      static TreeMap *createTreeMap(const UGraph &) {
	throw UninitializedParameter();
      }
    };
    ///\brief \ref named-templ-param "Named parameter" for setting TreeMap 
    ///
    ///\ref named-templ-param "Named parameter" for setting TreeMap
    ///
    template <class TM>
    struct DefTreeMap
      : public FredmanTarjan< UGraph, CostMap, DefTreeMapTraits<TM> > { 
      typedef FredmanTarjan< UGraph, CostMap, DefTreeMapTraits<TM> > Create;
    };

    ///@}


  protected:

    FredmanTarjan() {}

  private:

    template<class SrcGraph,class OrigMap,class Heap,class HeapCrossRef,
             class ProcessedMap,class PredMap>
    void processNextTree(const SrcGraph& graph,const OrigMap& orig,
                         Heap &heap, HeapCrossRef& crossref,
                         ProcessedMap& processed,PredMap& pred,
                         int& tree_counter,const int limit){
      std::vector<typename SrcGraph::Node> tree_nodes;
      int tree_index=tree_counter;
      bool stop=false;
      while(!heap.empty() && !stop){
        typename SrcGraph::Node v=heap.top();
        heap.pop();
	if(processed[v]!=-1){
	  heap.state(v,Heap::PRE_HEAP);
	  tree_index=processed[v];
	  _tree->set(orig[pred[v]],true);
	  stop=true;
	  break;
        }
	tree_nodes.push_back(v);
	for(typename SrcGraph::IncEdgeIt e(graph,v);e!=INVALID;++e){
	  typename SrcGraph::Node w=graph.oppositeNode(v,e);
	  switch(heap.state(w)){
	  case Heap::PRE_HEAP:
	    if(heap.size()>=limit){
	      stop=true;
	    }
	    else{
	      heap.push(w,(*cost)[orig[e]]);
	      pred.set(w,e);
	    }
	    break;
	  case Heap::IN_HEAP:
	    if ((*cost)[orig[e]]<heap[w]){
	      heap.decrease(w,(*cost)[orig[e]]); 
	      pred.set(w,e);
	    }
	    break;
	  case Heap::POST_HEAP:
	    break;
	  }
	}
      }
      for(int i=1;i<(int)tree_nodes.size();++i){
	_tree->set(orig[pred[tree_nodes[i]]],true);
        processed.set(tree_nodes[i],tree_index);
        crossref[tree_nodes[i]] = Heap::PRE_HEAP;
      }
      processed.set(tree_nodes[0],tree_index);
      crossref[tree_nodes[0]] = Heap::PRE_HEAP;
      while (!heap.empty()) {
        typename SrcGraph::Node v=heap.top();
	heap.pop();
        crossref[v] = Heap::PRE_HEAP;
      }
      heap.clear();
      if(!stop)++tree_counter;
    }

    template<class SrcGraph,class OrigMap,class ProcessedMap>
    void createTrees(const SrcGraph& graph, const OrigMap& orig, 
                     ProcessedMap& processed,
                     int edgenum,int& tree_counter){
      typedef typename SrcGraph::Node Node;
      typedef typename SrcGraph::UEdge UEdge;
      typedef typename SrcGraph::NodeIt NodeIt;
      typedef typename SrcGraph::template NodeMap<int> HeapCrossRef;
      typedef typename SrcGraph::template NodeMap<UEdge> PredMap;
      HeapCrossRef crossref(graph,-1);
      FibHeap<Value,HeapCrossRef> heap(crossref);
      PredMap pred(graph,INVALID);
      int rate=2*edgenum/countNodes(graph);
      int limit=(rate>std::numeric_limits<int>::digits)?
      std::numeric_limits<int>::max() : (1<<rate);
      for(NodeIt i(graph);i!=INVALID;++i){
	if(processed[i]==-1){
	  heap.push(i, Value());
	  processNextTree(graph,orig,heap,crossref,
                          processed,pred,tree_counter,limit);
	}
      }
    }

    template<class SrcGraph,class DestGraph,class SrcOrigMap,
             class DestOrigMap,class ProcessedMap>
    void collect(const SrcGraph& srcgraph,const SrcOrigMap& srcorig,
                 DestGraph& destgraph,DestOrigMap& destorig,
                 const ProcessedMap& processed,const int tree_counter){
      typedef typename SrcGraph::Node Node;
      typedef typename DestGraph::Node DNode;
      typedef typename SrcGraph::UEdge UEdge;
      typedef typename DestGraph::UEdge DUEdge;
      typedef typename SrcGraph::Edge Edge;
      typedef typename SrcGraph::EdgeIt EdgeIt;
      std::vector<Edge> edges;
      std::vector<DNode> nodes(tree_counter, INVALID);
      for(EdgeIt i(srcgraph);i!=INVALID;++i){
	if(processed[srcgraph.source(i)]<processed[srcgraph.target(i)]){
	  edges.push_back(i);
          if(nodes[processed[srcgraph.source(i)]]==INVALID) {
	    nodes[processed[srcgraph.source(i)]]=destgraph.addNode();
	  }
          if(nodes[processed[srcgraph.target(i)]]==INVALID) {
	    nodes[processed[srcgraph.target(i)]]=destgraph.addNode();
	  }
	}
      }
      
      radixSort(edges.begin(),edges.end(),
                mapFunctor(composeMap(processed,sourceMap(srcgraph))));
      counterSort(edges.begin(),edges.end(),
                  mapFunctor(composeMap(processed,targetMap(srcgraph))));
      for(int i=0;i!=(int)edges.size();++i){
	int srcproc=processed[srcgraph.source(edges[i])];
	int trgproc=processed[srcgraph.target(edges[i])];
        Value minval=(*cost)[srcorig[edges[i]]];
        UEdge minpos=edges[i];
	while (i+1!=(int)edges.size() && 
               srcproc==processed[srcgraph.source(edges[i+1])] &&
	  trgproc==processed[srcgraph.target(edges[i+1])]) {
	  if (minval>(*cost)[srcorig[edges[i+1]]]) {
            minval=(*cost)[srcorig[edges[i+1]]];
            minpos=edges[i+1];
	  }
          ++i;
	} 
	destorig[destgraph.addEdge(nodes[srcproc],nodes[trgproc])]=
          srcorig[minpos];
      }
    }

    template<class SrcGraph,class OrigMap>
    void phase(const SrcGraph& graph,const OrigMap& orig,int edgenum){
      int tree_counter = 0;
      typename SrcGraph::template NodeMap<int> processed(graph,-1);
      SmartUGraph destgraph;
      SmartUGraph::UEdgeMap<typename OrigMap::Value> destorig(destgraph);
      createTrees(graph,orig,processed,edgenum,tree_counter);
      collect(graph,orig,destgraph,destorig,processed,tree_counter);
      if (countNodes(destgraph)>1) {
        phase(destgraph,destorig,edgenum);
      }
    }

  public:      
    
    ///Constructor.
    
    ///\param _graph the graph the algorithm will run on.
    ///\param _cost the cost map used by the algorithm.
    FredmanTarjan(const UGraph& _graph, const CostMap& _cost) :
      graph(&_graph), cost(&_cost),
      _tree(0), local_tree(false)
    {
      checkConcept<concepts::UGraph, UGraph>();
    }
    
    ///Destructor.
    ~FredmanTarjan(){
      if(local_tree) delete _tree;
    }

    ///Sets the cost map.

    ///Sets the cost map.
    ///\return <tt> (*this) </tt>
    FredmanTarjan &costMap(const CostMap &m){
      cost = &m;
      return *this;
    }

    ///Sets the map storing the tree edges.

    ///Sets the map storing the tree edges.
    ///If you don't use this function before calling \ref run(),
    ///it will allocate one. The destuctor deallocates this
    ///automatically allocated map, of course.
    ///By default this is a BoolEdgeMap.
    ///\return <tt> (*this) </tt>
    FredmanTarjan &treeMap(TreeMap &m){
      if(local_tree) {
	delete _tree;
	local_tree=false;
      }
      _tree = &m;
      return *this;
    }

  public:
    ///\name Execution control
    ///The simplest way to execute the algorithm is to use
    ///one of the member functions called \c run(...).

    ///@{

    ///Initializes the internal data structures.

    ///Initializes the internal data structures.
    ///
    void init(){
      create_maps();
      for(typename Graph::UEdgeIt i(*graph);i!=INVALID;++i){
	_tree->set(i,false);
      }
    }

    ///Executes the algorithm.

    ///Executes the algorithm.
    ///
    ///\pre init() must be called and at least one node should be added
    ///with addSource() before using this function.
    ///
    ///This method runs the %FredmanTarjan algorithm from the node(s)
    ///in order to compute the
    ///minimum spanning tree.
    void start(){
	phase(*graph,identityMap<UEdge>(),countEdges(*graph));
    }
    
    ///Runs %FredmanTarjan algorithm.
    
    ///This method runs the %FredmanTarjan algorithm
    ///in order to compute the minimum spanning forest.
    ///
    ///\note ft.run() is just a shortcut of the following code.
    ///\code
    ///  ft.init();
    ///  ft.start();
    ///\endcode
    void run() {
      init();
      start();
    }

    ///@}

    ///\name Query Functions
    ///The result of the %FredmanTarjan algorithm can be obtained using these
    ///functions.\n
    ///Before the use of these functions,
    ///either run() or start() must be called.
    
    ///@{

    ///Returns a reference to the tree edges map.

    ///Returns a reference to the TreeEdgeMap of the edges of the
    ///minimum spanning tree. The value of the map is \c true only if the 
    ///edge is in the minimum spanning tree.
    ///
    ///\pre \ref run() or \ref start() must be called before using this 
    ///function.
    const TreeMap &treeMap() const { return *_tree;}
 
    ///Sets the tree edges map.

    ///Sets the TreeMap of the edges of the minimum spanning tree.
    ///The map values belonging to the edges of the minimum
    ///spanning tree are set to \c tree_edge_value or \c true by default 
    ///while the edge values not belonging to the minimum spanning tree are 
    ///set to
    ///\c tree_default_value or \c false by default.
    ///
    ///\pre \ref run() or \ref start() must be called before using this 
    ///function.

    template<class TreeMap>
    void treeEdges(
        TreeMap& tree,
        const typename TreeMap::Value& tree_edge_value=true,
        const typename TreeMap::Value& tree_default_value=false) const {
      for(typename UGraph::UEdgeIt i(*graph);i!=INVALID;++i){
	(*_tree)[i]?tree.set(i,tree_edge_value):tree.set(i,tree_default_value);
      }
    }

    ///\brief Checks if an edge is in the spanning tree or not.

    ///Checks if an edge is in the spanning tree or not.
    ///\param e is the edge that will be checked
    ///\return \c true if e is in the spanning tree, \c false otherwise
    bool tree(UEdge e){
      return (*_tree)[e];
    }
    ///@}
  };

  /// \ingroup spantree
  ///
  /// \brief Function type interface for FredmanTarjan algorithm.
  ///
  /// Function type interface for FredmanTarjan algorithm.
  /// \param graph the UGraph that the algorithm runs on
  /// \param cost the CostMap of the edges
  /// \retval tree the EdgeMap that contains whether an edge is in the 
  /// spanning tree or not
  ///
  /// \sa Prim
  template<class Graph,class CostMap,class TreeMap>
  void fredmanTarjan(const Graph& graph, const CostMap& cost,TreeMap& tree){
    typename FredmanTarjan<Graph,CostMap>::template DefTreeMap<TreeMap>::
      Create ft(graph,cost);
    ft.treeMap(tree);
    ft.run();
  }

} //END OF NAMESPACE LEMON

#endif
