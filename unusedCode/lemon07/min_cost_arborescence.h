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

#ifndef LEMON_MIN_COST_ARBORESCENCE_H
#define LEMON_MIN_COST_ARBORESCENCE_H

///\ingroup spantree
///\file
///\brief Minimum Cost Arborescence algorithm.

#include <vector>

#include <lemon/list_graph.h>
#include <lemon/bin_heap.h>

namespace lemon {


  /// \brief Default traits class of MinCostArborescence class.
  /// 
  /// Default traits class of MinCostArborescence class.
  /// \param _Graph Graph type.
  /// \param _CostMap Type of cost map.
  template <class _Graph, class _CostMap>
  struct MinCostArborescenceDefaultTraits{

    /// \brief The graph type the algorithm runs on. 
    typedef _Graph Graph;

    /// \brief The type of the map that stores the edge costs.
    ///
    /// The type of the map that stores the edge costs.
    /// It must meet the \ref concepts::ReadMap "ReadMap" concept.
    typedef _CostMap CostMap;

    /// \brief The value type of the costs.
    ///
    /// The value type of the costs.
    typedef typename CostMap::Value Value;

    /// \brief The type of the map that stores which edges are 
    /// in the arborescence.
    ///
    /// The type of the map that stores which edges are in the arborescence.
    /// It must meet the \ref concepts::WriteMap "WriteMap" concept.
    /// Initially it will be set to false on each edge. After it
    /// will set all arborescence edges once.
    typedef typename Graph::template EdgeMap<bool> ArborescenceMap; 

    /// \brief Instantiates a ArborescenceMap.
    ///
    /// This function instantiates a \ref ArborescenceMap. 
    /// \param _graph is the graph, to which we would like to define the 
    /// ArborescenceMap.
    static ArborescenceMap *createArborescenceMap(const Graph &_graph){
      return new ArborescenceMap(_graph);
    }

    /// \brief The type of the PredMap
    ///
    /// The type of the PredMap. It is a node map with an edge value type.
    typedef typename Graph::template NodeMap<typename Graph::Edge> PredMap;

    /// \brief Instantiates a PredMap.
    ///
    /// This function instantiates a \ref PredMap. 
    /// \param _graph is the graph, to which we would like to define the 
    /// PredMap.
    static PredMap *createPredMap(const Graph &_graph){
      return new PredMap(_graph);
    }
    
  };

  /// \ingroup spantree
  ///
  /// \brief %MinCostArborescence algorithm class.
  ///
  /// This class provides an efficient implementation of 
  /// %MinCostArborescence algorithm. The arborescence is a tree 
  /// which is directed from a given source node of the graph. One or
  /// more sources should be given for the algorithm and it will calculate
  /// the minimum cost subgraph which are union of arborescences with the
  /// given sources and spans all the nodes which are reachable from the
  /// sources. The time complexity of the algorithm is \f$ O(n^2+e) \f$.
  ///
  /// The algorithm provides also an optimal dual solution to arborescence
  /// that way the optimality of the solution can be proofed easily.
  ///
  /// \param _Graph The graph type the algorithm runs on. The default value
  /// is \ref ListGraph. The value of _Graph is not used directly by
  /// MinCostArborescence, it is only passed to 
  /// \ref MinCostArborescenceDefaultTraits.
  /// \param _CostMap This read-only EdgeMap determines the costs of the
  /// edges. It is read once for each edge, so the map may involve in
  /// relatively time consuming process to compute the edge cost if
  /// it is necessary. The default map type is \ref
  /// concepts::Graph::EdgeMap "Graph::EdgeMap<int>".  The value
  /// of _CostMap is not used directly by MinCostArborescence, 
  /// it is only passed to \ref MinCostArborescenceDefaultTraits.  
  /// \param _Traits Traits class to set various data types used 
  /// by the algorithm.  The default traits class is 
  /// \ref MinCostArborescenceDefaultTraits
  /// "MinCostArborescenceDefaultTraits<_Graph,_CostMap>".  See \ref
  /// MinCostArborescenceDefaultTraits for the documentation of a 
  /// MinCostArborescence traits class.
  ///
  /// \author Balazs Dezso
#ifndef DOXYGEN
  template <typename _Graph = ListGraph, 
            typename _CostMap = typename _Graph::template EdgeMap<int>,
            typename _Traits = 
            MinCostArborescenceDefaultTraits<_Graph, _CostMap> >
#else 
  template <typename _Graph, typename _CostMap, typedef _Traits>
#endif
  class MinCostArborescence {
  public:

    /// \brief \ref Exception for uninitialized parameters.
    ///
    /// This error represents problems in the initialization
    /// of the parameters of the algorithms.    
    class UninitializedParameter : public lemon::UninitializedParameter {
    public:
      virtual const char* what() const throw() {
	return "lemon::MinCostArborescence::UninitializedParameter";
      }
    };

    /// The traits.
    typedef _Traits Traits;
    /// The type of the underlying graph.
    typedef typename Traits::Graph Graph;
    /// The type of the map that stores the edge costs.
    typedef typename Traits::CostMap CostMap;
    ///The type of the costs of the edges.
    typedef typename Traits::Value Value;
    ///The type of the predecessor map.
    typedef typename Traits::PredMap PredMap;
    ///The type of the map that stores which edges are in the arborescence.
    typedef typename Traits::ArborescenceMap ArborescenceMap;

  protected:

    typedef typename Graph::Node Node;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::NodeIt NodeIt;
    typedef typename Graph::EdgeIt EdgeIt;
    typedef typename Graph::InEdgeIt InEdgeIt;
    typedef typename Graph::OutEdgeIt OutEdgeIt;

    struct CostEdge {

      Edge edge;
      Value value;

      CostEdge() {}
      CostEdge(Edge _edge, Value _value) : edge(_edge), value(_value) {}

    };

    const Graph *graph;
    const CostMap *cost;

    PredMap *_pred;
    bool local_pred;

    ArborescenceMap *_arborescence;
    bool local_arborescence;

    typedef typename Graph::template EdgeMap<int> EdgeOrder;
    EdgeOrder *_edge_order;
    
    typedef typename Graph::template NodeMap<int> NodeOrder;
    NodeOrder *_node_order;

    typedef typename Graph::template NodeMap<CostEdge> CostEdgeMap;
    CostEdgeMap *_cost_edges; 

    struct StackLevel {

      std::vector<CostEdge> edges;
      int node_level;

    };

    std::vector<StackLevel> level_stack;    
    std::vector<Node> queue;

    typedef std::vector<typename Graph::Node> DualNodeList;

    DualNodeList _dual_node_list;

    struct DualVariable {
      int begin, end;
      Value value;
      
      DualVariable(int _begin, int _end, Value _value)
        : begin(_begin), end(_end), value(_value) {}

    };

    typedef std::vector<DualVariable> DualVariables;

    DualVariables _dual_variables;

    typedef typename Graph::template NodeMap<int> HeapCrossRef;
    
    HeapCrossRef *_heap_cross_ref;

    typedef BinHeap<int, HeapCrossRef> Heap;

    Heap *_heap;

  public:

    /// \name Named template parameters

    /// @{

    template <class T>
    struct DefArborescenceMapTraits : public Traits {
      typedef T ArborescenceMap;
      static ArborescenceMap *createArborescenceMap(const Graph &)
      {
	throw UninitializedParameter();
      }
    };

    /// \brief \ref named-templ-param "Named parameter" for 
    /// setting ArborescenceMap type
    ///
    /// \ref named-templ-param "Named parameter" for setting 
    /// ArborescenceMap type
    template <class T>
    struct DefArborescenceMap 
      : public MinCostArborescence<Graph, CostMap,
                                   DefArborescenceMapTraits<T> > {
      typedef MinCostArborescence<Graph, CostMap, 
                                   DefArborescenceMapTraits<T> > Create;
    };

    template <class T>
    struct DefPredMapTraits : public Traits {
      typedef T PredMap;
      static PredMap *createPredMap(const Graph &)
      {
	throw UninitializedParameter();
      }
    };

    /// \brief \ref named-templ-param "Named parameter" for 
    /// setting PredMap type
    ///
    /// \ref named-templ-param "Named parameter" for setting 
    /// PredMap type
    template <class T>
    struct DefPredMap 
      : public MinCostArborescence<Graph, CostMap, DefPredMapTraits<T> > {
      typedef MinCostArborescence<Graph, CostMap, DefPredMapTraits<T> > Create;
    };
    
    /// @}

    /// \brief Constructor.
    ///
    /// \param _graph The graph the algorithm will run on.
    /// \param _cost The cost map used by the algorithm.
    MinCostArborescence(const Graph& _graph, const CostMap& _cost) 
      : graph(&_graph), cost(&_cost), _pred(0), local_pred(false),
        _arborescence(0), local_arborescence(false), 
        _edge_order(0), _node_order(0), _cost_edges(0), 
        _heap_cross_ref(0), _heap(0) {}

    /// \brief Destructor.
    ~MinCostArborescence() {
      destroyStructures();
    }

    /// \brief Sets the arborescence map.
    /// 
    /// Sets the arborescence map.
    /// \return \c (*this)
    MinCostArborescence& arborescenceMap(ArborescenceMap& m) {
      if (local_arborescence) {
        delete _arborescence;
      }
      local_arborescence = false;
      _arborescence = &m;
      return *this;
    }

    /// \brief Sets the arborescence map.
    /// 
    /// Sets the arborescence map.
    /// \return \c (*this)
    MinCostArborescence& predMap(PredMap& m) {
      if (local_pred) {
        delete _pred;
      }
      local_pred = false;
      _pred = &m;
      return *this;
    }

    /// \name Query Functions
    /// The result of the %MinCostArborescence algorithm can be obtained 
    /// using these functions.\n
    /// Before the use of these functions,
    /// either run() or start() must be called.
    
    /// @{

    /// \brief Returns a reference to the arborescence map.
    ///
    /// Returns a reference to the arborescence map.
    const ArborescenceMap& arborescenceMap() const {
      return *_arborescence;
    }

    /// \brief Returns true if the edge is in the arborescence.
    ///
    /// Returns true if the edge is in the arborescence.
    /// \param edge The edge of the graph.
    /// \pre \ref run() must be called before using this function.
    bool arborescence(Edge edge) const {
      return (*_pred)[graph->target(edge)] == edge;
    }

    /// \brief Returns a reference to the pred map.
    ///
    /// Returns a reference to the pred map.
    const PredMap& predMap() const {
      return *_pred;
    }

    /// \brief Returns the predecessor edge of the given node.
    ///
    /// Returns the predecessor edge of the given node.
    bool pred(Node node) const {
      return (*_pred)[node];
    }
 
    /// \brief Returns the cost of the arborescence.
    ///
    /// Returns the cost of the arborescence.
    Value arborescenceValue() const {
      Value sum = 0;
      for (EdgeIt it(*graph); it != INVALID; ++it) {
        if (arborescence(it)) {
          sum += (*cost)[it];
        }
      }
      return sum;
    }

    /// \brief Indicates that a node is reachable from the sources.
    ///
    /// Indicates that a node is reachable from the sources.
    bool reached(Node node) const {
      return (*_node_order)[node] != -3;
    }

    /// \brief Indicates that a node is processed.
    ///
    /// Indicates that a node is processed. The arborescence path exists 
    /// from the source to the given node.
    bool processed(Node node) const {
      return (*_node_order)[node] == -1;
    }

    /// \brief Returns the number of the dual variables in basis.
    ///
    /// Returns the number of the dual variables in basis.
    int dualSize() const {
      return _dual_variables.size();
    }

    /// \brief Returns the value of the dual solution.
    ///
    /// Returns the value of the dual solution. It should be
    /// equal to the arborescence value.
    Value dualValue() const {
      Value sum = 0;
      for (int i = 0; i < int(_dual_variables.size()); ++i) {
        sum += _dual_variables[i].value;
      }
      return sum;
    }

    /// \brief Returns the number of the nodes in the dual variable.
    ///
    /// Returns the number of the nodes in the dual variable.
    int dualSize(int k) const {
      return _dual_variables[k].end - _dual_variables[k].begin;
    }

    /// \brief Returns the value of the dual variable.
    ///
    /// Returns the the value of the dual variable.
    const Value& dualValue(int k) const {
      return _dual_variables[k].value;
    }

    /// \brief Lemon iterator for get a dual variable.
    ///
    /// Lemon iterator for get a dual variable. This class provides
    /// a common style lemon iterator which gives back a subset of
    /// the nodes.
    class DualIt {
    public:

      /// \brief Constructor.
      ///
      /// Constructor for get the nodeset of the variable. 
      DualIt(const MinCostArborescence& algorithm, int variable) 
        : _algorithm(&algorithm), _variable(variable) 
      {
        _index = _algorithm->_dual_variables[_variable].begin;
      }

      /// \brief Invalid constructor.
      ///
      /// Invalid constructor.
      DualIt(Invalid) : _algorithm(0) {}

      /// \brief Conversion to node.
      ///
      /// Conversion to node.
      operator Node() const { 
        return _algorithm ? _algorithm->_dual_node_list[_index] : INVALID;
      }

      /// \brief Increment operator.
      ///
      /// Increment operator.
      DualIt& operator++() {
        ++_index;
        if (_algorithm->_dual_variables[_variable].end == _index) {
          _algorithm = 0;
        }
        return *this; 
      }

      bool operator==(const DualIt& it) const { 
        return static_cast<Node>(*this) == static_cast<Node>(it); 
      }
      bool operator!=(const DualIt& it) const { 
        return static_cast<Node>(*this) != static_cast<Node>(it); 
      }
      bool operator<(const DualIt& it) const { 
        return static_cast<Node>(*this) < static_cast<Node>(it); 
      }
      
    private:
      const MinCostArborescence* _algorithm;
      int _variable;
      int _index;
    };

    /// @}
    
    /// \name Execution control
    /// The simplest way to execute the algorithm is to use
    /// one of the member functions called \c run(...). \n
    /// If you need more control on the execution,
    /// first you must call \ref init(), then you can add several 
    /// source nodes with \ref addSource().
    /// Finally \ref start() will perform the arborescence
    /// computation.

    ///@{

    /// \brief Initializes the internal data structures.
    ///
    /// Initializes the internal data structures.
    ///
    void init() {
      initStructures();
      _heap->clear();
      for (NodeIt it(*graph); it != INVALID; ++it) {
        (*_cost_edges)[it].edge = INVALID;
        _node_order->set(it, -3); 
        _heap_cross_ref->set(it, Heap::PRE_HEAP);
        _pred->set(it, INVALID);
      }
      for (EdgeIt it(*graph); it != INVALID; ++it) {
        _arborescence->set(it, false);
        _edge_order->set(it, -1);
      }
      _dual_node_list.clear();
      _dual_variables.clear();
    }

    /// \brief Adds a new source node.
    ///
    /// Adds a new source node to the algorithm.
    void addSource(Node source) {
      std::vector<Node> nodes;
      nodes.push_back(source);
      while (!nodes.empty()) {
        Node node = nodes.back();
        nodes.pop_back();
        for (OutEdgeIt it(*graph, node); it != INVALID; ++it) {
          Node target = graph->target(it);
          if ((*_node_order)[target] == -3) {
            (*_node_order)[target] = -2;
            nodes.push_back(target);
            queue.push_back(target);
          }
        }
      }
      (*_node_order)[source] = -1;
    }

    /// \brief Processes the next node in the priority queue.
    ///
    /// Processes the next node in the priority queue.
    ///
    /// \return The processed node.
    ///
    /// \warning The queue must not be empty!
    Node processNextNode() {
      Node node = queue.back();
      queue.pop_back();
      if ((*_node_order)[node] == -2) {
        Edge edge = prepare(node);
        Node source = graph->source(edge);
        while ((*_node_order)[source] != -1) {
          if ((*_node_order)[source] >= 0) {
            edge = contract(source);
          } else {
            edge = prepare(source);
          }
          source = graph->source(edge);
        }
        finalize(edge);
        level_stack.clear();        
      }
      return node;
    }

    /// \brief Returns the number of the nodes to be processed.
    ///
    /// Returns the number of the nodes to be processed.
    int queueSize() const {
      return queue.size();
    }

    /// \brief Returns \c false if there are nodes to be processed.
    ///
    /// Returns \c false if there are nodes to be processed.
    bool emptyQueue() const {
      return queue.empty();
    }

    /// \brief Executes the algorithm.
    ///
    /// Executes the algorithm.
    ///
    /// \pre init() must be called and at least one node should be added
    /// with addSource() before using this function.
    ///
    ///\note mca.start() is just a shortcut of the following code.
    ///\code
    ///while (!mca.emptyQueue()) {
    ///  mca.processNextNode();
    ///}
    ///\endcode
    void start() {
      while (!emptyQueue()) {
        processNextNode();
      }
    }

    /// \brief Runs %MinCostArborescence algorithm from node \c s.
    /// 
    /// This method runs the %MinCostArborescence algorithm from 
    /// a root node \c s.
    ///
    ///\note mca.run(s) is just a shortcut of the following code.
    ///\code
    ///mca.init();
    ///mca.addSource(s);
    ///mca.start();
    ///\endcode
    void run(Node node) {
      init();
      addSource(node);
      start();
    }

    ///@}

  protected:

    void initStructures() {
      if (!_pred) {
        local_pred = true;
        _pred = Traits::createPredMap(*graph);
      }
      if (!_arborescence) {
        local_arborescence = true;
        _arborescence = Traits::createArborescenceMap(*graph);
      }
      if (!_edge_order) {
        _edge_order = new EdgeOrder(*graph);
      }
      if (!_node_order) {
        _node_order = new NodeOrder(*graph);
      }
      if (!_cost_edges) {
        _cost_edges = new CostEdgeMap(*graph);
      }
      if (!_heap_cross_ref) {
        _heap_cross_ref = new HeapCrossRef(*graph, -1);
      }
      if (!_heap) {
        _heap = new Heap(*_heap_cross_ref);
      }
    }

    void destroyStructures() {
      if (local_arborescence) {
        delete _arborescence;
      }
      if (local_pred) {
        delete _pred;
      }
      if (!_edge_order) {
        delete _edge_order;
      }
      if (_node_order) {
        delete _node_order;
      }
      if (!_cost_edges) {
        delete _cost_edges;
      }
      if (!_heap) {
        delete _heap;
      }
      if (!_heap_cross_ref) {
        delete _heap_cross_ref;
      }
    }

    Edge prepare(Node node) {
      std::vector<Node> nodes;
      (*_node_order)[node] = _dual_node_list.size();
      StackLevel level;
      level.node_level = _dual_node_list.size();
      _dual_node_list.push_back(node);
      for (InEdgeIt it(*graph, node); it != INVALID; ++it) {
        Edge edge = it;
        Node source = graph->source(edge);
        Value value = (*cost)[it];
        if (source == node || (*_node_order)[source] == -3) continue;
        if ((*_cost_edges)[source].edge == INVALID) {
          (*_cost_edges)[source].edge = edge;
          (*_cost_edges)[source].value = value;
          nodes.push_back(source);
        } else {
          if ((*_cost_edges)[source].value > value) {
            (*_cost_edges)[source].edge = edge;
            (*_cost_edges)[source].value = value;
          }
        }      
      }
      CostEdge minimum = (*_cost_edges)[nodes[0]]; 
      for (int i = 1; i < int(nodes.size()); ++i) {
        if ((*_cost_edges)[nodes[i]].value < minimum.value) {
          minimum = (*_cost_edges)[nodes[i]];
        }
      }
      _edge_order->set(minimum.edge, _dual_variables.size());
      DualVariable var(_dual_node_list.size() - 1, 
                       _dual_node_list.size(), minimum.value);
      _dual_variables.push_back(var);
      for (int i = 0; i < int(nodes.size()); ++i) {
        (*_cost_edges)[nodes[i]].value -= minimum.value;
        level.edges.push_back((*_cost_edges)[nodes[i]]);
        (*_cost_edges)[nodes[i]].edge = INVALID;
      }
      level_stack.push_back(level);
      return minimum.edge;
    }
  
    Edge contract(Node node) {
      int node_bottom = bottom(node);
      std::vector<Node> nodes;
      while (!level_stack.empty() && 
             level_stack.back().node_level >= node_bottom) {
        for (int i = 0; i < int(level_stack.back().edges.size()); ++i) {
          Edge edge = level_stack.back().edges[i].edge;
          Node source = graph->source(edge);
          Value value = level_stack.back().edges[i].value;
          if ((*_node_order)[source] >= node_bottom) continue;
          if ((*_cost_edges)[source].edge == INVALID) {
            (*_cost_edges)[source].edge = edge;
            (*_cost_edges)[source].value = value;
            nodes.push_back(source);
          } else {
            if ((*_cost_edges)[source].value > value) {
              (*_cost_edges)[source].edge = edge;
              (*_cost_edges)[source].value = value;
            }
          }
        }
        level_stack.pop_back();
      }
      CostEdge minimum = (*_cost_edges)[nodes[0]]; 
      for (int i = 1; i < int(nodes.size()); ++i) {
        if ((*_cost_edges)[nodes[i]].value < minimum.value) {
          minimum = (*_cost_edges)[nodes[i]];
        }
      }
      _edge_order->set(minimum.edge, _dual_variables.size());
      DualVariable var(node_bottom, _dual_node_list.size(), minimum.value);
      _dual_variables.push_back(var);
      StackLevel level;
      level.node_level = node_bottom;
      for (int i = 0; i < int(nodes.size()); ++i) {
        (*_cost_edges)[nodes[i]].value -= minimum.value;
        level.edges.push_back((*_cost_edges)[nodes[i]]);
        (*_cost_edges)[nodes[i]].edge = INVALID;
      }
      level_stack.push_back(level);
      return minimum.edge;
    }

    int bottom(Node node) {
      int k = level_stack.size() - 1;
      while (level_stack[k].node_level > (*_node_order)[node]) {
        --k;
      }
      return level_stack[k].node_level;
    }

    void finalize(Edge edge) {
      Node node = graph->target(edge);
      _heap->push(node, (*_edge_order)[edge]);
      _pred->set(node, edge);
      while (!_heap->empty()) {
        Node source = _heap->top();
        _heap->pop();
        _node_order->set(source, -1);
        for (OutEdgeIt it(*graph, source); it != INVALID; ++it) {
          if ((*_edge_order)[it] < 0) continue;
          Node target = graph->target(it);
          switch(_heap->state(target)) {
          case Heap::PRE_HEAP:
            _heap->push(target, (*_edge_order)[it]); 
            _pred->set(target, it);
            break;
          case Heap::IN_HEAP:
            if ((*_edge_order)[it] < (*_heap)[target]) {
              _heap->decrease(target, (*_edge_order)[it]); 
              _pred->set(target, it);
            }
            break;
          case Heap::POST_HEAP:
            break;
          }
        }
        _arborescence->set((*_pred)[source], true);
      }
    }

  };

  /// \ingroup spantree
  ///
  /// \brief Function type interface for MinCostArborescence algorithm.
  ///
  /// Function type interface for MinCostArborescence algorithm.
  /// \param graph The Graph that the algorithm runs on.
  /// \param cost The CostMap of the edges.
  /// \param source The source of the arborescence.
  /// \retval arborescence The bool EdgeMap which stores the arborescence.
  /// \return The cost of the arborescence. 
  ///
  /// \sa MinCostArborescence
  template <typename Graph, typename CostMap, typename ArborescenceMap>
  typename CostMap::Value minCostArborescence(const Graph& graph, 
                                              const CostMap& cost,
                                              typename Graph::Node source,
                                              ArborescenceMap& arborescence) {
    typename MinCostArborescence<Graph, CostMap>
      ::template DefArborescenceMap<ArborescenceMap>
      ::Create mca(graph, cost);
    mca.arborescenceMap(arborescence);
    mca.run(source);
    return mca.arborescenceValue();
  }

}

#endif
