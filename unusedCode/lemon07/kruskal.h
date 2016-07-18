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

#ifndef LEMON_KRUSKAL_H
#define LEMON_KRUSKAL_H

#include <algorithm>
#include <vector>
#include <lemon/unionfind.h>
#include <lemon/graph_utils.h>
#include <lemon/maps.h>

#include <lemon/radix_sort.h>

#include <lemon/bits/utility.h>
#include <lemon/bits/traits.h>

///\ingroup spantree
///\file
///\brief Kruskal's algorithm to compute a minimum cost tree
///
///Kruskal's algorithm to compute a minimum cost tree.
///

namespace lemon {

  namespace _kruskal_bits {

    template <typename Map, typename Comp>
    struct MappedComp {

      typedef typename Map::Key Key;

      const Map& map;
      Comp comp;

      MappedComp(const Map& _map) 
        : map(_map), comp() {}

      bool operator()(const Key& left, const Key& right) {
        return comp(map[left], map[right]);
      }

    };
  
  }

  /// \brief Default traits class of Kruskal class.
  ///
  /// Default traits class of Kruskal class.
  /// \param _UGraph Undirected graph type.
  /// \param _CostMap Type of cost map.
  template <typename _UGraph, typename _CostMap>
  struct KruskalDefaultTraits{

    /// \brief The graph type the algorithm runs on. 
    typedef _UGraph UGraph;

    /// \brief The type of the map that stores the edge costs.
    ///
    /// The type of the map that stores the edge costs.
    /// It must meet the \ref concepts::ReadMap "ReadMap" concept.
    typedef _CostMap CostMap;

    /// \brief The type of the cost of the edges.
    typedef typename _CostMap::Value Value;

    /// \brief The type of the map that stores whether an edge is in the
    /// spanning tree or not.
    ///
    /// The type of the map that stores whether an edge is in the
    /// spanning tree or not.
    typedef typename _UGraph::template UEdgeMap<bool> TreeMap;

    /// \brief Instantiates a TreeMap.
    ///
    /// This function instantiates a \ref TreeMap.
    ///
    /// The first parameter is the graph, to which
    /// we would like to define the \ref TreeMap
    static TreeMap *createTreeMap(const _UGraph& graph){
      return new TreeMap(graph);
    }

    template <typename Iterator>
    static void sort(Iterator begin, Iterator end, const CostMap& cost) {
      _kruskal_bits::MappedComp<CostMap, std::less<Value> > comp(cost);
      std::sort(begin, end, comp);
    }

  };

  ///\ingroup spantree
  ///
  /// \brief Kruskal's algorithm to find a minimum cost tree of a graph.
  ///
  /// This class implements Kruskal's algorithm to find a minimum cost
  /// spanning tree. The 
  ///
  /// \param _UGraph Undirected graph type.
  /// \param _CostMap Type of cost map.
  template <typename _UGraph, typename _CostMap,
            typename _Traits = KruskalDefaultTraits<_UGraph, _CostMap> >
  class Kruskal {
  public:

    typedef _Traits Traits;

    typedef typename _Traits::UGraph UGraph;
    typedef typename _Traits::CostMap CostMap;

    typedef typename _Traits::TreeMap TreeMap;

    typedef typename _Traits::Value Value;

    template <typename Comp>
    struct DefSortCompareTraits : public Traits {
      template <typename Iterator>
      static void sort(Iterator begin, Iterator end, const CostMap& cost) {
        _kruskal_bits::MappedComp<CostMap, Comp> comp(cost, Comp());
        std::sort(begin, end, comp);
      }
    };

    /// \brief \ref named-templ-param "Named parameter" for setting the
    /// comparator object of the standard sort
    ///
    /// \ref named-templ-param "Named parameter" for setting the
    /// comparator object of the standard sort
    template <typename Comp>
    struct DefSortCompare
      : public Kruskal<UGraph, CostMap, DefSortCompareTraits<Comp> > {
      typedef Kruskal<UGraph, CostMap, DefSortCompareTraits<Comp> > Create;
    };    

    struct DefRadixSortTraits : public Traits {
      template <typename Iterator>
      static void sort(Iterator begin, Iterator end, const CostMap& cost) {
        radixSort(begin, end, cost);
      }
    };

    /// \brief \ref named-templ-param "Named parameter" for setting the
    /// sort function to radix sort
    ///
    /// \brief \ref named-templ-param "Named parameter" for setting the
    /// sort function to radix sort. The value type of the cost map should
    /// be integral, of course. 
    struct DefRadixSort
      : public Kruskal<UGraph, CostMap, DefRadixSortTraits> {
      typedef Kruskal<UGraph, CostMap, DefRadixSortTraits> Create;
    };    

    template <class TM>
    struct DefTreeMapTraits : public Traits {
      typedef TM TreeMap;
      static TreeMap *createTreeMap(const UGraph &) {
        throw UninitializedParameter();
      }
    };

    /// \brief \ref named-templ-param "Named parameter" for setting
    /// TreeMap
    ///
    /// \ref named-templ-param "Named parameter" for setting TreeMap
    ///
    template <class TM>
    struct DefTreeMap
      : public Kruskal< UGraph, CostMap, DefTreeMapTraits<TM> > {
      typedef Kruskal< UGraph, CostMap, DefTreeMapTraits<TM> > Create;
    };    
    

  private:

    typedef typename UGraph::Node Node;
    typedef typename UGraph::NodeIt NodeIt;

    typedef typename UGraph::UEdge UEdge;
    typedef typename UGraph::UEdgeIt UEdgeIt;

    const UGraph& graph;
    const CostMap& cost;
    
    std::vector<UEdge> edges;
    
    typedef typename UGraph::template NodeMap<int> UfIndex;
    typedef UnionFind<UfIndex> Uf;
    UfIndex *ufi;
    Uf *uf;

    int index;

    void initStructures() {
      if (!_tree) {
        _tree = Traits::createTreeMap(graph);
        local_tree = true;
      }
      if (!uf) {
        ufi = new typename UGraph::template NodeMap<int>(graph);
        uf = new UnionFind<typename UGraph::template NodeMap<int> >(*ufi);
      }
    }
    
    void initUnionFind() {
      uf->clear();
      for (NodeIt it(graph); it != INVALID; ++it) {
        uf->insert(it);
      }
    }

    bool local_tree;
    TreeMap* _tree;

  public:

    /// \brief Constructor
    ///
    /// Constructor of the algorithm.
    Kruskal(const UGraph& _graph, const CostMap& _cost) 
      : graph(_graph), cost(_cost),
        ufi(0), uf(0), local_tree(false), _tree(0) {}

    /// \brief Destructor
    ///
    /// Destructor
    ~Kruskal() {
      if (local_tree) {
        delete _tree;
      }
      if (uf) {
        delete uf;
        delete ufi;
      }
    }

    /// \brief Sets the map storing the tree edges.
    ///
    /// Sets the map storing the tree edges.
    /// If you don't use this function before calling \ref run(),
    /// it will allocate one. The destuctor deallocates this
    /// automatically allocated map, of course.
    /// \return \c *this </tt>
    Kruskal& treeMap(TreeMap &m){
      if (local_tree) {
	delete _tree;
	local_tree = false;
      }
      _tree = &m;
      return *this;
    }

    /// \brief Initialize the algorithm
    ///
    /// This member function initializes the unionfind data structure
    /// and sorts the edges into ascending order
    void init() {
      initStructures();
      initUnionFind();
      for (UEdgeIt e(graph); e != INVALID; ++e) {
        edges.push_back(e);
        _tree->set(e, false);
      }      
      Traits::sort(edges.begin(), edges.end(), cost); 
      index = 0;
    }


    /// \brief Initialize the algorithm
    ///
    /// This member function initializes the unionfind data structure
    /// and sets the edge order to the given sequence. The given
    /// sequence should be a valid STL range of undirected edges.
    template <typename Iterator>
    void initPresorted(Iterator begin, Iterator end) {
      initStructures();
      initUnionFind();
      edges.clear();
      std::copy(begin, end, std::back_inserter(edges));
      index = 0;
    }

    /// \brief Initialize the algorithm
    ///
    /// This member function initializes the unionfind data structure
    /// and sets the tree to empty. It does not change the order of
    /// the edges, it uses the order of the previous running.
    void reinit() {
      initStructures();
      initUnionFind();
      for (UEdgeIt e(graph); e != INVALID; ++e) {
        _tree->set(e, false);
      }
      index = 0;
    }


    /// \brief Executes the algorithm.
    ///
    /// Executes the algorithm.
    ///
    /// \pre init() must be called before using this function.
    ///
    /// This method runs the %Kruskal algorithm.
    void start() {
      while (index < int(edges.size())) {
        if (uf->join(graph.target(edges[index]), graph.source(edges[index]))) {
          _tree->set(edges[index], true);
        }
        ++index;
      }
    }

    /// \brief Runs the prim algorithm until it find a new tree edge
    ///
    /// Runs the prim algorithm until it find a new tree edge. If it
    /// does not next tree edge in the sequence it gives back \c INVALID.
    UEdge findNextTreeEdge() {
      while (index < int(edges.size())) {
        if (uf->join(graph.target(edges[index]), graph.source(edges[index]))) {
          _tree->set(edges[index], true);
          return edges[index++];
        }        
        ++index;
      }
      return INVALID;
    }
      
    /// \brief Processes the next edge in the sequence
    ///
    /// Processes the next edge in the sequence.
    ///
    /// \return The prcocessed edge.
    ///
    /// \warning The sequence must not be empty!
    UEdge processNextEdge() {
      UEdge edge = edges[index++];
      processEdge(edge);
      return edge;
    }

    /// \brief Processes an arbitrary edge
    ///
    /// Processes the next edge in the sequence.
    ///
    /// \return True when the edge is a tree edge.
    bool processEdge(const UEdge& edge) {
      if (uf->join(graph.target(edge), graph.source(edge))) {
        _tree->set(edge, true);
        return true;
      } else {
        return false;
      }    
    }

    /// \brief Returns \c false if there are edge to be processed in
    /// sequence
    ///
    /// Returns \c false if there are nodes to be processed in the
    /// sequence
    bool emptyQueue() {
      return index == int(edges.size());
    }

    /// \brief Returns the next edge to be processed
    ///
    /// Returns the next edge to be processed
    ///
    UEdge nextEdge() const {
      return edges[index];
    }

    /// \brief Runs %Kruskal algorithm.
    ///
    /// This method runs the %Kruskal algorithm in order to compute the
    /// minimum spanning tree (or minimum spanning forest).  The
    /// method also works on graphs that has more than one components.
    /// In this case it computes the minimum spanning forest.
    void run() {
      init();
      start();
    }

    /// \brief Returns a reference to the tree edges map
    ///
    /// Returns a reference to the TreeEdgeMap of the edges of the
    /// minimum spanning tree. The value of the map is \c true only if
    /// the edge is in the minimum spanning tree.
    ///
    const TreeMap &treeMap() const { return *_tree;}

    /// \brief Returns the total cost of the tree
    ///
    /// Returns the total cost of the tree
    Value treeValue() const {
      Value value = 0;
      for (UEdgeIt it(graph); it != INVALID; ++it) {
        if ((*_tree)[it]) {
          value += cost[it];
        }
      }
      return value;
    }

    /// \brief Returns true when the given edge is tree edge
    ///
    /// Returns true when the given edge is tree edge
    bool tree(UEdge e) const {
      return (*_tree)[e];
    }
    
    
  };


  namespace _kruskal_bits {

    template <typename Graph, typename In, typename Out>
    typename In::value_type::second_type
    kruskal(const Graph& graph, const In& in, Out& out) {
      typedef typename In::value_type::second_type Value;
      typedef typename Graph::template NodeMap<int> IndexMap;
      typedef typename Graph::Node Node;
      
      IndexMap index(graph);
      UnionFind<IndexMap> uf(index);
      for (typename Graph::NodeIt it(graph); it != INVALID; ++it) {
        uf.insert(it);
      }
      
      Value tree_value = 0;
      for (typename In::const_iterator it = in.begin(); it != in.end(); ++it) {
        if (uf.join(graph.target(it->first),graph.source(it->first))) {
          out.set(it->first, true);
          tree_value += it->second;
        }
        else {
          out.set(it->first, false);
        }
      }
      return tree_value;
    }


    template <typename Sequence>
    struct PairComp {
      typedef typename Sequence::value_type Value;
      bool operator()(const Value& left, const Value& right) {
	return left.second < right.second;
      }
    };

    template <typename In, typename Enable = void>
    struct SequenceInputIndicator {
      static const bool value = false;
    };

    template <typename In>
    struct SequenceInputIndicator<In, 
      typename exists<typename In::value_type::first_type>::type> {
      static const bool value = true;
    };

    template <typename In, typename Enable = void>
    struct MapInputIndicator {
      static const bool value = false;
    };

    template <typename In>
    struct MapInputIndicator<In, 
      typename exists<typename In::Value>::type> {
      static const bool value = true;
    };

    template <typename In, typename Enable = void>
    struct SequenceOutputIndicator {
      static const bool value = false;
    };
 
    template <typename Out>
    struct SequenceOutputIndicator<Out, 
      typename exists<typename Out::value_type>::type> {
      static const bool value = true;
    };

    template <typename Out, typename Enable = void>
    struct MapOutputIndicator {
      static const bool value = false;
    };

    template <typename Out>
    struct MapOutputIndicator<Out, 
      typename exists<typename Out::Value>::type> {
      static const bool value = true;
    };

    template <typename In, typename InEnable = void>
    struct KruskalValueSelector {};

    template <typename In>
    struct KruskalValueSelector<In,
      typename enable_if<SequenceInputIndicator<In>, void>::type> 
    {
      typedef typename In::value_type::second_type Value;
    };    

    template <typename In>
    struct KruskalValueSelector<In,
      typename enable_if<MapInputIndicator<In>, void>::type> 
    {
      typedef typename In::Value Value;
    };    
    
    template <typename Graph, typename In, typename Out,
              typename InEnable = void>
    struct KruskalInputSelector {};

    template <typename Graph, typename In, typename Out,
              typename InEnable = void>
    struct KruskalOutputSelector {};
    
    template <typename Graph, typename In, typename Out>
    struct KruskalInputSelector<Graph, In, Out,
      typename enable_if<SequenceInputIndicator<In>, void>::type > 
    {
      typedef typename In::value_type::second_type Value;

      static Value kruskal(const Graph& graph, const In& in, Out& out) {
        return KruskalOutputSelector<Graph, In, Out>::
          kruskal(graph, in, out);
      }

    };

    template <typename Graph, typename In, typename Out>
    struct KruskalInputSelector<Graph, In, Out,
      typename enable_if<MapInputIndicator<In>, void>::type > 
    {
      typedef typename In::Value Value;
      static Value kruskal(const Graph& graph, const In& in, Out& out) {
        typedef typename In::Key MapEdge;
        typedef typename In::Value Value;
        typedef typename ItemSetTraits<Graph, MapEdge>::ItemIt MapEdgeIt;
        typedef std::vector<std::pair<MapEdge, Value> > Sequence;
        Sequence seq;
        
        for (MapEdgeIt it(graph); it != INVALID; ++it) {
          seq.push_back(std::make_pair(it, in[it]));
        }

        std::sort(seq.begin(), seq.end(), PairComp<Sequence>());
        return KruskalOutputSelector<Graph, Sequence, Out>::
          kruskal(graph, seq, out);
      }
    };

    template <typename Graph, typename In, typename Out>
    struct KruskalOutputSelector<Graph, In, Out,
      typename enable_if<SequenceOutputIndicator<Out>, void>::type > 
    {
      typedef typename In::value_type::second_type Value;

      static Value kruskal(const Graph& graph, const In& in, Out& out) {
        typedef StoreBoolMap<Out> Map;
        Map map(out);
        return _kruskal_bits::kruskal(graph, in, map);
      }

    };

    template <typename Graph, typename In, typename Out>
    struct KruskalOutputSelector<Graph, In, Out,
      typename enable_if<MapOutputIndicator<Out>, void>::type > 
    {
      typedef typename In::value_type::second_type Value;

      static Value kruskal(const Graph& graph, const In& in, Out& out) {
        return _kruskal_bits::kruskal(graph, in, out);
      }
    };

  }

  /// \ingroup spantree
  ///
  /// \brief Kruskal's algorithm to find a minimum cost tree of a graph.
  ///
  /// This function runs Kruskal's algorithm to find a minimum cost tree.
  /// Due to hard C++ hacking, it accepts various input and output types.
  ///
  /// \param g The graph the algorithm runs on.
  /// It can be either \ref concepts::Graph "directed" or 
  /// \ref concepts::UGraph "undirected".
  /// If the graph is directed, the algorithm consider it to be 
  /// undirected by disregarding the direction of the edges.
  ///
  /// \param in This object is used to describe the edge costs. It can be one
  /// of the following choices.
  /// - An STL compatible 'Forward Container' with
  /// <tt>std::pair<GR::UEdge,X></tt> or
  /// <tt>std::pair<GR::Edge,X></tt> as its <tt>value_type</tt>, where
  /// \c X is the type of the costs. The pairs indicates the edges
  /// along with the assigned cost. <em>They must be in a
  /// cost-ascending order.</em>
  /// - Any readable Edge map. The values of the map indicate the edge costs.
  ///
  /// \retval out Here we also have a choise.
  /// - It can be a writable \c bool edge map.  After running the
  /// algorithm this will contain the found minimum cost spanning
  /// tree: the value of an edge will be set to \c true if it belongs
  /// to the tree, otherwise it will be set to \c false. The value of
  /// each edge will be set exactly once.
  /// - It can also be an iteraror of an STL Container with
  /// <tt>GR::UEdge</tt> or <tt>GR::Edge</tt> as its
  /// <tt>value_type</tt>.  The algorithm copies the elements of the
  /// found tree into this sequence.  For example, if we know that the
  /// spanning tree of the graph \c g has say 53 edges, then we can
  /// put its edges into an STL vector \c tree with a code like this.
  ///\code
  /// std::vector<Edge> tree(53);
  /// kruskal(g,cost,tree.begin());
  ///\endcode
  /// Or if we don't know in advance the size of the tree, we can
  /// write this.  
  ///\code std::vector<Edge> tree;
  /// kruskal(g,cost,std::back_inserter(tree)); 
  ///\endcode
  ///
  /// \return The total cost of the found tree.
  ///
  /// \warning If kruskal runs on an be consistent of using the same
  /// Edge type for input and output.
  ///

#ifdef DOXYGEN
  template <class Graph, class In, class Out>
  Value kruskal(GR const& g, const In& in, Out& out)
#else 
  template <class Graph, class In, class Out>
  inline typename _kruskal_bits::KruskalValueSelector<In>::Value 
  kruskal(const Graph& graph, const In& in, Out& out) 
#endif
  {
    return _kruskal_bits::KruskalInputSelector<Graph, In, Out>::
      kruskal(graph, in, out);
  }

 
  

  template <class Graph, class In, class Out>
  inline typename _kruskal_bits::KruskalValueSelector<In>::Value
  kruskal(const Graph& graph, const In& in, const Out& out)
  {
    return _kruskal_bits::KruskalInputSelector<Graph, In, const Out>::
      kruskal(graph, in, out);
  }  

} //namespace lemon

#endif //LEMON_KRUSKAL_H
