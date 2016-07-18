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

#ifndef LEMON_PR_BIPARTITE_MATCHING
#define LEMON_PR_BIPARTITE_MATCHING

#include <lemon/graph_utils.h>
#include <lemon/iterable_maps.h>
#include <iostream>
#include <queue>
#include <lemon/elevator.h>

///\ingroup matching
///\file
///\brief Push-prelabel maximum matching algorithms in bipartite graphs.
///
namespace lemon {

  ///Max cardinality matching algorithm based on push-relabel principle

  ///\ingroup matching
  ///Bipartite Max Cardinality Matching algorithm. This class uses the
  ///push-relabel principle which in several cases has better runtime
  ///performance than the augmenting path solutions.
  ///
  ///\author Alpar Juttner
  template<class Graph>
  class PrBipartiteMatching {
    typedef typename Graph::Node Node;
    typedef typename Graph::ANodeIt ANodeIt;
    typedef typename Graph::BNodeIt BNodeIt;
    typedef typename Graph::UEdge UEdge;
    typedef typename Graph::UEdgeIt UEdgeIt;
    typedef typename Graph::IncEdgeIt IncEdgeIt;
    
    const Graph &_g;
    int _node_num;
    int _matching_size;
    int _empty_level;
    
    typename Graph::template ANodeMap<typename Graph::UEdge> _matching;
    Elevator<Graph,typename Graph::BNode> _levels;
    typename Graph::template BNodeMap<int> _cov;

  public:

    /// Constructor

    /// Constructor
    ///
    PrBipartiteMatching(const Graph &g) :
      _g(g),
      _node_num(countBNodes(g)),
      _matching(g),
      _levels(g,_node_num),
      _cov(g,0)
    {
    }
    
    /// \name Execution control 
    /// The simplest way to execute the algorithm is to use one of the
    /// member functions called \c run(). \n
    /// If you need more control on the execution, first
    /// you must call \ref init() and then one variant of the start()
    /// member. 

    /// @{

    ///Initialize the data structures

    ///This function constructs a prematching first, which is a
    ///regular matching on the A-side of the graph, but on the B-side
    ///each node could cover more matching edges. After that, the
    ///B-nodes which multiple matched, will be pushed into the lowest
    ///level of the Elevator. The remaning B-nodes will be pushed to
    ///the consequent levels respect to a Bfs on following graph: the
    ///nodes are the B-nodes of the original bipartite graph and two
    ///nodes are adjacent if a node can pass over a matching edge to
    ///an other node. The source of the Bfs are the lowest level
    ///nodes. Last, the reached B-nodes without covered matching edge
    ///becomes active.
    void init() {
      _matching_size=0;
      _empty_level=_node_num;
      for(ANodeIt n(_g);n!=INVALID;++n)
	if((_matching[n]=IncEdgeIt(_g,n))!=INVALID)
	  ++_cov[_g.bNode(_matching[n])];

      std::queue<Node> q;
      _levels.initStart();
      for(BNodeIt n(_g);n!=INVALID;++n)
	if(_cov[n]>1) {
	  _levels.initAddItem(n);
	  q.push(n);
	}
      int hlev=0;
      while(!q.empty()) {
	Node n=q.front();
	q.pop();
	int nlev=_levels[n]+1;
	for(IncEdgeIt e(_g,n);e!=INVALID;++e) {
	  Node m=_g.runningNode(e);
	  if(e==_matching[m]) {
	    for(IncEdgeIt f(_g,m);f!=INVALID;++f) {
	      Node r=_g.runningNode(f);
	      if(_levels[r]>nlev) {
		for(;nlev>hlev;hlev++)
		  _levels.initNewLevel();
		_levels.initAddItem(r);
		q.push(r);
	      }
	    }
	  }
	}
      }
      _levels.initFinish();
      for(BNodeIt n(_g);n!=INVALID;++n)
	if(_cov[n]<1&&_levels[n]<_node_num)
	  _levels.activate(n);
    }

    ///Start the main phase of the algorithm
    
    ///This algorithm calculates the maximum matching with the
    ///push-relabel principle. This function should be called just
    ///after the init() function which already set the initial
    ///prematching, the level function on the B-nodes and the active,
    ///ie. unmatched B-nodes.
    ///
    ///The algorithm always takes highest active B-node, and it try to
    ///find a B-node which is eligible to pass over one of it's
    ///matching edge. This condition holds when the B-node is one
    ///level lower, and the opposite node of it's matching edge is
    ///adjacent to the highest active node. In this case the current
    ///node steals the matching edge and becomes inactive. If there is
    ///not eligible node then the highest active node should be lift
    ///to the next proper level.
    ///
    ///The nodes should not lift higher than the number of the
    ///B-nodes, if a node reach this level it remains unmatched. If
    ///during the execution one level becomes empty the nodes above it
    ///can be deactivated and lift to the highest level.
    void start() {
      Node act;
      Node bact=INVALID;
      Node last_activated=INVALID;
      while((act=_levels.highestActive())!=INVALID) {
	last_activated=INVALID;
	int actlevel=_levels[act];
	
	UEdge bedge=INVALID;
	int nlevel=_node_num;
	{
	  int nnlevel;
	  for(IncEdgeIt tbedge(_g,act);
	      tbedge!=INVALID && nlevel>=actlevel;
	      ++tbedge)
	    if((nnlevel=_levels[_g.bNode(_matching[_g.runningNode(tbedge)])])<
	       nlevel)
	      {
		nlevel=nnlevel;
		bedge=tbedge;
	      }
	}
	if(nlevel<_node_num) {
	  if(nlevel>=actlevel)
	    _levels.liftHighestActive(nlevel+1);
	  bact=_g.bNode(_matching[_g.aNode(bedge)]);
	  if(--_cov[bact]<1) {
	    _levels.activate(bact);
	    last_activated=bact;
	  }
	  _matching[_g.aNode(bedge)]=bedge;
	  _cov[act]=1;
	  _levels.deactivate(act);
	}
	else {
	  _levels.liftHighestActiveToTop();
	}

	if(_levels.emptyLevel(actlevel))
	  _levels.liftToTop(actlevel);
      }
      
      for(ANodeIt n(_g);n!=INVALID;++n) {
	if (_matching[n]==INVALID)continue;
	if (_cov[_g.bNode(_matching[n])]>1) {
	  _cov[_g.bNode(_matching[n])]--;
	  _matching[n]=INVALID;
	} else {
	  ++_matching_size;
	}
      }
    }

    ///Start the algorithm to find a perfect matching

    ///This function is close to identical to the simple start()
    ///member function but it calculates just perfect matching.
    ///However, the perfect property is only checked on the B-side of
    ///the graph
    ///
    ///The main difference between the two function is the handling of
    ///the empty levels. The simple start() function let the nodes
    ///above the empty levels unmatched while this variant if it find
    ///an empty level immediately terminates and gives back false
    ///return value.
    bool startPerfect() {
      Node act;
      Node bact=INVALID;
      Node last_activated=INVALID;
      while((act=_levels.highestActive())!=INVALID) {
	last_activated=INVALID;
	int actlevel=_levels[act];
	
	UEdge bedge=INVALID;
	int nlevel=_node_num;
	{
	  int nnlevel;
	  for(IncEdgeIt tbedge(_g,act);
	      tbedge!=INVALID && nlevel>=actlevel;
	      ++tbedge)
	    if((nnlevel=_levels[_g.bNode(_matching[_g.runningNode(tbedge)])])<
	       nlevel)
	      {
		nlevel=nnlevel;
		bedge=tbedge;
	      }
	}
	if(nlevel<_node_num) {
	  if(nlevel>=actlevel)
	    _levels.liftHighestActive(nlevel+1);
	  bact=_g.bNode(_matching[_g.aNode(bedge)]);
	  if(--_cov[bact]<1) {
	    _levels.activate(bact);
	    last_activated=bact;
	  }
	  _matching[_g.aNode(bedge)]=bedge;
	  _cov[act]=1;
	  _levels.deactivate(act);
	}
	else {
	  _levels.liftHighestActiveToTop();
	}

	if(_levels.emptyLevel(actlevel))
	  _empty_level=actlevel;
	  return false;
      }
      _matching_size = _node_num;
      return true;
    }
  
    ///Runs the algorithm
    
    ///Just a shortcut for the next code:
    ///\code
    /// init();
    /// start();
    ///\endcode
    void run() {
      init();
      start();
    }
    
    ///Runs the algorithm to find a perfect matching
    
    ///Just a shortcut for the next code:
    ///\code
    /// init();
    /// startPerfect();
    ///\endcode
    ///
    ///\note If the two nodesets of the graph have different size then
    ///this algorithm checks the perfect property on the B-side.
    bool runPerfect() {
      init();
      return startPerfect();
    }

    ///Runs the algorithm to find a perfect matching
    
    ///Just a shortcut for the next code:
    ///\code
    /// init();
    /// startPerfect();
    ///\endcode
    ///
    ///\note It checks that the size of the two nodesets are equal.
    bool checkedRunPerfect() {
      if (countANodes(_g) != _node_num) return false;
      init();
      return startPerfect();
    }

    ///@}

    /// \name Query Functions
    /// The result of the %Matching algorithm can be obtained using these
    /// functions.\n
    /// Before the use of these functions,
    /// either run() or start() must be called.
    ///@{

    ///Set true all matching uedge in the map.

    ///Set true all matching uedge in the map. It does not change the
    ///value mapped to the other uedges.
    ///\return The number of the matching edges.
    template <typename MatchingMap>
    int quickMatching(MatchingMap& mm) const {
      for (ANodeIt n(_g);n!=INVALID;++n) {
        if (_matching[n]!=INVALID) mm.set(_matching[n],true);
      }
      return _matching_size;
    }

    ///Set true all matching uedge in the map and the others to false.

    ///Set true all matching uedge in the map and the others to false.
    ///\return The number of the matching edges.
    template<class MatchingMap>
    int matching(MatchingMap& mm) const {
      for (UEdgeIt e(_g);e!=INVALID;++e) {
        mm.set(e,e==_matching[_g.aNode(e)]);
      }
      return _matching_size;
    }

    ///Gives back the matching in an ANodeMap.

    ///Gives back the matching in an ANodeMap. The parameter should
    ///be a write ANodeMap of UEdge values.
    ///\return The number of the matching edges.
    template<class MatchingMap>
    int aMatching(MatchingMap& mm) const {
      for (ANodeIt n(_g);n!=INVALID;++n) {
        mm.set(n,_matching[n]);
      }
      return _matching_size;
    }

    ///Gives back the matching in a BNodeMap.

    ///Gives back the matching in a BNodeMap. The parameter should
    ///be a write BNodeMap of UEdge values.
    ///\return The number of the matching edges.
    template<class MatchingMap>
    int bMatching(MatchingMap& mm) const {
      for (BNodeIt n(_g);n!=INVALID;++n) {
        mm.set(n,INVALID);
      }
      for (ANodeIt n(_g);n!=INVALID;++n) {
        if (_matching[n]!=INVALID)
	  mm.set(_g.bNode(_matching[n]),_matching[n]);
      }
      return _matching_size;
    }


    ///Returns true if the given uedge is in the matching.

    ///It returns true if the given uedge is in the matching.
    ///
    bool matchingEdge(const UEdge& e) const {
      return _matching[_g.aNode(e)]==e;
    }

    ///Returns the matching edge from the node.

    ///Returns the matching edge from the node. If there is not such
    ///edge it gives back \c INVALID.  
    ///\note If the parameter node is a B-node then the running time is
    ///propotional to the degree of the node.
    UEdge matchingEdge(const Node& n) const {
      if (_g.aNode(n)) {
        return _matching[n];
      } else {
	for (IncEdgeIt e(_g,n);e!=INVALID;++e)
	  if (e==_matching[_g.aNode(e)]) return e;
	return INVALID;
      }
    }

    ///Gives back the number of the matching edges.

    ///Gives back the number of the matching edges.
    int matchingSize() const {
      return _matching_size;
    }

    ///Gives back a barrier on the A-nodes
    
    ///The barrier is s subset of the nodes on the same side of the
    ///graph. If we tried to find a perfect matching and it failed
    ///then the barrier size will be greater than its neighbours. If
    ///the maximum matching searched then the barrier size minus its
    ///neighbours will be exactly the unmatched nodes on the A-side.
    ///\retval bar A WriteMap on the ANodes with bool value.
    template<class BarrierMap>
    void aBarrier(BarrierMap &bar) const 
    {
      for(ANodeIt n(_g);n!=INVALID;++n)
	bar.set(n,_matching[n]==INVALID ||
	  _levels[_g.bNode(_matching[n])]<_empty_level);  
    }  

    ///Gives back a barrier on the B-nodes
    
    ///The barrier is s subset of the nodes on the same side of the
    ///graph. If we tried to find a perfect matching and it failed
    ///then the barrier size will be greater than its neighbours. If
    ///the maximum matching searched then the barrier size minus its
    ///neighbours will be exactly the unmatched nodes on the B-side.
    ///\retval bar A WriteMap on the BNodes with bool value.
    template<class BarrierMap>
    void bBarrier(BarrierMap &bar) const
    {
      for(BNodeIt n(_g);n!=INVALID;++n) bar.set(n,_levels[n]>=_empty_level);  
    }

    ///Returns a minimum covering of the nodes.

    ///The minimum covering set problem is the dual solution of the
    ///maximum bipartite matching. It provides a solution for this
    ///problem what is proof of the optimality of the matching.
    ///\param covering NodeMap of bool values, the nodes of the cover
    ///set will set true while the others false.  
    ///\return The size of the cover set.
    ///\note This function can be called just after the algorithm have
    ///already found a matching. 
    template<class CoverMap>
    int coverSet(CoverMap& covering) const {
      int ret=0;
      for(BNodeIt n(_g);n!=INVALID;++n) {
	if (_levels[n]<_empty_level) { covering.set(n,true); ++ret; }
	else covering.set(n,false);
      }
      for(ANodeIt n(_g);n!=INVALID;++n) {
	if (_matching[n]!=INVALID &&
	    _levels[_g.bNode(_matching[n])]>=_empty_level) 
	  { covering.set(n,true); ++ret; }
	else covering.set(n,false);
      }
      return ret;
    }


    /// @}
    
  };
  
  
  ///Maximum cardinality of the matchings in a bipartite graph

  ///\ingroup matching
  ///This function finds the maximum cardinality of the matchings
  ///in a bipartite graph \c g.
  ///\param g An undirected bipartite graph.
  ///\return The cardinality of the maximum matching.
  ///
  ///\note The the implementation is based
  ///on the push-relabel principle.
  template<class Graph>
  int prBipartiteMatching(const Graph &g)
  {
    PrBipartiteMatching<Graph> bpm(g);
    bpm.run();
    return bpm.matchingSize();
  }

  ///Maximum cardinality matching in a bipartite graph

  ///\ingroup matching
  ///This function finds a maximum cardinality matching
  ///in a bipartite graph \c g.
  ///\param g An undirected bipartite graph.
  ///\retval matching A write ANodeMap of value type \c UEdge.
  /// The found edges will be returned in this map,
  /// i.e. for an \c ANode \c n the edge <tt>matching[n]</tt> is the one
  /// that covers the node \c n.
  ///\return The cardinality of the maximum matching.
  ///
  ///\note The the implementation is based
  ///on the push-relabel principle.
  template<class Graph,class MT>
  int prBipartiteMatching(const Graph &g,MT &matching) 
  {
    PrBipartiteMatching<Graph> bpm(g);
    bpm.run();
    bpm.aMatching(matching);
    return bpm.matchingSize();
  }

  ///Maximum cardinality matching in a bipartite graph

  ///\ingroup matching
  ///This function finds a maximum cardinality matching
  ///in a bipartite graph \c g.
  ///\param g An undirected bipartite graph.
  ///\retval matching A write ANodeMap of value type \c UEdge.
  /// The found edges will be returned in this map,
  /// i.e. for an \c ANode \c n the edge <tt>matching[n]</tt> is the one
  /// that covers the node \c n.
  ///\retval barrier A \c bool WriteMap on the BNodes. The map will be set
  /// exactly once for each BNode. The nodes with \c true value represent
  /// a barrier \e B, i.e. the cardinality of \e B minus the number of its
  /// neighbor is equal to the number of the <tt>BNode</tt>s minus the
  /// cardinality of the maximum matching.
  ///\return The cardinality of the maximum matching.
  ///
  ///\note The the implementation is based
  ///on the push-relabel principle.
  template<class Graph,class MT, class GT>
  int prBipartiteMatching(const Graph &g,MT &matching,GT &barrier) 
  {
    PrBipartiteMatching<Graph> bpm(g);
    bpm.run();
    bpm.aMatching(matching);
    bpm.bBarrier(barrier);
    return bpm.matchingSize();
  }  

  ///Perfect matching in a bipartite graph

  ///\ingroup matching
  ///This function checks whether the bipartite graph \c g
  ///has a perfect matching.
  ///\param g An undirected bipartite graph.
  ///\return \c true iff \c g has a perfect matching.
  ///
  ///\note The the implementation is based
  ///on the push-relabel principle.
  template<class Graph>
  bool prPerfectBipartiteMatching(const Graph &g)
  {
    PrBipartiteMatching<Graph> bpm(g);
    bpm.run();
    return bpm.checkedRunPerfect();
  }

  ///Perfect matching in a bipartite graph

  ///\ingroup matching
  ///This function finds a perfect matching in a bipartite graph \c g.
  ///\param g An undirected bipartite graph.
  ///\retval matching A write ANodeMap of value type \c UEdge.
  /// The found edges will be returned in this map,
  /// i.e. for an \c ANode \c n the edge <tt>matching[n]</tt> is the one
  /// that covers the node \c n.
  /// The values are unchanged if the graph
  /// has no perfect matching.
  ///\return \c true iff \c g has a perfect matching.
  ///
  ///\note The the implementation is based
  ///on the push-relabel principle.
  template<class Graph,class MT>
  bool prPerfectBipartiteMatching(const Graph &g,MT &matching) 
  {
    PrBipartiteMatching<Graph> bpm(g);
    bool ret = bpm.checkedRunPerfect();
    if (ret) bpm.aMatching(matching);
    return ret;
  }

  ///Perfect matching in a bipartite graph

  ///\ingroup matching
  ///This function finds a perfect matching in a bipartite graph \c g.
  ///\param g An undirected bipartite graph.
  ///\retval matching A write ANodeMap of value type \c UEdge.
  /// The found edges will be returned in this map,
  /// i.e. for an \c ANode \c n the edge <tt>matching[n]</tt> is the one
  /// that covers the node \c n.
  /// The values are unchanged if the graph
  /// has no perfect matching.
  ///\retval barrier A \c bool WriteMap on the BNodes. The map will only
  /// be set if \c g has no perfect matching. In this case it is set 
  /// exactly once for each BNode. The nodes with \c true value represent
  /// a barrier, i.e. a subset \e B a of BNodes with the property that
  /// the cardinality of \e B is greater than the number of its neighbors.
  ///\return \c true iff \c g has a perfect matching.
  ///
  ///\note The the implementation is based
  ///on the push-relabel principle.
  template<class Graph,class MT, class GT>
  bool prPerfectBipartiteMatching(const Graph &g,MT &matching,GT &barrier) 
  {
    PrBipartiteMatching<Graph> bpm(g);
    bool ret=bpm.checkedRunPerfect();
    if(ret)
      bpm.aMatching(matching);
    else
      bpm.bBarrier(barrier);
    return ret;
  }  
}

#endif
