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

#include<lemon/bits/invalid.h>
#include<lemon/topology.h>
#include <list>

/// \ingroup graph_prop
/// \file
/// \brief Euler tour
///
///This file provides an Euler tour iterator and ways to check
///if a graph is euler.


namespace lemon {

  ///Euler iterator for directed graphs.

  /// \ingroup graph_prop
  ///This iterator converts to the \c Edge type of the graph and using
  ///operator ++ it provides an Euler tour of a \e directed
  ///graph (if there exists).
  ///
  ///For example
  ///if the given graph if Euler (i.e it has only one nontrivial component
  ///and the in-degree is equal to the out-degree for all nodes),
  ///the following code will put the edges of \c g
  ///to the vector \c et according to an
  ///Euler tour of \c g.
  ///\code
  ///  std::vector<ListGraph::Edge> et;
  ///  for(EulerIt<ListGraph> e(g),e!=INVALID;++e)
  ///    et.push_back(e);
  ///\endcode
  ///If \c g is not Euler then the resulted tour will not be full or closed.
  ///\sa UEulerIt
  ///\todo Test required
  template<class Graph>
  class EulerIt 
  {
    typedef typename Graph::Node Node;
    typedef typename Graph::NodeIt NodeIt;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::EdgeIt EdgeIt;
    typedef typename Graph::OutEdgeIt OutEdgeIt;
    typedef typename Graph::InEdgeIt InEdgeIt;
    
    const Graph &g;
    typename Graph::template NodeMap<OutEdgeIt> nedge;
    std::list<Edge> euler;

  public:
    
    ///Constructor

    ///\param _g A directed graph.
    ///\param start The starting point of the tour. If it is not given
    ///       the tour will start from the first node.
    EulerIt(const Graph &_g,typename Graph::Node start=INVALID)
      : g(_g), nedge(g)
    {
      if(start==INVALID) start=NodeIt(g);
      for(NodeIt n(g);n!=INVALID;++n) nedge[n]=OutEdgeIt(g,n);
      while(nedge[start]!=INVALID) {
	euler.push_back(nedge[start]);
	Node next=g.target(nedge[start]);
	++nedge[start];
	start=next;
      }
    }
    
    ///Edge Conversion
    operator Edge() { return euler.empty()?INVALID:euler.front(); }
    bool operator==(Invalid) { return euler.empty(); }
    bool operator!=(Invalid) { return !euler.empty(); }
    
    ///Next edge of the tour
    EulerIt &operator++() {
      Node s=g.target(euler.front());
      euler.pop_front();
      //This produces a warning.Strange.
      //std::list<Edge>::iterator next=euler.begin();
      typename std::list<Edge>::iterator next=euler.begin();
      while(nedge[s]!=INVALID) {
	euler.insert(next,nedge[s]);
	Node n=g.target(nedge[s]);
	++nedge[s];
	s=n;
      }
      return *this;
    }
    ///Postfix incrementation
    
    ///\warning This incrementation
    ///returns an \c Edge, not an \ref EulerIt, as one may
    ///expect.
    Edge operator++(int) 
    {
      Edge e=*this;
      ++(*this);
      return e;
    }
  };

  ///Euler iterator for undirected graphs.

  /// \ingroup graph_prop
  ///This iterator converts to the \c Edge (or \c UEdge)
  ///type of the graph and using
  ///operator ++ it provides an Euler tour of an undirected
  ///graph (if there exists).
  ///
  ///For example
  ///if the given graph if Euler (i.e it has only one nontrivial component
  ///and the degree of each node is even),
  ///the following code will print the edge IDs according to an
  ///Euler tour of \c g.
  ///\code
  ///  for(UEulerIt<ListUGraph> e(g),e!=INVALID;++e) {
  ///    std::cout << g.id(UEdge(e)) << std::eol;
  ///  }
  ///\endcode
  ///Although the iterator provides an Euler tour of an undirected graph,
  ///in order to indicate the direction of the tour, UEulerIt
  ///returns directed edges (that convert to the undirected ones, of course).
  ///
  ///If \c g is not Euler then the resulted tour will not be full or closed.
  ///\sa EulerIt
  ///\todo Test required
  template<class Graph>
  class UEulerIt
  {
    typedef typename Graph::Node Node;
    typedef typename Graph::NodeIt NodeIt;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::UEdge UEdge;
    typedef typename Graph::EdgeIt EdgeIt;
    typedef typename Graph::OutEdgeIt OutEdgeIt;
    typedef typename Graph::InEdgeIt InEdgeIt;
    
    const Graph &g;
    typename Graph::template NodeMap<OutEdgeIt> nedge;
    typename Graph::template UEdgeMap<bool> visited;
    std::list<Edge> euler;

  public:
    
    ///Constructor

    ///\param _g An undirected graph.
    ///\param start The starting point of the tour. If it is not given
    ///       the tour will start from the first node.
    UEulerIt(const Graph &_g,typename Graph::Node start=INVALID)
      : g(_g), nedge(g), visited(g,false)
    {
      if(start==INVALID) start=NodeIt(g);
      for(NodeIt n(g);n!=INVALID;++n) nedge[n]=OutEdgeIt(g,n);
      while(nedge[start]!=INVALID) {
	euler.push_back(nedge[start]);
	visited[nedge[start]]=true;
	Node next=g.target(nedge[start]);
	++nedge[start];
	start=next;
	while(nedge[start]!=INVALID && visited[nedge[start]]) ++nedge[start];
      }
    }
    
    ///Edge Conversion
    operator Edge() const { return euler.empty()?INVALID:euler.front(); }
    ///Edge Conversion
    operator UEdge() const { return euler.empty()?INVALID:euler.front(); }
    ///\e
    bool operator==(Invalid) const { return euler.empty(); }
    ///\e
    bool operator!=(Invalid) const { return !euler.empty(); }
    
    ///Next edge of the tour
    UEulerIt &operator++() {
      Node s=g.target(euler.front());
      euler.pop_front();
      typename std::list<Edge>::iterator next=euler.begin();

      while(nedge[s]!=INVALID) {
	while(nedge[s]!=INVALID && visited[nedge[s]]) ++nedge[s];
	if(nedge[s]==INVALID) break;
	else {
	  euler.insert(next,nedge[s]);
	  visited[nedge[s]]=true;
	  Node n=g.target(nedge[s]);
	  ++nedge[s];
	  s=n;
	}
      }
      return *this;
    }
    
    ///Postfix incrementation
    
    ///\warning This incrementation
    ///returns an \c Edge, not an \ref UEulerIt, as one may
    ///expect.
    Edge operator++(int) 
    {
      Edge e=*this;
      ++(*this);
      return e;
    }
  };


  ///Checks if the graph is Euler

  /// \ingroup graph_prop
  ///Checks if the graph is Euler. It works for both directed and
  ///undirected graphs.
  ///\note By definition, a directed graph is called \e Euler if
  ///and only if connected and the number of it is incoming and outgoing edges
  ///are the same for each node.
  ///Similarly, an undirected graph is called \e Euler if
  ///and only if it is connected and the number of incident edges is even
  ///for each node. <em>Therefore, there are graphs which are not Euler, but
  ///still have an Euler tour</em>.
  ///\todo Test required
  template<class Graph>
#ifdef DOXYGEN
  bool
#else
  typename enable_if<UndirectedTagIndicator<Graph>,bool>::type
  euler(const Graph &g) 
  {
    for(typename Graph::NodeIt n(g);n!=INVALID;++n)
      if(countIncEdges(g,n)%2) return false;
    return connected(g);
  }
  template<class Graph>
  typename disable_if<UndirectedTagIndicator<Graph>,bool>::type
#endif
  euler(const Graph &g) 
  {
    for(typename Graph::NodeIt n(g);n!=INVALID;++n)
      if(countInEdges(g,n)!=countOutEdges(g,n)) return false;
    return connected(g);
  }
  
}
