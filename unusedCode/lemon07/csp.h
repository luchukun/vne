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

#ifndef LEMON_CSP_H
#define LEMON_CSP_H

///\ingroup approx
///\file
///\brief Algorithm for the Resource Constrained Shortest Path problem.
///

#include <lemon/list_graph.h>
#include <lemon/graph_utils.h>
#include <lemon/error.h>
#include <lemon/maps.h>
#include <lemon/tolerance.h>
#include <lemon/dijkstra.h>
#include <lemon/path.h>
#include <lemon/counter.h>
namespace lemon {

  ///\ingroup approx
  ///
  ///\brief Algorithms for the Resource Constrained Shortest Path Problem
  ///
  ///The Resource Constrained Shortest (Least Cost) Path problem is the
  ///following. We are given a directed graph with two additive weightings
  ///on the edges, referred as \e cost and \e delay. In addition,
  ///a source and a destination node \e s and \e t and a delay
  ///constraint \e D is given. A path \e p is called \e feasible
  ///if <em>delay(p)\<=D</em>. Then, the task is to find the least cost
  ///feasible path.
  ///
  template<class Graph,
	   class CM=typename Graph:: template EdgeMap<double>,
	   class DM=CM>
  class ConstrainedShortestPath 
  {
  public:
    
    GRAPH_TYPEDEFS(typename Graph);
    
    typedef SimplePath<Graph> Path;
    
  private:
    
    const Graph &_g;
    Tolerance<double> tol;

    const CM &_cost;
    const DM &_delay;

    class CoMap 
    {
      const CM &_cost;
      const DM &_delay;
      double _lambda;
    public:
      typedef typename CM::Key Key;
      typedef double Value;
      CoMap(const CM &c, const DM &d) :_cost(c), _delay(d), _lambda(0) {}
      double lambda() const { return _lambda; }
      void lambda(double l)  { _lambda=l; }
      Value operator[](Key &e) const 
      {
	return _cost[e]+_lambda*_delay[e];
      }
    };

    CoMap _co_map;
    
    
    Dijkstra<Graph, CoMap> _dij;

  public:
    
    /// \brief Constructor
    
    ///Constructor
    ///
    ConstrainedShortestPath(const Graph &g, const CM &ct, const DM &dl)
      : _g(g), _cost(ct), _delay(dl),
	_co_map(ct, dl), _dij(_g,_co_map) {}
    

    ///Compute the cost of a path
    double cost(const Path &p) const
    {
      double s=0;
      //      Path r;  
      for(typename Path::EdgeIt e(p);e!=INVALID;++e) s+=_cost[e];
      return s;
    }

    ///Compute the delay of a path
    double delay(const Path &p) const
    {
      double s=0;
      for(typename Path::EdgeIt e(p);e!=INVALID;++e) s+=_delay[e];
      return s;
    }
    
    ///Runs the LARAC algorithm
    
    ///This function runs a Lagrange relaxation based heuristic to find
    ///a delay constrained least cost path.
    ///\param s source node
    ///\param t target node
    ///\param delta upper bound on the delta
    ///\retval lo_bo a lower bound on the optimal solution
    ///\return the found path or an empty 
    Path larac(Node s, Node t, double delta, double &lo_bo) 
    {
      double lambda=0;
      double cp,cq,dp,dq,cr,dr;
      Path p;
      Path q;
      Path r;
      {
	Dijkstra<Graph,CM> dij(_g,_cost);
	dij.run(s,t);
	if(!dij.reached(t)) return Path();
	p=dij.path(t);
	cp=cost(p);
	dp=delay(p);
      }
      if(delay(p)<=delta) return p;
      {
	Dijkstra<Graph,DM> dij(_g,_delay);
	dij.run(s,t);
	q=dij.path(t);
	cq=cost(q);
	dq=delay(q);
      }
      if(delay(q)>delta) return Path();
      while (true) {
	lambda=(cp-cq)/(dq-dp);
	_co_map.lambda(lambda);
	_dij.run(s,t);
	r=_dij.path(t);
	cr=cost(r);
	dr=delay(r);
	if(!tol.less(cr+lambda*dr,cp+lambda*dp)) {
	  lo_bo=cq+lambda*(dq-delta);
	  return q;
	}
	else if(tol.less(dr,delta)) 
	  {
	    q=r;
	    cq=cr;
	    dq=dr;
	  }
	else if(tol.less(delta,dr))
	  {
	    p=r;
	    cp=cr;
	    dp=dr;
	  }
	else return r;
      }
    }
  };
  

} //END OF NAMESPACE LEMON

#endif
