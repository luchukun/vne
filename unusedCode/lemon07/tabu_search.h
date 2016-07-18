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

#ifndef LEMON_TABU_SEARCH_H
#define LEMON_TABU_SEARCH_H

/// \ingroup metah
/// \file
/// \brief TabuSearch algorithm.
///
/// \author Szabadkai Mark

#include <lemon/bits/utility.h>
#include <lemon/error.h>
#include <lemon/time_measure.h>
#include <functional>
#include <deque>


namespace lemon {

  /// \brief Default Traits for TabuSearch class.
  /// 
  /// This template defines the needed types for the \ref TabuSearch class.
  /// Is main purpos is to simplify the main class's template interface,
  /// but it provides the EdgeIt type, passing to the concrete graph wheter
  /// it is directed or undirected.
#ifdef DOXYGEN
  template< typename GRAPH, typename VALUE, 
            typename HEIGHTMAP, typename BETTER, bool UNDIR >
#else
  template< typename GRAPH, typename VALUE,
            typename HEIGHTMAP = typename GRAPH::template NodeMap<VALUE>,
            typename BETTER = std::less<VALUE>,
            bool UNDIR = UndirectedTagIndicator<GRAPH>::value >
#endif
  struct TabuSearchDefaultTraits {
    typedef  VALUE  Value; 
    typedef  BETTER  Better;

    typedef  GRAPH  Graph;
    typedef  typename GRAPH::Node  Node;
    typedef  HEIGHTMAP  HeightMap;

    typedef  typename GRAPH::IncEdgeIt  EdgeIt;
  };

  template< typename GRAPH, typename VALUE, 
            typename HEIGHTMAP, typename BETTER >
  struct TabuSearchDefaultTraits< GRAPH, VALUE, HEIGHTMAP, BETTER, false > {
    typedef  VALUE  Value;
    typedef  BETTER  Better;

    typedef  GRAPH  Graph;
    typedef  typename GRAPH::Node  Node;
    typedef  HEIGHTMAP  HeightMap;

    typedef  typename GRAPH::OutEdgeIt  EdgeIt;
  };



  /// \brief Policy hierarchy to controll the search algorithm.
  ///
  /// The fallowing template hierarchy offers a clean interface to define own
  /// policies, and combine existing ones.
  template< typename TS >
  struct TabuSearchPolicyConcept {
    void  target( TS *ts ) {}

    void  reset()  {}
    bool  onStep() { return false; }
    bool  onStick() { return false; }
    bool  onImprove( const typename TS::Value &best ) { return false; }
  };

  template< typename TS >
  struct YesPolicy {
    void  target( TS *ts ) {}

    void  reset()  {}
    bool  onStep() { return true; }
    bool  onStick() { return true; }
    bool  onImprove( const typename TS::Value &best ) { return true; }
  };

  template< typename TS >
  struct NoPolicy : public TabuSearchPolicyConcept<TS> {};

  /// \brief Some basic methode, how tow Policies can be combined
  struct PolicyAndCombination {
    static bool  evaluate( const bool r1, const bool r2 ) {
      return r1 && r2;
    }
  };

  struct PolicyOrCombination {
    static bool  evaluate( const bool r1, const bool r2 ) {
      return r1 || r2;
    }
  };

  /// \brief CombinePolicies
  ///
  /// It combines tow policies using the given combination methode (mainly
  /// some of the basic logical methodes) to create a new one.
#ifdef DOXYGEN
  template< template<typename> class CP1, template<typename> class CP2, 
            typename COMBINATION >
#else
  template< template<typename> class CP1, template<typename> class CP2,
            typename COMBINATION = PolicyAndCombination >
#endif
  struct CombinePolicies {
    template< typename TS >
    struct Policy {
      typedef CP1<TS>  Policy1;
      typedef CP2<TS>  Policy2;
      
      Policy1  policy1;
      Policy2  policy2;

      inline Policy() : policy1(), policy2() {}
      inline Policy( const Policy1 &cp1, const Policy2 &cp2 ) 
        : policy1(cp1), policy2(cp2) {}

      void  target( TS *ts ) {
        policy1.target(ts), policy2.target(ts);
      };

      void  reset() {
        policy1.reset(), policy2.reset();
      }

      bool  onStep() {
        return cmb.evaluate( policy1.onStep(), policy2.onStep() );
      }

      bool  onStick() {
        return cmb.evaluate( policy1.onStick(), policy2.onStick() );
      }

      bool  onImprove( const typename TS::Value &best ) {
        return cmb.evaluate( policy1.onImprove(best), 
                             policy2.onImprove(best) );
      }

    private:
      COMBINATION cmb;
    };
  };


  /// \brief IterationPolicy limits the number of iterations and the
  /// number of iterations without improvement
  template< typename TS >
  struct IterationPolicy {
    IterationPolicy() : _it_lim(100000), _noimpr_it_lim(5000) {}
    IterationPolicy( const long int itl, const long int noimpritl )
      : _it_lim(itl), _noimpr_it_lim(noimpritl)
    {}

    void  target( TS *ts ) {}

    void  reset() {
      _it = _noimpr_it = 0;
    }

    bool  onStep() {
      ++_it; ++_noimpr_it;
      return (_it <= _it_lim) && (_noimpr_it <= _noimpr_it_lim);
    }
		
    bool  onStick() {
      return false;
    }

    bool  onImprove( const typename TS::Value &best ) {
      _noimpr_it = 0;
      return true;
    }

    long int  iterationLimit() const {
      return _it_lim;
    }

    void  iterationLimit( const long int itl ) {
      _it_lim = itl;
    }

    long int  noImprovementIterationLimit() const {
      return _noimpr_it_lim;
    }

    void  noImprovementIterationLimit( const long int noimpritl ) {
      _noimpr_it_lim = noimpritl;
    }

  private:
    long int  _it_lim, _noimpr_it_lim;
    long int  _it, _noimpr_it;
  };

  /// \brief HeightPolicy stops the search when a given height is reached or
  /// exceeds
  template< typename TS >
  struct HeightPolicy {
    typedef typename TS::Value  Value;

    HeightPolicy() : _height_lim(), _found(false) {}
    HeightPolicy( const Value &hl ) : _height_lim(hl), _found(false) {}

    void  target( TS *ts ) {}

    void  reset() {
      _found = false;
    }

    bool  onStep() {
      return !_found;
    }

    bool  onStick() {
      return false;
    }

    bool  onImprove( const Value &best ) {
      typename TS::Better  better;
      _found = better(best, _height_lim) || (best == _height_lim);
      return !_found;
    }

    Value  heightLimi() const {
      return _height_lim;
    }

    void  heightLimi( const Value &hl ) {
      _height_lim = hl;
    }

  private:
    Value  _height_lim;
    bool  _found;
  };

  /// \brief TimePolicy limits the time for searching.
  template< typename TS >
  struct TimePolicy {
    TimePolicy() : _time_lim(60.0), _timeisup(false) {}
    TimePolicy( const double tl ) : _time_lim(tl), _timeisup(false) {}

    void  target( TS *ts ) {}

    void  reset() {
      _timeisup = false;
      _t.reset();
    }

    bool  onStep() {
      update();
      return !_timeisup;
    }

    bool  onStick() {
      return false;
    }

    bool  onImprove( const typename TS::Value &best ) {
      update();
      return !_timeisup;
    }

    double timeLimit() const {
      return _time_lim;
    }

    void  setTimeLimit( const double tl ) {
      _time_lim = tl;
      update();
    }

  private:
    lemon::Timer  _t;
    double  _time_lim;
    bool  _timeisup;

    inline void  update() {
      _timeisup = _t.realTime() > _time_lim;
    }
  };



  /// \ingroup metah
  ///
  /// \brief TabuSearch main class
  ///
  /// This class offers the implementation of tabu-search algorithm. The
  /// tabu-serach is a local-search. It starts from a specified point of the
  /// problem's graph representation, and in every step it goes to the localy
  /// best next Node except those in tabu set. The maximum size of this tabu
  /// set defines how many Node will be remembered. The best Node ever found
  /// will also stored, so we wont lose it, even is the search continues.
  /// The class can be used on any kind of Graph and with any kind of Value
  /// with a total-settlement on it.
  ///
  /// \param _Graph The graph type the algorithm runs on.
  /// \param _Value The values' type associated to the nodes.
  /// \param _Policy Controlls the search. Determinates when to stop, or how
  /// manage stuck search. Default value is \ref IterationPolicy .
  /// \param _Traits Collection of needed types. Default value is
  /// \ref TabuSearchDefaultTraits .
  ///
  /// \author Szabadkai Mark
#ifdef DOXYGEN
  template< typename GRAPH, typename VALUE, template<typename> class POLICY, typename TRAITS >
#else
  template< typename GRAPH, typename VALUE,
            template<typename> class POLICY = IterationPolicy,
            typename TRAITS = TabuSearchDefaultTraits<GRAPH, VALUE> >
#endif
  class TabuSearch
  {
  public:

    /// \brief Thrown by setting the size of the tabu-set and the given size
    /// is less than 2.
    class BadParameterError : public lemon::LogicError {
    public:
      virtual const char* what() const throw() {
        return "lemon::TabuSearch::BadParameterError";
      }
    };

    ///Public types
    typedef  TabuSearch<GRAPH,VALUE,POLICY,TRAITS>  SelfType;

    typedef  typename TRAITS::Graph  Graph;
    typedef  typename TRAITS::Node  Node;
    typedef  typename TRAITS::Value  Value;
    typedef  typename TRAITS::HeightMap  HeightMap;
    typedef  typename TRAITS::Better  Better;
    typedef  typename std::deque< Node >::const_iterator  TabuIterator;

    typedef  POLICY<SelfType>  Policy;

  protected:
    typedef  typename TRAITS::EdgeIt  EdgeIt;

    const Graph  &gr;
    const HeightMap  &height;
    /// The tabu set. Teh current node is the first
    std::deque< Node >  tabu;
    /// Maximal tabu size
    unsigned int  mts;
    /// The best Node found
    Node  b;

    Better  better;
    Policy  pol;

  public:
    /// \brief Constructor
    ///
    /// \param graph the graph the algorithm will run on.
    /// \param hm the height map used by the algorithm.
    /// \param tabusz the maximal size of the tabu set. Default value is 3
    /// \param p the Policy controlling the search.
    TabuSearch( const Graph &graph, const HeightMap &hm, 
                const int tabusz = 3, Policy p = Policy() )
      : gr(graph), height(hm), mts(tabusz), pol(p)
    {
      pol.target(this);
    }

    /// \brief Destructor
    ~TabuSearch() {
      pol.target(NULL);
    }

    /// Set/Get the size of the tabu set
    void  tabuSize( const unsigned int size )
    {
      if( size < 2 )
      throw BadParameterError( "Tabu size must be at least 2!" );
      mts = size;
      while( mts < tabu.size() )
      tabu.pop_back();
    }

    unsigned int  tabuSize() const {
      return mts;
    }

    /// Set/Get Policy
    void  policy( Policy p ) {
      pol.target(NULL);
      pol = p;
      pol.target(this);
    }
		
    Policy& policy()  {
      return pol;
    }

    /// \name Execution control
    /// The simplest way to execute the algorithm is to use the member
    /// functions called \c run( 'startnode' ).
    ///@{

    /// \brief Initializes the internal data.
    ///
    /// \param startn The start node where the search begins.
    void  init( const Node startn ) {
      tabu.clear();
      tabu.push_front( startn );
      b = startn;
      pol.reset();
    }

    /// \brief Does one iteration
    ///
    /// If the Policy allows it searches for the best next node, then steps
    /// onto it.
    /// \return %False if one Policy condition wants to stop the search.
    bool  step()
    {
      ///Request premmision from ControllPolicy
      if( !pol.onStep() )
      return false;
	
      ///Find the best next potential node
      Node n; bool found = false;
      for( EdgeIt e(gr,tabu[0]); e != INVALID; ++e )
      {
        Node m = (gr.source(e) == tabu[0]) ? gr.target(e) : gr.source(e);
        bool wrong = false;
        for( int i = 1; i != (signed int)tabu.size(); ++i )
          if( m == tabu[i] ) {
            wrong = true;
            break;
          }
        if( wrong )
          continue;

        if( !found ) {
          n = m;
          found = true;
        } else
          if( better(height[m], height[n]) ) {
            n = m;
          }
      }

      ///Handle stuck search
      if( !found ) {
        return pol.onStick();
      }

      ///Move on...
      tabu.push_front(n);
      while( mts < tabu.size() ) {
        tabu.pop_back();
      }
      if( better(height[n], height[b]) ) {
        b = n;
        if( !pol.onImprove(height[b]) )
        return false;
      }

      return true;
    }

    /// \brief Runs a search while the Policy stops it.
    ///
    /// \param startn The start node where the search begins.
    inline void  run( const Node startn ) {
      std::cin.unsetf( std::ios_base::skipws );
      char c;
      init( startn );
      while( step() )
      std::cin >> c;
      std::cin.setf( std::ios_base::skipws );
    }

    ///@}

    /// \name Query Functions
    /// The result of the TabuSearch algorithm can be obtained using these
    /// functions.\n
    ///@{

    /// \brief The node, the search is standing on.
    inline Node  current() const {
      return tabu[0];
    }

    /// \brief The best node found until now.
    inline Node  best() const {
      return b;
    }

    /// \brief Beginning to iterate on the current tabu set.
    inline TabuIterator  tabu_begin() const {
      return tabu.begin();
    }

    /// \brief Ending to iterate on the current tabu set.
    inline TabuIterator  tabu_end() const {
      return tabu.end();
    }

    ///@}
  };
}
#endif
