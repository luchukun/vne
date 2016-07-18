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

#ifndef LEMON_MIN_MEAN_CYCLE_H
#define LEMON_MIN_MEAN_CYCLE_H

/// \ingroup shortest_path
///
/// \file
/// \brief Howard's algorithm for finding a minimum mean directed cycle.

#include <vector>
#include <lemon/graph_utils.h>
#include <lemon/path.h>
#include <lemon/tolerance.h>
#include <lemon/topology.h>

namespace lemon {

  /// \addtogroup shortest_path
  /// @{

  /// \brief Implementation of Howard's algorithm for finding a minimum
  /// mean directed cycle.
  ///
  /// \ref MinMeanCycle implements Howard's algorithm for finding a
  /// minimum mean directed cycle.
  ///
  /// \tparam Graph The directed graph type the algorithm runs on.
  /// \tparam LengthMap The type of the length (cost) map.
  ///
  /// \warning \c LengthMap::Value must be convertible to \c double.
  ///
  /// \author Peter Kovacs

  template < typename Graph,
             typename LengthMap = typename Graph::template EdgeMap<int> >
  class MinMeanCycle
  {
    GRAPH_TYPEDEFS(typename Graph);

    typedef typename LengthMap::Value Length;
    typedef Path<Graph> Path;

  private:

    // The directed graph the algorithm runs on
    const Graph &_graph;
    // The length of the edges
    const LengthMap &_length;

    // The total length of the found cycle
    Length _cycle_length;
    // The number of edges on the found cycle
    int _cycle_size;
    // The found cycle
    Path *_cycle_path;

    bool _local_path;
    bool _cycle_found;
    Node _cycle_node;

    typename Graph::template NodeMap<bool> _reached;
    typename Graph::template NodeMap<double> _dist;
    typename Graph::template NodeMap<Edge> _policy;

    typename Graph::template NodeMap<int> _component;
    int _component_num;

    std::vector<Node> _nodes;
    std::vector<Edge> _edges;
    Tolerance<double> _tolerance;

  public:

    /// \brief The constructor of the class.
    ///
    /// The constructor of the class.
    ///
    /// \param graph The directed graph the algorithm runs on.
    /// \param length The length (cost) of the edges.
    MinMeanCycle( const Graph &graph,
                  const LengthMap &length ) :
      _graph(graph), _length(length), _cycle_length(0), _cycle_size(-1),
      _cycle_path(NULL), _local_path(false), _reached(graph),
      _dist(graph), _policy(graph), _component(graph)
    {}

    /// The destructor of the class.
    ~MinMeanCycle() {
      if (_local_path) delete _cycle_path;
    }

    /// \brief Sets the \ref Path "path" structure for storing the found
    /// cycle.
    ///
    /// Sets an external \ref Path "path" structure for storing the
    /// found cycle.
    ///
    /// If you don't call this function before calling \ref run() or
    /// \ref init(), it will allocate a local \ref Path "path"
    /// structure.
    /// The destuctor deallocates this automatically allocated map,
    /// of course.
    ///
    /// \note The algorithm calls only the \ref lemon::Path::addBack()
    /// "addBack()" function of the given \ref Path "path" structure.
    ///
    /// \return <tt>(*this)</tt>
    ///
    /// \sa cycle()
    MinMeanCycle& cyclePath(Path &path) {
      if (_local_path) {
        delete _cycle_path;
        _local_path = false;
      }
      _cycle_path = &path;
      return *this;
    }

    /// \name Execution control
    /// The simplest way to execute the algorithm is to call the run()
    /// function.
    /// \n
    /// If you only need the minimum mean value, you may call init()
    /// and findMinMean().
    /// \n
    /// If you would like to run the algorithm again (e.g. the
    /// underlaying graph and/or the edge costs were modified), you may
    /// not create a new instance of the class, rather call reset(),
    /// findMinMean(), and findCycle() instead.

    /// @{

    /// \brief Runs the algorithm.
    ///
    /// Runs the algorithm.
    ///
    /// \return Returns \c true if a directed cycle exists in the graph.
    ///
    /// \note Apart from the return value, <tt>mmc.run()</tt> is just a
    /// shortcut of the following code.
    /// \code
    ///   mmc.init();
    ///   mmc.findMinMean();
    ///   mmc.findCycle();
    /// \endcode
    bool run() {
      init();
      return findMinMean() && findCycle();
    }

    /// \brief Initializes the internal data structures.
    ///
    /// Initializes the internal data structures.
    ///
    /// \sa reset()
    void init() {
      _tolerance.epsilon(1e-6);
      if (!_cycle_path) {
        _local_path = true;
        _cycle_path = new Path;
      }
      _cycle_found = false;
      _component_num = stronglyConnectedComponents(_graph, _component);
    }

    /// \brief Resets the internal data structures.
    ///
    /// Resets the internal data structures so that \ref findMinMean()
    /// and \ref findCycle() can be called again (e.g. when the
    /// underlaying graph has been modified).
    ///
    /// \sa init()
    void reset() {
      if (_cycle_path) _cycle_path->clear();
      _cycle_found = false;
      _component_num = stronglyConnectedComponents(_graph, _component);
    }

    /// \brief Finds the minimum cycle mean length in the graph.
    ///
    /// Computes all the required data and finds the minimum cycle mean
    /// length in the graph.
    ///
    /// \return Returns \c true if a directed cycle exists in the graph.
    ///
    /// \pre \ref init() must be called before using this function.
    bool findMinMean() {
      // Finding the minimum mean cycle in the components
      for (int comp = 0; comp < _component_num; ++comp) {
        if (!initCurrentComponent(comp)) continue;
        while (true) {
          if (!findPolicyCycles()) break;
          contractPolicyGraph(comp);
          if (!computeNodeDistances(comp)) break;
        }
      }
      return _cycle_found;
    }

    /// \brief Finds a critical (minimum mean) directed cycle.
    ///
    /// Finds a critical (minimum mean) directed cycle using the data
    /// computed in the \ref findMinMean() function.
    ///
    /// \return Returns \c true if a directed cycle exists in the graph.
    ///
    /// \pre \ref init() and \ref findMinMean() must be called before
    /// using this function.
    bool findCycle() {
      if (!_cycle_found) return false;
      _cycle_path->addBack(_policy[_cycle_node]);
      for ( Node v = _cycle_node;
            (v = _graph.target(_policy[v])) != _cycle_node; ) {
        _cycle_path->addBack(_policy[v]);
      }
      return true;
    }
    
    /// @}

    /// \name Query Functions
    /// The result of the algorithm can be obtained using these
    /// functions.
    /// \n The algorithm should be executed before using them.

    /// @{
    
    /// \brief Returns the total length of the found cycle.
    ///
    /// Returns the total length of the found cycle.
    ///
    /// \pre \ref run() or \ref findMinMean() must be called before
    /// using this function.
    Length cycleLength() const {
      return _cycle_length;
    }

    /// \brief Returns the number of edges on the found cycle.
    ///
    /// Returns the number of edges on the found cycle.
    ///
    /// \pre \ref run() or \ref findMinMean() must be called before
    /// using this function.
    int cycleEdgeNum() const {
      return _cycle_size;
    }

    /// \brief Returns the mean length of the found cycle.
    ///
    /// Returns the mean length of the found cycle.
    ///
    /// \pre \ref run() or \ref findMinMean() must be called before
    /// using this function.
    ///
    /// \note <tt>mmc.cycleMean()</tt> is just a shortcut of the
    /// following code.
    /// \code
    ///   return double(mmc.cycleLength()) / mmc.cycleEdgeNum();
    /// \endcode
    double cycleMean() const {
      return double(_cycle_length) / _cycle_size;
    }

    /// \brief Returns a const reference to the \ref Path "path"
    /// structure storing the found cycle.
    ///
    /// Returns a const reference to the \ref Path "path"
    /// structure storing the found cycle.
    ///
    /// \pre \ref run() or \ref findCycle() must be called before using
    /// this function.
    ///
    /// \sa cyclePath()
    const Path& cycle() const {
      return *_cycle_path;
    }
    
    ///@}
    
  private:

    // Initializes the internal data structures for the current strongly
    // connected component and creating the policy graph.
    // The policy graph can be represented by the _policy map because
    // the out degree of every node is 1.
    bool initCurrentComponent(int comp) {
      // Finding the nodes of the current component
      _nodes.clear();
      for (NodeIt n(_graph); n != INVALID; ++n) {
        if (_component[n] == comp) _nodes.push_back(n);
      }
      if (_nodes.size() <= 1) return false;
      // Finding the edges of the current component
      _edges.clear();
      for (EdgeIt e(_graph); e != INVALID; ++e) {
        if ( _component[_graph.source(e)] == comp &&
             _component[_graph.target(e)] == comp )
          _edges.push_back(e);
      }
      // Initializing _reached, _dist, _policy maps
      for (int i = 0; i < int(_nodes.size()); ++i) {
        _reached[_nodes[i]] = false;
        _policy[_nodes[i]] = INVALID;
      }
      Node u; Edge e;
      for (int j = 0; j < int(_edges.size()); ++j) {
        e = _edges[j];
        u = _graph.source(e);
        if (!_reached[u] || _length[e] < _dist[u]) {
          _dist[u] = _length[e];
          _policy[u] = e;
          _reached[u] = true;
        }
      }
      return true;
    }

    // Finds all cycles in the policy graph.
    // Sets _cycle_found to true if a cycle is found and sets
    // _cycle_length, _cycle_size, _cycle_node to represent the minimum
    // mean cycle in the policy graph.
    bool findPolicyCycles() {
      typename Graph::template NodeMap<int> level(_graph, -1);
      bool curr_cycle_found = false;
      Length clength;
      int csize;
      int path_cnt = 0;
      Node u, v;
      // Searching for cycles
      for (int i = 0; i < int(_nodes.size()); ++i) {
        if (level[_nodes[i]] < 0) {
          u = _nodes[i];
          level[u] = path_cnt;
          while (level[u = _graph.target(_policy[u])] < 0)
            level[u] = path_cnt;
          if (level[u] == path_cnt) {
            // A cycle is found
            curr_cycle_found = true;
            clength = _length[_policy[u]];
            csize = 1;
            for (v = u; (v = _graph.target(_policy[v])) != u; ) {
              clength += _length[_policy[v]];
              ++csize;
            }
            if ( !_cycle_found ||
                 clength * _cycle_size < _cycle_length * csize ) {
              _cycle_found = true;
              _cycle_length = clength;
              _cycle_size = csize;
              _cycle_node = u;
            }
          }
          ++path_cnt;
        }
      }
      return curr_cycle_found;
    }

    // Contracts the policy graph to be connected by cutting all cycles
    // except for the main cycle (i.e. the minimum mean cycle).
    void contractPolicyGraph(int comp) {
      // Finding the component of the main cycle using
      // reverse BFS search
      typename Graph::template NodeMap<int> found(_graph, false);
      std::deque<Node> queue;
      queue.push_back(_cycle_node);
      found[_cycle_node] = true;
      Node u, v;
      while (!queue.empty()) {
        v = queue.front(); queue.pop_front();
        for (InEdgeIt e(_graph, v); e != INVALID; ++e) {
          u = _graph.source(e);
          if (_component[u] == comp && !found[u] && _policy[u] == e) {
            found[u] = true;
            queue.push_back(u);
          }
        }
      }
      // Connecting all other nodes to this component using
      // reverse BFS search
      queue.clear();
      for (int i = 0; i < int(_nodes.size()); ++i)
        if (found[_nodes[i]]) queue.push_back(_nodes[i]);
      int found_cnt = queue.size();
      while (found_cnt < int(_nodes.size()) && !queue.empty()) {
        v = queue.front(); queue.pop_front();
        for (InEdgeIt e(_graph, v); e != INVALID; ++e) {
          u = _graph.source(e);
          if (_component[u] == comp && !found[u]) {
            found[u] = true;
            ++found_cnt;
            _policy[u] = e;
            queue.push_back(u);
          }
        }
      }
    }

    // Computes node distances in the policy graph and updates the
    // policy graph if the node distances can be improved.
    bool computeNodeDistances(int comp) {
      // Computing node distances using reverse BFS search
      double cycle_mean = double(_cycle_length) / _cycle_size;
      typename Graph::template NodeMap<int> found(_graph, false);
      std::deque<Node> queue;
      queue.push_back(_cycle_node);
      found[_cycle_node] = true;
      _dist[_cycle_node] = 0;
      Node u, v;
      while (!queue.empty()) {
        v = queue.front(); queue.pop_front();
        for (InEdgeIt e(_graph, v); e != INVALID; ++e) {
          u = _graph.source(e);
          if (_component[u] == comp && !found[u] && _policy[u] == e) {
            found[u] = true;
            _dist[u] = _dist[v] + _length[e] - cycle_mean;
            queue.push_back(u);
          }
        }
      }
      // Improving node distances
      bool improved = false;
      for (int j = 0; j < int(_edges.size()); ++j) {
        Edge e = _edges[j];
        u = _graph.source(e); v = _graph.target(e);
        double delta = _dist[v] + _length[e] - cycle_mean;
        if (_tolerance.less(delta, _dist[u])) {
          improved = true;
          _dist[u] = delta;
          _policy[u] = e;
        }
      }
      return improved;
    }

  }; //class MinMeanCycle

  ///@}

} //namespace lemon

#endif //LEMON_MIN_MEAN_CYCLE_H
