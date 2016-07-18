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

#ifndef LEMON_DIMACS_H
#define LEMON_DIMACS_H

#include <iostream>
#include <string>
#include <vector>
#include <lemon/maps.h>
#include <lemon/bits/invalid.h>

/// \ingroup dimacs_group
/// \file
/// \brief DIMACS file format reader.

namespace lemon {

  ///@defgroup dimacs_group DIMACS format
  ///\brief Read and write files in DIMACS format
  ///
  ///Tools to read a graph from or write it to a file in DIMACS format
  ///data
  ///\ingroup io_group

  /// \addtogroup dimacs_group
  /// @{
  
  /// DIMACS min cost flow reader function.
  ///
  /// This function reads a min cost flow instance from DIMACS format,
  /// i.e. from DIMACS files having a line starting with
  /// \code
  ///   p min
  /// \endcode
  /// At the beginning \c g is cleared by \c g.clear(). The supply 
  /// amount of the nodes are written to \c supply (signed). The 
  /// lower bounds, capacities and costs of the edges are written to 
  /// \c lower, \c capacity and \c cost.
  ///
  /// \author Marton Makai and Peter Kovacs
  template <typename Graph, typename LowerMap, 
    typename CapacityMap, typename CostMap, 
    typename SupplyMap>
  void readDimacs( std::istream& is,
		   Graph &g,
		   LowerMap& lower, 
		   CapacityMap& capacity, 
		   CostMap& cost,
		   SupplyMap& supply )
  {
    g.clear();
    std::vector<typename Graph::Node> nodes;
    typename Graph::Edge e;
    std::string problem, str;
    char c;
    int n, m; 
    int i, j;
    typename SupplyMap::Value sup;
    typename CapacityMap::Value low;
    typename CapacityMap::Value cap;
    typename CostMap::Value co;
    while (is >> c) {
      switch (c) {
      case 'c': // comment line
	getline(is, str);
	break;
      case 'p': // problem definition line
	is >> problem >> n >> m;
	getline(is, str);
	if (problem != "min") return;
	nodes.resize(n + 1);
	for (int k = 1; k <= n; ++k) {
	  nodes[k] = g.addNode();
	  supply.set(nodes[k], 0);
	}
	break;
      case 'n': // node definition line
	is >> i >> sup;
	getline(is, str);
	supply.set(nodes[i], sup);
	break;
      case 'a': // edge (arc) definition line
	is >> i >> j >> low >> cap >> co;
	getline(is, str);
	e = g.addEdge(nodes[i], nodes[j]);
	lower.set(e, low);
	if (cap >= 0)
	  capacity.set(e, cap);
	else
	  capacity.set(e, -1);
	cost.set(e, co);
	break;
      }
    }
  }

  /// DIMACS max flow reader function.
  ///
  /// This function reads a max flow instance from DIMACS format,
  /// i.e. from DIMACS files having a line starting with
  /// \code
  ///   p max
  /// \endcode
  /// At the beginning \c g is cleared by \c g.clear(). The edge 
  /// capacities are written to \c capacity and \c s and \c t are
  /// set to the source and the target nodes.
  ///
  /// \author Marton Makai
  template<typename Graph, typename CapacityMap>
  void readDimacs(std::istream& is, Graph &g, CapacityMap& capacity, 
		  typename Graph::Node &s, typename Graph::Node &t) {
    g.clear();
    std::vector<typename Graph::Node> nodes;
    typename Graph::Edge e;
    std::string problem;
    char c, d;
    int n, m; 
    int i, j;
    typename CapacityMap::Value _cap;
    std::string str;
    while (is >> c) {
      switch (c) {
      case 'c': // comment line
	getline(is, str);
	break;
      case 'p': // problem definition line
	is >> problem >> n >> m;
	getline(is, str);
	nodes.resize(n + 1);
	for (int k = 1; k <= n; ++k)
	  nodes[k] = g.addNode();
	break;
      case 'n': // node definition line
	if (problem == "sp") { // shortest path problem
	  is >> i;
	  getline(is, str);
	  s = nodes[i];
	}
	if (problem == "max") { // max flow problem
	  is >> i >> d;
	  getline(is, str);
	  if (d == 's') s = nodes[i];
	  if (d == 't') t = nodes[i];
	}
	break;
      case 'a': // edge (arc) definition line
	if (problem == "max" || problem == "sp") {
	  is >> i >> j >> _cap;
	  getline(is, str);
	  e = g.addEdge(nodes[i], nodes[j]);
	  capacity.set(e, _cap);
	} else {
	  is >> i >> j;
	  getline(is, str);
	  g.addEdge(nodes[i], nodes[j]);
	}
	break;
      }
    }
  }

  /// DIMACS shortest path reader function.
  ///
  /// This function reads a shortest path instance from DIMACS format,
  /// i.e. from DIMACS files having a line starting with
  /// \code
  ///   p sp
  /// \endcode
  /// At the beginning \c g is cleared by \c g.clear(). The edge
  /// capacities are written to \c capacity and \c s is set to the
  /// source node.
  ///
  /// \author Marton Makai
  template<typename Graph, typename CapacityMap>
  void readDimacs(std::istream& is, Graph &g, CapacityMap& capacity, 
		  typename Graph::Node &s) {
    readDimacs(is, g, capacity, s, s);
  }

  /// DIMACS capacitated graph reader function.
  ///
  /// This function reads an edge capacitated graph instance from
  /// DIMACS format. At the beginning \c g is cleared by \c g.clear()
  /// and the edge capacities are written to \c capacity.
  ///
  /// \author Marton Makai
  template<typename Graph, typename CapacityMap>
  void readDimacs(std::istream& is, Graph &g, CapacityMap& capacity) {
    typename Graph::Node u;
    readDimacs(is, g, capacity, u, u);
  }

  /// DIMACS plain graph reader function.
  ///
  /// This function reads a graph without any designated nodes and
  /// maps from DIMACS format, i.e. from DIMACS files having a line
  /// starting with
  /// \code
  ///   p mat
  /// \endcode
  /// At the beginning \c g is cleared by \c g.clear().
  ///
  /// \author Marton Makai
  template<typename Graph>
  void readDimacs(std::istream& is, Graph &g) {
    typename Graph::Node u;
    NullMap<typename Graph::Edge, int> n;
    readDimacs(is, g, n, u, u);
  }
  
  /// DIMACS plain graph writer function.
  ///
  /// This function writes a graph without any designated nodes and
  /// maps into DIMACS format, i.e. into DIMACS file having a line 
  /// starting with
  /// \code
  ///   p mat
  /// \endcode
  ///
  /// \author Marton Makai
  template<typename Graph>
  void writeDimacs(std::ostream& os, const Graph &g) {
    typedef typename Graph::NodeIt NodeIt;
    typedef typename Graph::EdgeIt EdgeIt;  
    
    os << "c matching problem" << std::endl;
    os << "p mat " << g.nodeNum() << " " << g.edgeNum() << std::endl;

    typename Graph::template NodeMap<int> nodes(g);
    int i = 1;
    for(NodeIt v(g); v != INVALID; ++v) {
      nodes.set(v, i);
      ++i;
    }    
    for(EdgeIt e(g); e != INVALID; ++e) {
      os << "a " << nodes[g.source(e)] << " " << nodes[g.target(e)] << std::endl; 
    }
  }

  /// @}

} //namespace lemon

#endif //LEMON_DIMACS_H
