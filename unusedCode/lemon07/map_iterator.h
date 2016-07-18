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

#ifndef LEMON_MAP_ITERATOR_H
#define LEMON_MAP_ITERATOR_H

#include <lemon/bits/traits.h>
#include <lemon/bits/utility.h>

/// \ingroup gutils
/// \file
/// \brief Iterators on the maps.

namespace lemon {

  /// \ingroup gutils
  ///
  /// \brief Iterator for maps with possibility of changing values.
  ///
  /// Iterator for maps with possibility of changing values.
  template <typename Graph, typename Item, typename Map>
  class MapIt : public ItemSetTraits<Graph, Item>::ItemIt {
  public:
      
    typedef typename ItemSetTraits<Graph, Item>::ItemIt Parent;
    
    typedef typename Map::Value Value;
    
    /// \brief Creates an iterator
    ///
    /// Creates an iterator for the map, which iterates on the
    /// given graph item set.
    MapIt(const Graph& _graph, Map& _map) : Parent(_graph), map(_map) {}

    /// \brief Gives back the map's value on the current position.
    ///
    /// Gives back the map's value on the current position.
    typename MapTraits<Map>::ConstReturnValue operator*() const {
      return map[*this];
    }

    /// \brief Gives back a reference to the map's value.
    ///
    /// Gives back a reference to the map's value on the current position.
    typename MapTraits<Map>::ReturnValue operator*() {
      return map[*this];
    }
    
    /// \brief Sets the value on the current position
    ///
    /// Sets the value on the current position.
    void set(const Value& value) {
      map.set(*this, value);
    }
    
  protected:
    Map& map;
      
  };

  /// \ingroup gutils
  ///
  /// \brief Iterator for maps with possibility of getting values.
  ///
  /// Iterator for maps with possibility of getting values.
  template <typename Graph, typename Item, typename Map>
  class ConstMapIt : public ItemSetTraits<Graph, Item>::ItemIt {
  public:
    
    typedef typename ItemSetTraits<Graph, Item>::ItemIt Parent;

    typedef typename Map::Value Value;
    
    /// \brief Creates an iterator
    ///
    /// Creates an iterator for the map, which iterates on the
    /// given graph item set.
    ConstMapIt(const Graph& _graph, const Map& _map) 
      : Parent(_graph), map(_map) {}
    
    /// \brief Gives back the map's value on the current position.
    ///
    /// Gives back the map's value on the current position.
    typename MapTraits<Map>::ConstReturnValue operator*() const {
      return map[*this];
    }
    
  protected:
    const Map& map;
  };


  /// \ingroup gutils
  ///
  /// \brief Iterator for maps which filters items by the values.
  ///
  /// Iterator for maps which gives back only that items which mapped
  /// to an given value.
  template <typename Graph, typename Item, typename Map>
  class FilterMapIt 
    : public ItemSetTraits<Graph, Item>::ItemIt {
  public:
    
    typedef typename ItemSetTraits<Graph, Item>::ItemIt Parent;

    typedef typename Map::Value Value;
    
    /// \brief Creates an iterator
    ///
    /// Creates an iterator for the map, which iterates on the
    /// given graph item set and filters all items which mapped value
    /// is not equal to the \c _value.
    FilterMapIt(const Graph& _graph, const Map& _map, const Value& _value) 
      : Parent(_graph), map(_map), value(_value) {}
    
    /// \brief Increment operator
    ///
    /// Skips items which has not mapped to the given value.
    FilterMapIt& operator++() {
      Parent::operator++();
      while (*this != INVALID && map[*this] != value) Parent::operator++();
    }
    
  protected:
    const Map& map;
    Value value;
  };

  
}

#endif
