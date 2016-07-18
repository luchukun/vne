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

#ifndef LEMON_ITERABLE_MAPS_H
#define LEMON_ITERABLE_MAPS_H

#include <lemon/bits/traits.h>
#include <lemon/bits/invalid.h>

#include <lemon/bits/default_map.h>
#include <lemon/bits/map_extender.h>

#include <vector>
#include <map>

#include <iterator>
#include <limits>

///\ingroup maps
///\file
///\brief Maps that makes it possible to iterate through the keys having
///a certain value
///
///


namespace lemon {

  ///\ingroup graph_maps
  ///
  /// \brief Dynamic iterable bool map.
  ///
  /// This class provides a special graph map type which can store
  /// for each graph item(node, edge, etc.) a bool value. For both 
  /// the true and the false it is possible to iterate on the keys which
  /// mapped to the given value.
  /// 
  /// \param _Graph The graph type.
  /// \param _Item One of the graph's item type, the key of the map.
  template <typename _Graph, typename _Item>
  class IterableBoolMap : protected DefaultMap<_Graph, _Item, int> {
  private:
    typedef _Graph Graph;
    
    typedef typename ItemSetTraits<Graph, _Item>::ItemIt KeyIt;
    typedef DefaultMap<_Graph, _Item, int> Parent;
    
    std::vector<_Item> array;
    int sep;

    const Graph& graph;

  public:

    /// Indicates that the map if reference map.
    typedef True ReferenceMapTag;

    /// The key type
    typedef _Item Key;
    /// The value type
    typedef bool Value;
    /// The const reference type.    
    typedef const Value& ConstReference;

  private:

    int position(const Key& key) const {
      return Parent::operator[](key);
    }

  public:

    /// \brief Refernce to the value of the map.
    ///
    /// This class is near to similar to the bool type. It can
    /// be converted to bool and it has the same operators.
    class Reference {
      friend class IterableBoolMap;
    private:
      Reference(IterableBoolMap& map, const Key& key) 
	: _key(key), _map(map) {} 
    public:

      Reference& operator=(const Reference& value) {
	_map.set(_key, static_cast<bool>(value));
 	return *this;
      }

      operator bool() const { 
	return static_cast<const IterableBoolMap&>(_map)[_key]; 
      }

      Reference& operator=(bool value) { 
	_map.set(_key, value); 
	return *this; 
      }
      Reference& operator&=(bool value) {
	_map.set(_key, _map[_key] & value); 
	return *this; 	
      }
      Reference& operator|=(bool value) {
	_map.set(_key, _map[_key] | value); 
	return *this; 	
      }
      Reference& operator^=(bool value) {
	_map.set(_key, _map[_key] ^ value); 
	return *this; 	
      }
    private:
      Key _key;
      IterableBoolMap& _map; 
    };
    
    /// \brief Constructor of the Map with a default value.
    ///
    /// Constructor of the Map with a default value.
    explicit IterableBoolMap(const Graph& _graph, bool def = false) 
      : Parent(_graph), graph(_graph) {
      for (KeyIt it(graph); it != INVALID; ++it) {
        Parent::set(it, array.size());
        array.push_back(it);
      }
      sep = (def ? array.size() : 0);
    }

    /// \brief Const subscript operator of the map.
    ///
    /// Const subscript operator of the map.
    bool operator[](const Key& key) const {
      return position(key) < sep;
    }

    /// \brief Subscript operator of the map.
    ///
    /// Subscript operator of the map.
    Reference operator[](const Key& key) {
      return Reference(*this, key);
    }

    /// \brief Set operation of the map.
    ///
    /// Set operation of the map.
    void set(const Key& key, bool value) {
      int pos = position(key);
      if (value) {
        if (pos < sep) return;
        Key tmp = array[sep];
        array[sep] = key;
        Parent::set(key, sep);
        array[pos] = tmp;
        Parent::set(tmp, pos); 
        ++sep;
      } else {
        if (pos >= sep) return;
        --sep;
        Key tmp = array[sep];
        array[sep] = key;
        Parent::set(key, sep);
        array[pos] = tmp;
        Parent::set(tmp, pos);
      }
    }

    /// \brief Set all items.
    ///
    /// Set all items in the map.
    /// \note Constant time operation.
    void setAll(bool value) {
      sep = (value ? array.size() : 0);      
    }

    /// \brief Returns the number of the keys mapped to true.
    ///
    /// Returns the number of the keys mapped to true.
    int trueNum() const {
      return sep;
    } 
    
    /// \brief Returns the number of the keys mapped to false.
    ///
    /// Returns the number of the keys mapped to false.
    int falseNum() const {
      return array.size() - sep;
    }

    /// \brief Iterator for the keys mapped to true.
    ///
    /// Iterator for the keys mapped to true. It works
    /// like a graph item iterator in the map, it can be converted
    /// the key type of the map, incremented with \c ++ operator, and
    /// if the iterator leave the last valid key it will be equal to 
    /// \c INVALID.
    class TrueIt : public Key {
    public:
      typedef Key Parent;
      
      /// \brief Creates an iterator.
      ///
      /// Creates an iterator. It iterates on the 
      /// keys which mapped to true.
      /// \param _map The IterableIntMap
      explicit TrueIt(const IterableBoolMap& _map) 
        : Parent(_map.sep > 0 ? _map.array[_map.sep - 1] : INVALID), 
          map(&_map) {}

      /// \brief Invalid constructor \& conversion.
      ///
      /// This constructor initializes the key to be invalid.
      /// \sa Invalid for more details.
      TrueIt(Invalid) : Parent(INVALID), map(0) {}

      /// \brief Increment operator.
      ///
      /// Increment Operator.
      TrueIt& operator++() {
        int pos = map->position(*this);
        Parent::operator=(pos > 0 ? map->array[pos - 1] : INVALID);
        return *this;
      }

      
    private:
      const IterableBoolMap* map;
    };

    /// \brief Iterator for the keys mapped to false.
    ///
    /// Iterator for the keys mapped to false. It works
    /// like a graph item iterator in the map, it can be converted
    /// the key type of the map, incremented with \c ++ operator, and
    /// if the iterator leave the last valid key it will be equal to 
    /// \c INVALID.
    class FalseIt : public Key {
    public:
      typedef Key Parent;
      
      /// \brief Creates an iterator.
      ///
      /// Creates an iterator. It iterates on the 
      /// keys which mapped to false.
      /// \param _map The IterableIntMap
      explicit FalseIt(const IterableBoolMap& _map) 
        : Parent(_map.sep < int(_map.array.size()) ? 
                 _map.array.back() : INVALID), map(&_map) {}

      /// \brief Invalid constructor \& conversion.
      ///
      /// This constructor initializes the key to be invalid.
      /// \sa Invalid for more details.
      FalseIt(Invalid) : Parent(INVALID), map(0) {}

      /// \brief Increment operator.
      ///
      /// Increment Operator.
      FalseIt& operator++() {
        int pos = map->position(*this);
        Parent::operator=(pos > map->sep ? map->array[pos - 1] : INVALID);
        return *this;
      }

    private:
      const IterableBoolMap* map;
    };

    /// \brief Iterator for the keys mapped to a given value.
    ///
    /// Iterator for the keys mapped to a given value. It works
    /// like a graph item iterator in the map, it can be converted
    /// the key type of the map, incremented with \c ++ operator, and
    /// if the iterator leave the last valid key it will be equal to 
    /// \c INVALID.
    class ItemIt : public Key {
    public:
      typedef Key Parent;
      
      /// \brief Creates an iterator.
      ///
      /// Creates an iterator. It iterates on the 
      /// keys which mapped to false.
      /// \param _map The IterableIntMap
      /// \param value Which elements should be iterated.
      ItemIt(const IterableBoolMap& _map, bool value) 
        : Parent(value ? (_map.sep > 0 ? _map.array[_map.sep - 1] : INVALID) :
                 (_map.sep < int(_map.array.size()) ? 
                  _map.array.back() : INVALID)), map(&_map) {}

      /// \brief Invalid constructor \& conversion.
      ///
      /// This constructor initializes the key to be invalid.
      /// \sa Invalid for more details.
      ItemIt(Invalid) : Parent(INVALID), map(0) {}

      /// \brief Increment operator.
      ///
      /// Increment Operator.
      ItemIt& operator++() {
        int pos = map->position(*this);
        int sep = pos >= map->sep ? map->sep : 0;
        Parent::operator=(pos > sep ? map->array[pos - 1] : INVALID);
        return *this;
      }

    private:
      const IterableBoolMap* map;
    };

  protected:
    
    virtual void add(const Key& key) {
      Parent::add(key);
      Parent::set(key, array.size());
      array.push_back(key);
    }

    virtual void add(const std::vector<Key>& keys) {
      Parent::add(keys);
      for (int i = 0; i < int(keys.size()); ++i) {
        Parent::set(keys[i], array.size());
        array.push_back(keys[i]);
      }
    }

    virtual void erase(const Key& key) {
      int pos = position(key); 
      if (pos < sep) {
        --sep;
        Parent::set(array[sep], pos);
        array[pos] = array[sep];
        Parent::set(array.back(), sep);
        array[sep] = array.back();
        array.pop_back();
      } else {
        Parent::set(array.back(), pos);
        array[pos] = array.back();
        array.pop_back();
      }
      Parent::erase(key);
    }

    virtual void erase(const std::vector<Key>& keys) {
      for (int i = 0; i < int(keys.size()); ++i) {
        int pos = position(keys[i]); 
        if (pos < sep) {
          --sep;
          Parent::set(array[sep], pos);
          array[pos] = array[sep];
          Parent::set(array.back(), sep);
          array[sep] = array.back();
          array.pop_back();
        } else {
          Parent::set(array.back(), pos);
          array[pos] = array.back();
          array.pop_back();
        }
      }
      Parent::erase(keys);
    }    

    virtual void build() {
      Parent::build();
      for (KeyIt it(graph); it != INVALID; ++it) {
        Parent::set(it, array.size());
        array.push_back(it);
      }
      sep = 0;      
    }

    virtual void clear() {
      array.clear();
      sep = 0;
      Parent::clear();
    }
    
  };
  

  namespace _iterable_maps_bits {
    template <typename Item>
    struct IterableIntMapNode {
      IterableIntMapNode() : value(-1) {}
      IterableIntMapNode(int _value) : value(_value) {}
      Item prev, next;
      int value;
    };
  }

  ///\ingroup graph_maps
  ///
  /// \brief Dynamic iterable integer map.
  ///
  /// This class provides a special graph map type which can store
  /// for each graph item(node, edge, etc.) an integer value. For each
  /// non negative value it is possible to iterate on the keys which
  /// mapped to the given value.
  /// 
  /// \param _Graph The graph type.
  /// \param _Item One of the graph's item type, the key of the map.
  template <typename _Graph, typename _Item>
  class IterableIntMap 
    : protected MapExtender<DefaultMap<_Graph, _Item, _iterable_maps_bits::
                                       IterableIntMapNode<_Item> > >{
  public:
    typedef MapExtender<DefaultMap<_Graph, _Item, _iterable_maps_bits::
                                   IterableIntMapNode<_Item> > > Parent;

    /// The key type
    typedef _Item Key;
    /// The value type
    typedef int Value;
    /// The graph type
    typedef _Graph Graph;

    /// \brief Constructor of the Map.
    ///
    /// Constructor of the Map. It set all values -1.
    explicit IterableIntMap(const Graph& graph) 
      : Parent(graph) {}

    /// \brief Constructor of the Map with a given value.
    ///
    /// Constructor of the Map with a given value.
    explicit IterableIntMap(const Graph& graph, int value) 
      : Parent(graph, _iterable_maps_bits::IterableIntMapNode<_Item>(value)) {
      if (value >= 0) {
	for (typename Parent::ItemIt it(*this); it != INVALID; ++it) {
	  lace(it);
	}
      }
    }

  private:
    
    void unlace(const Key& key) {
      typename Parent::Value& node = Parent::operator[](key);
      if (node.value < 0) return;
      if (node.prev != INVALID) {
	Parent::operator[](node.prev).next = node.next;	
      } else {
	first[node.value] = node.next;
      }
      if (node.next != INVALID) {
	Parent::operator[](node.next).prev = node.prev;
      }
      while (!first.empty() && first.back() == INVALID) {
	first.pop_back();
      }
    }

    void lace(const Key& key) {
      typename Parent::Value& node = Parent::operator[](key);
      if (node.value < 0) return;
      if (node.value >= int(first.size())) {
	first.resize(node.value + 1, INVALID);
      } 
      node.prev = INVALID;
      node.next = first[node.value];
      if (node.next != INVALID) {
	Parent::operator[](node.next).prev = key;	
      }
      first[node.value] = key;
    }

  public:

    /// Indicates that the map if reference map.
    typedef True ReferenceMapTag;

    /// \brief Refernce to the value of the map.
    ///
    /// This class is near to similar to the int type. It can
    /// be converted to int and it has the same operators.
    class Reference {
      friend class IterableIntMap;
    private:
      Reference(IterableIntMap& map, const Key& key) 
	: _key(key), _map(map) {} 
    public:

      Reference& operator=(const Reference& value) {
	_map.set(_key, static_cast<const int&>(value));
 	return *this;
      }

      operator const int&() const { 
	return static_cast<const IterableIntMap&>(_map)[_key]; 
      }

      Reference& operator=(int value) { 
	_map.set(_key, value); 
	return *this; 
      }
      Reference& operator++() {
	_map.set(_key, _map[_key] + 1); 
	return *this; 	
      }
      int operator++(int) {
	int value = _map[_key];
	_map.set(_key, value + 1); 
	return value; 	
      }
      Reference& operator--() {
	_map.set(_key, _map[_key] - 1); 
	return *this; 	
      }
      int operator--(int) {
	int value = _map[_key];
	_map.set(_key, value - 1); 
	return value; 	
      }
      Reference& operator+=(int value) { 
	_map.set(_key, _map[_key] + value); 
	return *this;
      }
      Reference& operator-=(int value) { 
	_map.set(_key, _map[_key] - value); 
	return *this;
      }
      Reference& operator*=(int value) { 
	_map.set(_key, _map[_key] * value); 
	return *this;
      }
      Reference& operator/=(int value) { 
	_map.set(_key, _map[_key] / value); 
	return *this;
      }
      Reference& operator%=(int value) { 
	_map.set(_key, _map[_key] % value); 
	return *this;
      }
      Reference& operator&=(int value) { 
	_map.set(_key, _map[_key] & value); 
	return *this;
      }
      Reference& operator|=(int value) { 
	_map.set(_key, _map[_key] | value); 
	return *this;
      }
      Reference& operator^=(int value) { 
	_map.set(_key, _map[_key] ^ value); 
	return *this;
      }
      Reference& operator<<=(int value) { 
	_map.set(_key, _map[_key] << value); 
	return *this;
      }
      Reference& operator>>=(int value) { 
	_map.set(_key, _map[_key] >> value); 
	return *this;
      }

    private:
      Key _key;
      IterableIntMap& _map; 
    };

    /// The const reference type.    
    typedef const Value& ConstReference;

    /// \brief Gives back the maximal value plus one.
    ///
    /// Gives back the maximal value plus one.
    unsigned int size() const {
      return first.size();
    }
    
    /// \brief Set operation of the map.
    ///
    /// Set operation of the map.
    void set(const Key& key, const Value& value) {
      unlace(key);
      Parent::operator[](key).value = value;
      lace(key);
    }

    /// \brief Const subscript operator of the map.
    ///
    /// Const subscript operator of the map.
    const Value& operator[](const Key& key) const {
      return Parent::operator[](key).value;
    }

    /// \brief Subscript operator of the map.
    ///
    /// Subscript operator of the map.
    Reference operator[](const Key& key) {
      return Reference(*this, key);
    }

    /// \brief Iterator for the keys with the same value.
    ///
    /// Iterator for the keys with the same value. It works
    /// like a graph item iterator in the map, it can be converted
    /// the item type of the map, incremented with \c ++ operator, and
    /// if the iterator leave the last valid item it will be equal to 
    /// \c INVALID.
    class ItemIt : public _Item {
    public:
      typedef _Item Parent;

      /// \brief Invalid constructor \& conversion.
      ///
      /// This constructor initializes the item to be invalid.
      /// \sa Invalid for more details.
      ItemIt(Invalid) : Parent(INVALID), _map(0) {}

      /// \brief Creates an iterator with a value.
      ///
      /// Creates an iterator with a value. It iterates on the 
      /// keys which have the given value.
      /// \param map The IterableIntMap
      /// \param value The value
      ItemIt(const IterableIntMap& map, int value) : _map(&map) {
	if (value < 0 || value >= int(_map->first.size())) {	  
	  Parent::operator=(INVALID);
	} else {
	  Parent::operator=(_map->first[value]);
	}
      } 

      /// \brief Increment operator.
      ///
      /// Increment Operator.
      ItemIt& operator++() {
	Parent::operator=(_map->IterableIntMap::Parent::
			  operator[](static_cast<Parent&>(*this)).next);
	return *this;
      }


    private:
      const IterableIntMap* _map;
    };

  protected:
    
    virtual void erase(const Key& key) {
      unlace(key);
      Parent::erase(key);
    }

    virtual void erase(const std::vector<Key>& keys) {
      for (int i = 0; i < int(keys.size()); ++i) {
        unlace(keys[i]);
      }
      Parent::erase(keys);
    }

    virtual void clear() {
      first.clear();
      Parent::clear();
    }

  private:
    std::vector<_Item> first;
  };

  namespace _iterable_maps_bits {
    template <typename Item, typename Value>
    struct IterableValueMapNode {
      IterableValueMapNode(Value _value = Value()) : value(_value) {}
      Item prev, next;
      Value value;
    };
  }

  ///\ingroup graph_maps
  ///
  /// \brief Dynamic iterable map for comparable values.
  ///
  /// This class provides a special graph map type which can store
  /// for each graph item(node, edge, etc.) a value. For each
  /// value it is possible to iterate on the keys which mapped to the 
  /// given value. The type stores for each value a linked list with
  /// the items which mapped to the value, and the values are stored
  /// in balanced binary tree. The values of the map can be accessed
  /// with stl compatible forward iterator.
  ///
  /// This type is not reference map so it cannot be modified with
  /// the subscription operator.
  ///
  /// \see InvertableMap
  /// 
  /// \param _Graph The graph type.
  /// \param _Item One of the graph's item type, the key of the map.
  /// \param _Value Any comparable value type.
  template <typename _Graph, typename _Item, typename _Value>
  class IterableValueMap 
    : protected MapExtender<DefaultMap<_Graph, _Item, _iterable_maps_bits::
                                       IterableValueMapNode<_Item, _Value> > >{
  public:
    typedef MapExtender<DefaultMap<_Graph, _Item, _iterable_maps_bits::
                                   IterableValueMapNode<_Item, _Value> > >
    Parent;

    /// The key type
    typedef _Item Key;
    /// The value type
    typedef _Value Value;
    /// The graph type
    typedef _Graph Graph;

  public:

    /// \brief Constructor of the Map with a given value.
    ///
    /// Constructor of the Map with a given value.
    explicit IterableValueMap(const Graph& graph, 
                              const Value& value = Value()) 
      : Parent(graph, _iterable_maps_bits::
               IterableValueMapNode<_Item, _Value>(value)) {
      for (typename Parent::ItemIt it(*this); it != INVALID; ++it) {
        lace(it);
      }
    }

  protected:
    
    void unlace(const Key& key) {
      typename Parent::Value& node = Parent::operator[](key);
      if (node.prev != INVALID) {
	Parent::operator[](node.prev).next = node.next;	
      } else {
        if (node.next != INVALID) {
          first[node.value] = node.next;
        } else {
          first.erase(node.value);
        }
      }
      if (node.next != INVALID) {
	Parent::operator[](node.next).prev = node.prev;
      }
    }

    void lace(const Key& key) {
      typename Parent::Value& node = Parent::operator[](key);
      typename std::map<Value, Key>::iterator it = first.find(node.value);
      if (it == first.end()) {
        node.prev = node.next = INVALID;
        if (node.next != INVALID) {
          Parent::operator[](node.next).prev = key;	
        }
        first.insert(make_pair(node.value, key));
      } else {
        node.prev = INVALID;
        node.next = it->second;
        if (node.next != INVALID) {
          Parent::operator[](node.next).prev = key;	
        }
        it->second = key;
      }
    }

  public:

    /// \brief Forward iterator for values.
    ///
    /// This iterator is an stl compatible forward
    /// iterator on the values of the map. The values can
    /// be accessed in the [beginValue, endValue) range.
    ///
    class ValueIterator 
      : public std::iterator<std::forward_iterator_tag, Value> {
      friend class IterableValueMap;
    private:
      ValueIterator(typename std::map<Value, Key>::const_iterator _it) 
        : it(_it) {}
    public:
      
      ValueIterator() {}

      ValueIterator& operator++() { ++it; return *this; }
      ValueIterator operator++(int) { 
        ValueIterator tmp(*this); 
        operator++();
        return tmp; 
      }

      const Value& operator*() const { return it->first; }
      const Value* operator->() const { return &(it->first); }

      bool operator==(ValueIterator jt) const { return it == jt.it; }
      bool operator!=(ValueIterator jt) const { return it != jt.it; }
      
    private:
      typename std::map<Value, Key>::const_iterator it;
    };

    /// \brief Returns an iterator to the first value.
    ///
    /// Returns an stl compatible iterator to the 
    /// first value of the map. The values of the
    /// map can be accessed in the [beginValue, endValue)
    /// range.
    ValueIterator beginValue() const {
      return ValueIterator(first.begin());
    }

    /// \brief Returns an iterator after the last value.
    ///
    /// Returns an stl compatible iterator after the 
    /// last value of the map. The values of the
    /// map can be accessed in the [beginValue, endValue)
    /// range.
    ValueIterator endValue() const {
      return ValueIterator(first.end());
    }

    /// \brief Set operation of the map.
    ///
    /// Set operation of the map.
    void set(const Key& key, const Value& value) {
      unlace(key);
      Parent::operator[](key).value = value;
      lace(key);
    }

    /// \brief Const subscript operator of the map.
    ///
    /// Const subscript operator of the map.
    const Value& operator[](const Key& key) const {
      return Parent::operator[](key).value;
    }

    /// \brief Iterator for the keys with the same value.
    ///
    /// Iterator for the keys with the same value. It works
    /// like a graph item iterator in the map, it can be converted
    /// the item type of the map, incremented with \c ++ operator, and
    /// if the iterator leave the last valid item it will be equal to 
    /// \c INVALID.
    class ItemIt : public _Item {
    public:
      typedef _Item Parent;

      /// \brief Invalid constructor \& conversion.
      ///
      /// This constructor initializes the item to be invalid.
      /// \sa Invalid for more details.
      ItemIt(Invalid) : Parent(INVALID), _map(0) {}

      /// \brief Creates an iterator with a value.
      ///
      /// Creates an iterator with a value. It iterates on the 
      /// keys which have the given value.
      /// \param map The IterableValueMap
      /// \param value The value
      ItemIt(const IterableValueMap& map, const Value& value) : _map(&map) {
        typename std::map<Value, Key>::const_iterator it = 
          map.first.find(value); 
	if (it == map.first.end()) {	  
	  Parent::operator=(INVALID);
	} else {
	  Parent::operator=(it->second);
	}
      } 

      /// \brief Increment operator.
      ///
      /// Increment Operator.
      ItemIt& operator++() {
	Parent::operator=(_map->IterableValueMap::Parent::
			  operator[](static_cast<Parent&>(*this)).next);
	return *this;
      }


    private:
      const IterableValueMap* _map;
    };

  protected:
    
    virtual void add(const Key& key) {
      Parent::add(key);
      unlace(key);
    }

    virtual void add(const std::vector<Key>& keys) {
      Parent::add(keys);
      for (int i = 0; i < int(keys.size()); ++i) {
        lace(keys[i]);
      }
    }

    virtual void erase(const Key& key) {
      unlace(key);
      Parent::erase(key);
    }

    virtual void erase(const std::vector<Key>& keys) {
      for (int i = 0; i < int(keys.size()); ++i) {
        unlace(keys[i]);
      }
      Parent::erase(keys);
    }

    virtual void build() {
      Parent::build();
      for (typename Parent::ItemIt it(*this); it != INVALID; ++it) {
        lace(it);
      }
    }

    virtual void clear() {
      first.clear();
      Parent::clear();
    }

  private:
    std::map<Value, Key> first;
  };

  /// @}
}

#endif
