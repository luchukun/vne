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

#ifndef LEMON_MATRIX_MAPS_H
#define LEMON_MATRIX_MAPS_H


#include <vector>
#include <lemon/bits/utility.h>
#include <lemon/bits/invalid.h>
#include <lemon/maps.h>

#include <lemon/concepts/matrix_maps.h>

/// \file
/// \ingroup matrices
/// \brief Maps indexed with pairs of items.
///
/// \todo This file has the same name as the concept file in concepts/,
///  and this is not easily detectable in docs...
namespace lemon {

  /// \brief Map for the coloumn view of the matrix
  ///
  /// \ingroup matrices
  /// Map for the coloumn view of the matrix.
  ///
  template <typename _MatrixMap>
  class MatrixRowMap : public MatrixMapTraits<_MatrixMap> {
  public:
    typedef _MatrixMap MatrixMap;
    typedef typename MatrixMap::SecondKey Key;
    typedef typename MatrixMap::Value Value;


    /// \brief Constructor of the row map
    ///
    /// Constructor of the row map.
    MatrixRowMap(MatrixMap& _matrix, typename MatrixMap::FirstKey _row) 
      : matrix(_matrix), row(_row) {}

    /// \brief Subscription operator
    ///
    /// Subscription operator.
    typename MatrixMapTraits<MatrixMap>::ReturnValue
    operator[](Key col) {
      return matrix(row, col);
    }

    /// \brief Setter function
    ///
    /// Setter function.
    void set(Key col, const Value& val) {
      matrix.set(row, col, val);
    }
      
    /// \brief Subscription operator
    ///
    /// Subscription operator.
    typename MatrixMapTraits<MatrixMap>::ConstReturnValue
    operator[](Key col) const {
      return matrix(row, col);
    }

  private:
    MatrixMap& matrix;
    typename MatrixMap::FirstKey row;
  };

  /// \brief Map for the row view of the matrix
  ///
  /// \ingroup matrices
  /// Map for the row view of the matrix.
  ///
  template <typename _MatrixMap>
  class ConstMatrixRowMap : public MatrixMapTraits<_MatrixMap> {
  public:
    typedef _MatrixMap MatrixMap;
    typedef typename MatrixMap::SecondKey Key;
    typedef typename MatrixMap::Value Value;


    /// \brief Constructor of the row map
    ///
    /// Constructor of the row map.
    ConstMatrixRowMap(const MatrixMap& _matrix, 
		      typename MatrixMap::FirstKey _row) 
      : matrix(_matrix), row(_row) {}

    /// \brief Subscription operator
    ///
    /// Subscription operator.
    typename MatrixMapTraits<MatrixMap>::ConstReturnValue
    operator[](Key col) const {
      return matrix(row, col);
    }

  private:
    const MatrixMap& matrix;
    typename MatrixMap::FirstKey row;
  };

  /// \ingroup matrices
  ///
  /// \brief Gives back a row view of the matrix map
  ///
  /// Gives back a row view of the matrix map.
  ///
  /// \sa MatrixRowMap
  /// \sa ConstMatrixRowMap
  template <typename MatrixMap>
  MatrixRowMap<MatrixMap> matrixRowMap(MatrixMap& matrixMap,
				       typename MatrixMap::FirstKey row) {
    return MatrixRowMap<MatrixMap>(matrixMap, row);
  }

  template <typename MatrixMap>
  ConstMatrixRowMap<MatrixMap>
  matrixRowMap(const MatrixMap& matrixMap, typename MatrixMap::FirstKey row) {
    return ConstMatrixRowMap<MatrixMap>(matrixMap, row);
  }

  /// \brief Map for the column view of the matrix
  ///
  /// \ingroup matrices
  /// Map for the column view of the matrix.
  ///
  template <typename _MatrixMap>
  class MatrixColMap : public MatrixMapTraits<_MatrixMap> {
  public:
    typedef _MatrixMap MatrixMap;
    typedef typename MatrixMap::FirstKey Key;
    typedef typename MatrixMap::Value Value;

    /// \brief Constructor of the column map
    ///
    /// Constructor of the column map.
    MatrixColMap(MatrixMap& _matrix, typename MatrixMap::SecondKey _col) 
      : matrix(_matrix), col(_col) {}

    /// \brief Subscription operator
    ///
    /// Subscription operator.
    typename MatrixMapTraits<MatrixMap>::ReturnValue
    operator[](Key row) {
      return matrix(row, col);
    }

    /// \brief Setter function
    ///
    /// Setter function.
    void set(Key row, const Value& val) {
      matrix.set(row, col, val);
    }
      
    /// \brief Subscription operator
    ///
    /// Subscription operator.
    typename MatrixMapTraits<MatrixMap>::ConstReturnValue
    operator[](Key row) const {
      return matrix(row, col);
    }

  private:
    MatrixMap& matrix;
    typename MatrixMap::SecondKey col;
  };

  /// \brief Map for the column view of the matrix
  ///
  /// \ingroup matrices
  /// Map for the column view of the matrix.
  ///
  template <typename _MatrixMap>
  class ConstMatrixColMap : public MatrixMapTraits<_MatrixMap> {
  public:
    typedef _MatrixMap MatrixMap;
    typedef typename MatrixMap::FirstKey Key;
    typedef typename MatrixMap::Value Value;

    /// \brief Constructor of the column map
    ///
    /// Constructor of the column map.
    ConstMatrixColMap(const MatrixMap& _matrix, 
		      typename MatrixMap::SecondKey _col) 
      : matrix(_matrix), col(_col) {}

    /// \brief Subscription operator
    ///
    /// Subscription operator.
    typename MatrixMapTraits<MatrixMap>::ConstReturnValue
    operator[](Key row) const {
      return matrix(row, col);
    }

  private:
    const MatrixMap& matrix;
    typename MatrixMap::SecondKey col;
  };

  /// \ingroup matrices
  ///
  /// \brief Gives back a column view of the matrix map
  ///
  /// Gives back a column view of the matrix map.
  ///
  /// \sa MatrixColMap
  /// \sa ConstMatrixColMap
  template <typename MatrixMap>
  MatrixColMap<MatrixMap> matrixColMap(MatrixMap& matrixMap,
				       typename MatrixMap::SecondKey col) {
    return MatrixColMap<MatrixMap>(matrixMap, col);
  }

  template <typename MatrixMap>
  ConstMatrixColMap<MatrixMap> 
  matrixColMap(const MatrixMap& matrixMap, typename MatrixMap::SecondKey col) {
    return ConstMatrixColMap<MatrixMap>(matrixMap, col);
  }

  /// \brief Container for store values for each ordered pair of graph items
  ///
  /// \ingroup matrices
  /// This data structure can strore for each pair of the same item
  /// type a value. It increase the size of the container when the 
  /// associated graph modified, so it updated automaticly whenever
  /// it is needed.
  template <typename _Graph, typename _Item, typename _Value>
  class DynamicMatrixMap 
    : protected ItemSetTraits<_Graph, _Item>::ItemNotifier::ObserverBase {
  public:
    typedef typename ItemSetTraits<_Graph, _Item>::ItemNotifier::ObserverBase 
    Parent;

    typedef _Graph Graph;
    typedef _Item Key;

    typedef _Item FirstKey;
    typedef _Item SecondKey;
    typedef _Value Value;

    typedef True ReferenceMapTag;

  private:
		
    typedef std::vector<Value> Container;

  public:

    typedef typename Container::reference Reference;
    typedef typename Container::const_reference ConstReference;

    /// \brief Creates an item matrix for the given graph
    ///
    /// Creates an item matrix for the given graph.
    DynamicMatrixMap(const Graph& _graph) 
      : values(size(_graph.maxId(Key()) + 1)) {
      Parent::attach(_graph.notifier(Key()));
    }

    /// \brief Creates an item matrix for the given graph
    ///
    /// Creates an item matrix for the given graph and assigns for each
    /// pairs of keys the given parameter.
    DynamicMatrixMap(const Graph& _graph, const Value& _val) 
      : values(size(_graph.maxId(Key()) + 1), _val) {
      Parent::attach(_graph.notifier(Key()));
    }

    ///\brief The assignement operator.
    ///
    ///It allow to assign a map to an other.
    DynamicMatrixMap& operator=(const DynamicMatrixMap& _cmap){
      return operator=<DynamicMatrixMap>(_cmap);
    }
      
    ///\brief Template assignement operator.
    ///
    ///It copy the element of the given map to its own container.  The
    ///type of the two map shall be the same.
    template <typename CMap>
    DynamicMatrixMap& operator=(const CMap& _cmap){
      checkConcept<concepts::ReadMatrixMap<FirstKey, SecondKey, Value>, CMap>();
      typename Parent::Notifier* notifier = Parent::notifier();
      Key first, second;
      for(notifier->first(first); first != INVALID; 
          notifier->next(first)){
        for(notifier->first(second); second != INVALID; 
            notifier->next(second)){
          set(first, second, _cmap(first, second));
        }
      }
      return *this;
    }

    /// \brief Gives back the value assigned to the \c first - \c second
    /// ordered pair.
    ///
    /// Gives back the value assigned to the \c first - \c second ordered pair.
    ConstReference operator()(const Key& first, const Key& second) const {
      return values[index(Parent::notifier()->id(first), 
                          Parent::notifier()->id(second))];
    }
    
    /// \brief Gives back the value assigned to the \c first - \c second
    /// ordered pair.
    ///
    /// Gives back the value assigned to the \c first - \c second ordered pair.
    Reference operator()(const Key& first, const Key& second) {
      return values[index(Parent::notifier()->id(first), 
                          Parent::notifier()->id(second))];
    }

    /// \brief Setter function for the matrix map.
    ///
    /// Setter function for the matrix map.
    void set(const Key& first, const Key& second, const Value& val) {
      values[index(Parent::notifier()->id(first), 
                   Parent::notifier()->id(second))] = val;
    }

  protected:

    static int index(int i, int j) {
      if (i < j) {
	return j * j + i;
      } else {
	return i * i + i + j;
      }
    }

    static int size(int s) {
      return s * s;
    }

    virtual void add(const Key& key) {
      if (size(Parent::notifier()->id(key) + 1) >= int(values.size())) {
	values.resize(size(Parent::notifier()->id(key) + 1));	
      }
    }

    virtual void add(const std::vector<Key>& keys) {
      int new_size = 0;
      for (int i = 0; i < int(keys.size()); ++i) {
        if (size(Parent::notifier()->id(keys[i]) + 1) >= new_size) {
          new_size = size(Parent::notifier()->id(keys[i]) + 1);	
        }
      }
      if (new_size > int(values.size())) {
        values.resize(new_size);
      }
    }

    virtual void erase(const Key&) {}

    virtual void erase(const std::vector<Key>&) {}

    virtual void build() {
      values.resize(size(Parent::notifier()->maxId() + 1));
    }

    virtual void clear() {
      values.clear();
    }   
    
  private:
    std::vector<Value> values;
  };

  /// \brief Container for store values for each unordered pair of graph items
  ///
  /// \ingroup matrices
  /// This data structure can strore for each pair of the same item
  /// type a value. It increase the size of the container when the 
  /// associated graph modified, so it updated automaticly whenever
  /// it is needed. 
  template <typename _Graph, typename _Item, typename _Value>
  class DynamicSymMatrixMap 
    : protected ItemSetTraits<_Graph, _Item>::ItemNotifier::ObserverBase {
  public:
    typedef typename ItemSetTraits<_Graph, _Item>::ItemNotifier::ObserverBase 
    Parent;

    typedef _Graph Graph;
    typedef _Item Key;

    typedef _Item FirstKey;
    typedef _Item SecondKey;
    typedef _Value Value;

    typedef True ReferenceMapTag;

  private:
		
    typedef std::vector<Value> Container;

  public:

    typedef typename Container::reference Reference;
    typedef typename Container::const_reference ConstReference;

    /// \brief Creates an item matrix for the given graph
    ///
    /// Creates an item matrix for the given graph.
    DynamicSymMatrixMap(const Graph& _graph) 
      : values(size(_graph.maxId(Key()) + 1)) {
      Parent::attach(_graph.notifier(Key()));
    }

    /// \brief Creates an item matrix for the given graph
    ///
    /// Creates an item matrix for the given graph and assigns for each
    /// pairs of keys the given parameter.
    DynamicSymMatrixMap(const Graph& _graph, const Value& _val) 
      : values(size(_graph.maxId(Key()) + 1), _val) {
      Parent::attach(_graph.notifier(Key()));
    }


    ///\brief The assignement operator.
    ///
    ///It allow to assign a map to an other.
    ///
    DynamicSymMatrixMap& operator=(const DynamicSymMatrixMap& _cmap){
      return operator=<DynamicSymMatrixMap>(_cmap);
    }
      
    ///\brief Template assignement operator.
    ///
    ///It copy the element of the given map to its own container.  The
    ///type of the two map shall be the same.
    template <typename CMap>
    DynamicSymMatrixMap& operator=(const CMap& _cmap){
      checkConcept<concepts::ReadMatrixMap<FirstKey, SecondKey, Value>, CMap>();
      typename Parent::Notifier* notifier = Parent::notifier();
      Key first, second;
      for(notifier->first(first); first != INVALID; 
          notifier->next(first)){
        for(notifier->first(second); second != first; 
            notifier->next(second)){
          set(first, second, _cmap(first, second));
        }
        set(first, first, _cmap(first, first));        
      }
      return *this;
    }

    /// \brief Gives back the value assigned to the \c first - \c second
    /// unordered pair.
    ///
    /// Gives back the value assigned to the \c first - \c second unordered 
    /// pair.
    ConstReference operator()(const Key& first, const Key& second) const {
      return values[index(Parent::notifier()->id(first), 
                          Parent::notifier()->id(second))];
    }
    
    /// \brief Gives back the value assigned to the \c first - \c second
    /// unordered pair.
    ///
    /// Gives back the value assigned to the \c first - \c second unordered 
    /// pair.
    Reference operator()(const Key& first, const Key& second) {
      return values[index(Parent::notifier()->id(first), 
                          Parent::notifier()->id(second))];
    }

    /// \brief Setter function for the matrix map.
    ///
    /// Setter function for the matrix map.
    ///
    void set(const Key& first, const Key& second, const Value& val) {
      values[index(Parent::notifier()->id(first), 
                   Parent::notifier()->id(second))] = val;
    }

  protected:

    static int index(int i, int j) {
      if (i < j) {
	return j * (j + 1) / 2 + i;
      } else {
	return i * (i + 1) / 2 + j;
      }
    }

    static int size(int s) {
      return s * (s + 1) / 2;
    }

    virtual void add(const Key& key) {
      if (size(Parent::notifier()->id(key) + 1) >= int(values.size())) {
	values.resize(size(Parent::notifier()->id(key) + 1));	
      }
    }

    virtual void add(const std::vector<Key>& keys) {
      int new_size = 0;
      for (int i = 0; i < int(keys.size()); ++i) {
        if (size(Parent::notifier()->id(keys[i]) + 1) >= new_size) {
          new_size = size(Parent::notifier()->id(keys[i]) + 1);	
        }
      }
      if (new_size > int(values.size())) {
        values.resize(new_size);
      }
    }

    virtual void erase(const Key&) {}

    virtual void erase(const std::vector<Key>&) {}

    virtual void build() {
      values.resize(size(Parent::notifier()->maxId() + 1));
    }

    virtual void clear() {
      values.clear();
    }   
    
  private:
    std::vector<Value> values;
  };
  
  ///\brief Dynamic Asymmetric Matrix Map.
  ///
  ///\ingroup matrices
  ///Dynamic Asymmetric Matrix Map.  Container for store values for each
  ///ordered pair of containers items.  This data structure can store
  ///data with different key types from different container types. It
  ///increases the size of the container if the linked containers
  ///content change, so it is updated automaticly whenever it is
  ///needed.
  ///
  ///This map meet with the concepts::ReferenceMatrixMap<typename K1,
  ///typename K2, typename V, typename R, typename CR> called as
  ///"ReferenceMatrixMap".
  ///
  ///\warning Do not use this map type when the two item sets are
  ///equal or based on the same item set.
  ///
  ///\param _FirstContainer the desired type of first container. It is
  ///ususally a Graph type, but can be any type with alteration
  ///property.
  ///  
  ///\param _FirstContainerItem the nested type of the
  ///FirstContainer. It is usually a graph item as Node, Edge,
  ///etc. This type will be the FirstKey type.
  ///
  ///\param _SecondContainer the desired type of the second
  ///container. It is usualy a Graph type, but can be any type with
  ///alteration property.
  ///
  ///\param _SecondContainerItem the nested type of the
  ///SecondContainer. It is usually a graph item such as Node, Edge,
  ///UEdge, etc. This type will be the SecondKey type.
  ///
  ///\param _Value the type of the strored values in the container.
  ///
  /// \author Janos Nagy
  template <typename _FirstContainer, typename _FirstContainerItem, 
            typename _SecondContainer, typename _SecondContainerItem, 
            typename _Value>
  class DynamicAsymMatrixMap{
  public:

    ///The first key type.
    typedef _FirstContainerItem FirstKey;
      
    ///The second key type.
    typedef _SecondContainerItem SecondKey;
      
    ///The value type of the map.
    typedef _Value Value;
      
    ///Indicates it is a reference map.
    typedef True ReferenceMapTag;
    
  protected:
      
    ///\brief Proxy class for the first key type.
    ///
    ///The proxy class belongs to the FirstKey type. It is necessary because
    ///if one want use the same conatainer types and same nested types but on
    ///other instances of containers than due to the type equiality of nested
    ///types it requires a proxy mechanism. 
    class FirstKeyProxy 
      : protected 
    ItemSetTraits<_FirstContainer,_FirstContainerItem>::
    ItemNotifier::ObserverBase 
    {
        
    public:

      friend class DynamicAsymMatrixMap;
          
      ///Constructor.
      FirstKeyProxy(DynamicAsymMatrixMap& _map) : _owner(_map) { }
    protected:

      ///\brief Add a new FirstKey to the map.
      ///
      ///It adds a new FirstKey to the map. It is called by the
      ///observer notifier and it is ovverride the add() virtual
      ///member function in the observer base. It will call the
      ///maps addFirstKey() function.
      virtual void add(const FirstKey& _firstKey){
        _owner.addFirstKey(_firstKey);
      }
          
      ///\brief Add more new FirstKey to the map.
      ///
      ///It adds more new FirstKey to the map. It is called by the
      ///observer notifier and it is ovverride the add() virtual
      ///member function in the observer base. It will call the
      ///map's addFirstKeys() function.
      virtual void add(const std::vector<FirstKey>& _firstKeys){
        _owner.addFirstKeys(_firstKeys);
      }
          
      ///\brief Erase a FirstKey from the map.
      ///
      ///Erase a FirstKey from the map. It called by the observer
      ///notifier and it overrides the erase() virtual member
      ///function of the observer base. It will call the map's
      ///eraseFirstKey() function.
      virtual void erase(const FirstKey& _firstKey){
        _owner.eraseFirstKey(_firstKey);
      }
          
      ///\brief Erase more FirstKey from the map.
      ///
      ///Erase more FirstKey from the map. It called by the
      ///observer notifier and it overrides the erase() virtual
      ///member function of the observer base. It will call the
      ///map's eraseFirstKeys() function.
      virtual void erase(const std::vector<FirstKey>& _firstKeys){
        _owner.eraseFirstKeys(_firstKeys);
      }
          
      ///\brief Builds the map.
      ///
      ///It buildes the map. It called by the observer notifier
      ///and it overrides the build() virtual member function of
      ///the observer base.  It will call the map's build()
      ///function.
      virtual void build() {
        _owner.build();
        //_owner.buildFirst();
      }
          
      ///\brief Clear the map.
      ///
      ///It erases all items from the map. It called by the
      ///observer notifier and it overrides the clear() virtual
      ///memeber function of the observer base. It will call the
      ///map's clear() function.
      virtual void clear() {
        _owner.clear();
        //_owner.clearFirst();
      }
    private:
          
      ///The map type for it is linked.
      DynamicAsymMatrixMap& _owner;
    };//END OF FIRSTKEYPROXY
      
      ///\brief Proxy class for the second key type.
      ///
      ///The proxy class belongs to the SecondKey type. It is
      ///necessary because if one want use the same conatainer types
      ///and same nested types but on other instances of containers
      ///than due to the type equiality of nested types it requires a
      ///proxy mechanism.
    class SecondKeyProxy
      : protected 
    ItemSetTraits<_SecondContainer, _SecondContainerItem>::
    ItemNotifier::ObserverBase {
        
    public:

      friend class DynamicAsymMatrixMap;
      ///Constructor.
      SecondKeyProxy(DynamicAsymMatrixMap& _map) : _owner(_map) { }

    protected:
          
      ///\brief Add a new SecondKey to the map.
      ///
      ///It adds a new SecondKey to the map. It is called by the
      ///observer notifier and it is ovverride the add() virtual
      ///member function in the observer base. It will call the
      ///maps addSecondKey() function.
      virtual void add(const SecondKey& _secondKey){
        _owner.addSecondKey(_secondKey);
      }
    
      ///\brief Add more new SecondKey to the map.
      ///
      ///It adds more new SecondKey to the map. It is called by
      ///the observer notifier and it is ovverride the add()
      ///virtual member function in the observer base. It will
      ///call the maps addSecondKeys() function.
      virtual void add(const std::vector<SecondKey>& _secondKeys){
        _owner.addSecondKeys(_secondKeys);
      }
          
      ///\brief Erase a SecondKey from the map.
      ///
      ///Erase a SecondKey from the map. It called by the observer
      ///notifier and it overrides the erase() virtual member
      ///function of the observer base. It will call the map's
      ///eraseSecondKey() function.
      virtual void erase(const SecondKey& _secondKey){
        _owner.eraseSecondKey(_secondKey);
      }
          
      ///\brief Erase more SecondKeys from the map.
      ///
      ///Erase more SecondKey from the map. It called by the
      ///observer notifier and it overrides the erase() virtual
      ///member function of the observer base. It will call the
      ///map's eraseSecondKeys() function.
      virtual void erase(const std::vector<SecondKey>& _secondKeys){
        _owner.eraseSecondKeys(_secondKeys);
      }
          
      ///\brief Builds the map.
      ///
      ///It buildes the map. It called by the observer notifier
      ///and it overrides the build() virtual member function of
      ///the observer base.  It will call the map's build()
      ///function.
      virtual void build() {
        _owner.build();
      }
          
      ///\brief Clear the map.
      ///
      ///It erases all items from the map. It called by the
      ///observer notifier and it overrides the clear() virtual
      ///memeber function of the observer base. It will call the
      ///map's clear() function.
      virtual void clear() {
        _owner.clear();
        //_owner.clearFirst();
      }
    private:
          
      ///The type of map for which it is attached.
      DynamicAsymMatrixMap& _owner;
    };//END OF SECONDKEYPROXY
      
  private:
    
    /// \e
    typedef std::vector<Value> Container;
      
    ///The type of constainer which stores the values of the map.
    typedef std::vector<Container> DContainer;

    ///The std:vector type which contains the data
    DContainer values;
      
    ///Member for the first proxy class
    FirstKeyProxy _first_key_proxy;
      
    ///Member for the second proxy class
    SecondKeyProxy _second_key_proxy;

  public:
    
    ///The refernce type of the map.
    typedef typename Container::reference Reference;
      
    ///The const reference type of the constainer.
    typedef typename Container::const_reference ConstReference;

    ///\brief Constructor what create the map for the two containers type.
    ///
    ///Creates the matrix map and initialize the values with Value()
    DynamicAsymMatrixMap(const _FirstContainer& _firstContainer, 
                  const _SecondContainer& _secondContainer)
      : values(DContainer(_firstContainer.maxId(FirstKey())+1,
                          Container(_secondContainer.maxId(SecondKey())+1))),
        _first_key_proxy(*this),
        _second_key_proxy(*this)
    {
      _first_key_proxy.attach(_firstContainer.notifier(FirstKey()));
      _second_key_proxy.attach(_secondContainer.notifier(SecondKey()));
    }

    ///\brief Constructor what create the map for the two containers type.
    ///
    ///Creates the matrix map and initialize the values with the given _value
    DynamicAsymMatrixMap(const _FirstContainer& _firstContainer, 
                  const _SecondContainer& _secondContainer, 
                  const Value& _value)
      : values(DContainer(_firstContainer.maxId(FirstKey())+1,
                          Container(_secondContainer.maxId(SecondKey())+1,
                                    _value))),
        _first_key_proxy(*this),
        _second_key_proxy(*this)
    {
      _first_key_proxy.attach(_firstContainer.notifier(FirstKey()));
      _second_key_proxy.attach(_secondContainer.notifier(SecondKey()));
    }
      
    ///\brief Copy constructor.
    ///
    ///The copy constructor of the map.
    DynamicAsymMatrixMap(const DynamicAsymMatrixMap& _copy) 
      : _first_key_proxy(*this), _second_key_proxy(*this) {
      if(_copy._first_key_proxy.attached() && 
         _copy._second_key_proxy.attached()){
        _first_key_proxy.attach(*_copy._first_key_proxy.notifier());
        _second_key_proxy.attach(*_copy._second_key_proxy.notifier());
        values = _copy.values;
      }
    }
      
    ///\brief Destructor
    ///
    ///Destructor what detach() from the attached objects.  May this
    ///function is not necessary because the destructor of
    ///ObserverBase do the same.
    ~DynamicAsymMatrixMap() {
      if(_first_key_proxy.attached()){
        _first_key_proxy.detach();
      }
      if(_second_key_proxy.attached()){
        _second_key_proxy.detach();
      }
    }
      
    ///\brief Gives back the value assigned to the \c first - \c
    ///second ordered pair.
    ///
    ///Gives back the value assigned to the \c first - \c second
    ///ordered pair.
    Reference operator()(const FirstKey& _first, const SecondKey& _second) {
      return values[_first_key_proxy.notifier()->id(_first)]
        [_second_key_proxy.notifier()->id(_second)];
    }

    ///\brief Gives back the value assigned to the \c first - \c
    ///second ordered pair.
    ///
    ///Gives back the value assigned to the \c first - \c second
    ///ordered pair.
    ConstReference operator()(const FirstKey& _first, 
                              const SecondKey& _second) const {
      return values[_first_key_proxy.notifier()->id(_first)]
        [_second_key_proxy.notifier()->id(_second)];
    }

    ///\brief Setter function for this matrix map.
    ///
    ///Setter function for this matrix map.
    void set(const FirstKey& first, const SecondKey& second, 
             const Value& value){
      values[_first_key_proxy.notifier()->id(first)]
        [_second_key_proxy.notifier()->id(second)] = value;
    }

    ///\brief The assignement operator.
    ///
    ///It allow to assign a map to an other. It
    DynamicAsymMatrixMap& operator=(const DynamicAsymMatrixMap& _cmap){
      return operator=<DynamicAsymMatrixMap>(_cmap);
    }
      
    ///\brief Template assignement operator.
    ///
    ///It copy the element of the given map to its own container.  The
    ///type of the two map shall be the same.
    template <typename CMap>
    DynamicAsymMatrixMap& operator=(const CMap& _cdmap){
      checkConcept<concepts::ReadMatrixMap<FirstKey, SecondKey, Value>, CMap>();
      const typename FirstKeyProxy::Notifier* notifierFirstKey = 
        _first_key_proxy.notifier();
      const typename SecondKeyProxy::Notifier* notifierSecondKey = 
        _second_key_proxy.notifier();
      FirstKey itemFirst;
      SecondKey itemSecond;
      for(notifierFirstKey->first(itemFirst); itemFirst != INVALID; 
          notifierFirstKey->next(itemFirst)){
        for(notifierSecondKey->first(itemSecond); itemSecond != INVALID; 
            notifierSecondKey->next(itemSecond)){
          set(itemFirst, itemSecond, _cdmap(itemFirst,itemSecond));
        }
      }
      return *this;
    }
      
  protected:
    
    ///\brief Add a new FirstKey to the map.
    ///
    ///It adds a new FirstKey to the map. It is called by the observer
    ///class belongs to the FirstKey type.
    void addFirstKey(const FirstKey& firstKey) {
      int size = int(values.size());
      if( _first_key_proxy.notifier()->id(firstKey)+1 >= size ){
        values.resize(_first_key_proxy.notifier()->id(firstKey)+1);
        if( int(values[0].size()) != 0 ){
          int innersize = int(values[0].size());
          for(int i = size; i < int(values.size());++i){
            (values[i]).resize(innersize);
          }
        }else if(_second_key_proxy.notifier()->maxId() >= 0){
          int innersize = _second_key_proxy.notifier()->maxId();
          for(int i = 0; i < int(values.size()); ++i){
            values[0].resize(innersize);
          }
        }
      }
    }

    ///\brief Adds more new FirstKeys to the map.
    ///
    ///It adds more new FirstKeys to the map. It called by the
    ///observer class belongs to the FirstKey type.
    void addFirstKeys(const std::vector<FirstKey>& firstKeys){
      int max = values.size() - 1;
      for(int i = 0; i < int(firstKeys.size()); ++i){
        int id = _first_key_proxy.notifier()->id(firstKeys[i]);
        if(max < id){
          max = id;
        }
      }
      int size = int(values.size());
      if(max >= size){
        values.resize(max + 1);
        if( int(values[0].size()) != 0){
          int innersize = int(values[0].size());
          for(int i = size; i < (max + 1); ++i){
            values[i].resize(innersize);
          }
        }else if(_second_key_proxy.notifier()->maxId() >= 0){
          int innersize = _second_key_proxy.notifier()->maxId();
          for(int i = 0; i < int(values.size()); ++i){
            values[i].resize(innersize);
          }
        }
      }
    }

    ///\brief Add a new SecondKey to the map.
    ///
    ///It adds a new SecondKey to the map. It is called by the
    ///observer class belongs to the SecondKey type.
    void addSecondKey(const SecondKey& secondKey) {
      if(values.size() == 0){
        return;
      }
      int id = _second_key_proxy.notifier()->id(secondKey);
      if(id >= int(values[0].size())){
        for(int i = 0; i < int(values.size());++i){
          values[i].resize(id+1);
        }
      }
    }
        
    ///\brief Adds more new SecondKeys to the map.
    ///
    ///It adds more new SecondKeys to the map. It called by the
    ///observer class belongs to the SecondKey type.
    void addSecondKeys(const std::vector<SecondKey>& secondKeys){
      if(values.size() == 0){
        return;
      }
      int max = values[0].size();
      for(int i = 0; i < int(secondKeys.size()); ++i){
        int id = _second_key_proxy.notifier()->id(secondKeys[i]);
        if(max < id){
          max = id;
        }
      }
      if(max > int(values[0].size())){
        for(int i = 0; i < int(values.size()); ++i){
          values[i].resize(max + 1);
        }
      }
    }
    
    ///\brief Erase a FirstKey from the map.
    ///
    ///Erase a FirstKey from the map. It called by the observer
    ///class belongs to the FirstKey type.
    void eraseFirstKey(const FirstKey& first) {
      int id = _first_key_proxy.notifier()->id(first);
      for(int i = 0; i < int(values[id].size()); ++i){
        values[id][i] = Value();
      }
    }
        
    ///\brief Erase more FirstKey from the map.
    ///
    ///Erase more FirstKey from the map. It called by the observer
    ///class belongs to the FirstKey type.
    void eraseFirstKeys(const std::vector<FirstKey>& firstKeys) {
      for(int j = 0; j < int(firstKeys.size()); ++j){
        int id = _first_key_proxy.notifier()->id(firstKeys[j]);
        for(int i = 0; i < int(values[id].size()); ++i){
          values[id][i] = Value();
        }
      }
    }

    ///\brief Erase a SecondKey from the map.
    ///
    ///Erase a SecondKey from the map. It called by the observer class
    ///belongs to the SecondKey type.
    void eraseSecondKey(const SecondKey& second) {
      if(values.size() == 0){
        return;
      }
      int id = _second_key_proxy.notifier()->id(second);
      for(int i = 0; i < int(values.size()); ++i){
        values[i][id] = Value();
      }
    }
        
    ///\brief Erase more SecondKey from the map.
    ///
    ///Erase more SecondKey from the map. It called by the observer
    ///class belongs to the SecondKey type.
    void eraseSecondKeys(const std::vector<SecondKey>& secondKeys) {
      if(values.size() == 0){
        return;
      }
      for(int j = 0; j < int(secondKeys.size()); ++j){
        int id = _second_key_proxy.notifier()->id(secondKeys[j]);
        for(int i = 0; i < int(values.size()); ++i){
          values[i][id] = Value();
        }
      }
    }

    ///\brief Builds the map.
    ///
    ///It buildes the map. It is called by the observer class belongs
    ///to the FirstKey or SecondKey type.
    void build() {
      values.resize(_first_key_proxy.notifier()->maxId());
      for(int i = 0; i< int(values.size()); ++i){
        values[i].resize(_second_key_proxy.notifier()->maxId());
      }
    }
    
    ///\brief Clear the map.
    ///
    ///It erases all items from the map. It is called by the observer class
    ///belongs to the FirstKey or SecondKey type.
    void clear() {
      for(int i = 0; i < int(values.size()); ++i) {
        values[i].clear();
      }
      values.clear();
    }
 
  };



}

#endif
