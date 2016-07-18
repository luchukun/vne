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

///\ingroup lemon_io
///\file
///\brief Lemon Format reader.


#ifndef LEMON_LEMON_READER_H
#define LEMON_LEMON_READER_H


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>

#include <lemon/error.h>
#include <lemon/graph_utils.h>
#include <lemon/bits/utility.h>
#include <lemon/bits/item_reader.h>

#include <lemon/dim2.h>

#include <lemon/concept_check.h>
#include <lemon/concepts/maps.h>

namespace lemon {

  namespace _reader_bits {

    template <typename T>
    bool operator<(T, T) {
      throw DataFormatError("Label is not comparable");
    }

    template <typename T>
    struct Less {
      bool operator()(const T& p, const T& q) const {
	return p < q;
      }
    };
    
    template <typename Item>
    class ItemLabelReader {
    public:

      bool isLabelReader() { return true; }

      void readLabel(std::istream&, Item&) {}
      
      template <class _ItemLabelReader>
      struct Constraints {
	void constraints() {
	  bool b = reader.isLabelReader();
	  ignore_unused_variable_warning(b);
	  Item item;
	  reader.readLabel(is, item);
	}
	_ItemLabelReader& reader;
	std::istream& is;
      };

    };

    template <typename Item>
    class ItemReader {
    public:
      void read(std::istream&, Item&) {}
      
      template <class _ItemReader>
      struct Constraints {
	void constraints() {
	  Item item;
	  reader.read(is, item);
	}
	_ItemReader& reader;
	std::istream& is;
      };

    };

    template <typename Map>
    struct Ref { typedef Map& Type; };
    template <typename Map>
    struct Arg { typedef Map& Type; };

    template <typename Graph, typename Map>
    class ForwardComposeMap {
    public:
      typedef typename Graph::UEdge Key;
      typedef typename Map::Value Value;

      ForwardComposeMap(const Graph& _graph, typename Arg<Map>::Type _map) 
	: graph(_graph), map(_map) {}
      
      void set(const Key& key, const Value& val) {
	map.set(graph.direct(key, true), val);
      }

    private:
      const Graph& graph;
      typename Ref<Map>::Type map;
    };

    template <typename Graph, typename Map>
    ForwardComposeMap<Graph, Map>
    forwardComposeMap(const Graph& graph, const Map& map) {
      return ForwardComposeMap<Graph, Map>(graph, map);
    }

    template <typename Graph, typename Map>
    ForwardComposeMap<Graph, Map>
    forwardComposeMap(const Graph& graph, Map& map) {
      return ForwardComposeMap<Graph, Map>(graph, map);
    }

    template <typename Graph, typename Map>
    class BackwardComposeMap {
    public:
      typedef typename Graph::UEdge Key;
      typedef typename Map::Value Value;

      BackwardComposeMap(const Graph& _graph, typename Arg<Map>::Type _map) 
	: graph(_graph), map(_map) {}
      
      void set(const Key& key, const Value& val) {
	map.set(graph.direct(key, false), val);
      }

    private:
      const Graph& graph;
      typename Ref<Map>::Type map;
    };


    template <typename Graph, typename Map>
    BackwardComposeMap<Graph, Map>
    backwardComposeMap(const Graph& graph, const Map& map) {
      return BackwardComposeMap<Graph, Map>(graph, map);
    }

    template <typename Graph, typename Map>
    BackwardComposeMap<Graph, Map>
    backwardComposeMap(const Graph& graph, Map& map) {
      return BackwardComposeMap<Graph, Map>(graph, map);
    }

    template <typename Graph, typename Map>
    struct Ref<ForwardComposeMap<Graph, Map> > { 
      typedef ForwardComposeMap<Graph, Map> Type;
    };
    template <typename Graph, typename Map>
    struct Arg<ForwardComposeMap<Graph, Map> > { 
      typedef const ForwardComposeMap<Graph, Map>& Type;
    };

    template <typename Graph, typename Map>
    struct Ref<BackwardComposeMap<Graph, Map> > { 
      typedef BackwardComposeMap<Graph, Map> Type; 
    };
    template <typename Graph, typename Map>
    struct Arg<BackwardComposeMap<Graph, Map> > { 
      typedef const BackwardComposeMap<Graph, Map>& Type; 
    };

    template <typename Map>
    struct Ref<dim2::XMap<Map> > { 
      typedef dim2::XMap<Map> Type;
    };
    template <typename Map>
    struct Arg<dim2::XMap<Map> > { 
      typedef const dim2::XMap<Map>& Type;
    };

    template <typename Map>
    struct Ref<dim2::YMap<Map> > { 
      typedef dim2::YMap<Map> Type;
    };
    template <typename Map>
    struct Arg<dim2::YMap<Map> > { 
      typedef const dim2::YMap<Map>& Type;
    };


    template <typename _Item>
    class MapReaderBase;
    
    template <typename _Item>
    class MapInverterBase : public MapReaderBase<_Item> {
    public:
      typedef _Item Item;
      virtual void read(std::istream&, const Item&) = 0;
      virtual Item read(std::istream&) const = 0;

      virtual MapInverterBase<_Item>* getInverter() {
	return this;
      }
    };

    template <typename _Item, typename _Map, typename _Reader>
    class MapReaderInverter : public MapInverterBase<_Item> {
    public:
      typedef _Item Item;
      typedef _Reader Reader;
      typedef typename Reader::Value Value;
      typedef _Map Map;
      typedef std::map<Value, Item, _reader_bits::Less<Value> > Inverse;

      typename _reader_bits::Ref<Map>::Type map;
      Reader reader;
      Inverse inverse;

      MapReaderInverter(typename _reader_bits::Arg<Map>::Type _map,
			const Reader& _reader) 
	: map(_map), reader(_reader) {}

      virtual ~MapReaderInverter() {}

      virtual void read(std::istream& is, const Item& item) {
	Value value;
	reader.read(is, value);
	map.set(item, value);
	typename Inverse::iterator it = inverse.find(value);
	if (it == inverse.end()) {
	  inverse.insert(std::make_pair(value, item));
	} else {
	  throw DataFormatError("Multiple label occurence");
	}
      }

      virtual Item read(std::istream& is) const {
	Value value;
	reader.read(is, value);	
	typename Inverse::const_iterator it = inverse.find(value);
	if (it != inverse.end()) {
	  return it->second;
	} else {
	  ErrorMessage msg;
	  msg << "Invalid label error"; 
	  throw DataFormatError(msg.message());
	}
      }      
    };

    template <typename _Item, typename _Reader>
    class SkipReaderInverter : public MapInverterBase<_Item> {
    public:
      typedef _Item Item;
      typedef _Reader Reader;
      typedef typename Reader::Value Value;
      typedef std::map<Value, Item, _reader_bits::Less<Value> > Inverse;

      Reader reader;

      SkipReaderInverter(const Reader& _reader) 
	: reader(_reader) {}

      virtual ~SkipReaderInverter() {}

      virtual void read(std::istream& is, const Item& item) {
	Value value;
	reader.read(is, value);
	typename Inverse::iterator it = inverse.find(value);
	if (it == inverse.end()) {
	  inverse.insert(std::make_pair(value, item));
	} else {
	  throw DataFormatError("Multiple label occurence error");
	}
      }

      virtual Item read(std::istream& is) const {
	Value value;
	reader.read(is, value);
	typename Inverse::const_iterator it = inverse.find(value);
	if (it != inverse.end()) {
	  return it->second;
	} else {
	  ErrorMessage msg;
	  msg << "Invalid label error: " << value; 
	  throw DataFormatError(msg.message());
	}
      }

    private:
      Inverse inverse;
    };

    template <typename _Item>    
    class MapReaderBase {
    public:
      typedef _Item Item;

      MapReaderBase() { _touched = false; }
      
      void touch(bool value = true) { _touched = value; }
      bool touched() const { return _touched; }

      virtual ~MapReaderBase() {}

      virtual void read(std::istream& is, const Item& item) = 0;
      virtual MapInverterBase<_Item>* getInverter() = 0;

    private:      
      bool _touched;

    };

    template <typename _Item, typename _Map, typename _Reader>
    class MapReader : public MapReaderBase<_Item> {
    public:
      typedef _Map Map;
      typedef _Reader Reader;
      typedef typename Reader::Value Value;
      typedef _Item Item;
      
      typename _reader_bits::Ref<Map>::Type map;
      Reader reader;

      MapReader(typename _reader_bits::Arg<Map>::Type _map, 
		const Reader& _reader) 
	: map(_map), reader(_reader) {}

      virtual ~MapReader() {}

      virtual void read(std::istream& is, const Item& item) {
	Value value;
	reader.read(is, value);
	map.set(item, value);
      }

      virtual MapInverterBase<_Item>* getInverter() {
	return new MapReaderInverter<Item, Map, Reader>(map, reader);
      }
    };


    template <typename _Item, typename _Reader>
    class SkipReader : public MapReaderBase<_Item> {
    public:
      typedef _Reader Reader;
      typedef typename Reader::Value Value;
      typedef _Item Item;

      Reader reader;
      SkipReader(const Reader& _reader) : reader(_reader) {}

      virtual ~SkipReader() {}

      virtual void read(std::istream& is, const Item&) {
	Value value;
	reader.read(is, value);
      }      

      virtual MapInverterBase<Item>* getInverter() {
	return new SkipReaderInverter<Item, Reader>(reader);
      }
    };

    template <typename _Item>
    class LabelReaderBase {
    public:
      typedef _Item Item;
      virtual ~LabelReaderBase() {}
      virtual Item read(std::istream& is) const = 0;
      virtual bool isLabelReader() const = 0;
      virtual LabelReaderBase<_Item>* clone() const = 0;
    };

    template <typename _Item, typename _BoxedLabelReader>
    class LabelReader : public LabelReaderBase<_Item> {
    public:
      typedef _Item Item;
      typedef _BoxedLabelReader BoxedLabelReader;
      
      const BoxedLabelReader& labelReader;

      LabelReader(const BoxedLabelReader& _labelReader) 
	: labelReader(_labelReader) {}

      virtual Item read(std::istream& is) const {
	Item item;
	labelReader.readLabel(is, item);
	return item;
      }

      virtual bool isLabelReader() const {
	return labelReader.isLabelReader();
      }
      
      LabelReader<Item, BoxedLabelReader>* clone() const {
	return new LabelReader<Item, BoxedLabelReader>(labelReader);
      }
    };

    template <typename _Item>
    class ItemStore {
    public:

      typedef _Item Item;

      ItemStore(Item& _item) : item(&_item) { 
	_touched = false; 
      }
      
      void touch() { _touched = true; }
      bool touched() const { return _touched; }

      void read(const Item& _item) {
	*item = _item;
      }
      
    private:
      Item* item;
      bool _touched;
    };

    class ValueReaderBase {
    public:
      virtual void read(std::istream&) {};
      ValueReaderBase() { _touched = false; }

      void touch() { _touched = true; }
      bool touched() const { return _touched; }

      virtual ~ValueReaderBase() {}
    private:
      bool _touched;
    };

    template <typename _Value, typename _Reader>
    class ValueReader : public ValueReaderBase {
    public:
      typedef _Value Value;
      typedef _Reader Reader;

      ValueReader(Value& _value, const Reader& _reader)
 	: value(_value), reader(_reader) {}

      virtual void read(std::istream& is) {
	reader.read(is, value);
      }
    private:
      Value& value;
      Reader reader;
    };

  }

  /// \ingroup lemon_io
  /// \brief Lemon Format reader class.
  /// 
  /// The Lemon Format contains several sections. We do not want to
  /// determine what sections are in a lemon file we give only a framework
  /// to read a section oriented format.
  ///
  /// In the Lemon Format each section starts with a line containing a
  /// \c \@ character on the first not white space position. This line
  /// is the header line of the section. Each of the next lines belong
  /// to this section until a line starting with \c \@ character is
  /// found. This line can start a new section or it can close the
  /// file with the \c \@end line.  The file format ignores the empty
  /// and comment lines. The line is comment line if it starts with a
  /// \c # character.
  ///
  /// The framework provides an abstract LemonReader::SectionReader class
  /// that defines the interface of a SectionReader. The SectionReader
  /// has the \c header() member function that gets a header line string and
  /// decides if it wants to process the next section. Several SectionReaders
  /// can be attached to a LemonReader and the first attached that can
  /// process the section will be used. Its \c read() member will be called
  /// with a stream containing the section. From this stream the empty and
  /// comment lines are filtered out.
  ///
  /// \relates GraphReader
  /// \relates NodeSetReader
  /// \relates EdgeSetReader
  /// \relates NodesReader
  /// \relates EdgesReader
  /// \relates AttributeReader
  class LemonReader {
  private:
    
    class FilterStreamBuf : public std::streambuf {
    public:

      typedef std::streambuf Parent;
      typedef Parent::char_type char_type;
      FilterStreamBuf(std::istream& is, int& num) 
	: _is(is), _base(0), _eptr(0), 
	  _num(num), skip_state(after_endl) {}

    protected:

      enum skip_state_type {
	no_skip,
	after_endl,
	comment_line
      };

      char_type small_buf[1];


      std::istream& _is;

      char_type* _base;
      char_type* _eptr;

      int& _num;

      skip_state_type skip_state;


      char_type* base() { return _base; }

      char_type* eptr() { return _eptr; }

      int_type blen() { return _eptr - _base; }

      void setb(char_type* buf, int_type len) {
	_base = buf;
	_eptr = buf + len;
      }
  
      virtual std::streambuf* setbuf(char *buf, std::streamsize len) {
	if (base()) return 0;
	if (buf != 0 && len >= int(sizeof(small_buf))) {
	  setb(buf, len);
	} else {
	  setb(small_buf, sizeof(small_buf));
	}
	setg(0, 0, 0);
	return this;
      }

      bool put_char(char c) {
	switch (skip_state) {
	case no_skip:
	  switch (c) {
	  case '\n': 
	    skip_state = after_endl;
	    return true;
	  default:
	    return true;
	  }
	case after_endl:
	  switch (c) {
	  case '@':
	    return false;
	  case '\n': 
	    return false;
	  case '#':
	    skip_state = comment_line;
	    return false;
	  default:
	    if (!isspace(c)) {
	      skip_state = no_skip;
	      return true;
	    } else {
	      return false;
	    }
	  }
	  break;
	case comment_line:
	  switch (c) {
	  case '\n': 
	    skip_state = after_endl;
	    return false;
	  default:
	    return false;
	  }
	}
	return false;
      }

      virtual int_type underflow() {
	char c;
	if ((c = _is.peek()) != EOF) {
	  if (c == '@') {
	    return EOF;
	  }
	} else {
	  return EOF;
	}
	char_type *ptr;
	for (ptr = base(); ptr != eptr(); ++ptr) {
	  if ((c = _is.get()) != EOF) {
	    if (c == '\n') ++_num;
	    if (put_char(c)) {
	      *ptr = c;
	    } else {
	      if (skip_state == after_endl && c == '@') {
		_is.putback(c);
		break;
	      }
	      --ptr;
	    }
	  } else {
	    break;
	  }
	}
	setg(base(), base(), ptr);
	return *base();
      }

      virtual int_type sync() {
	return EOF;
      }

    public:

      int line_num() const {
	int r = _num;
	for (char_type* p = gptr(); p != egptr(); ++p) {
	  if (*p == '\n') --r;
	}
	return r;
      }

    };

    static void skipPreSection(std::istream& is, int& line_num) {
      enum skip_state_type { skip, after_endl };

      skip_state_type skip_state = after_endl;
      char c;
      
      while ((c = is.get()) != EOF) {
	if (c == '\n') ++line_num;

	switch (skip_state) {
	case skip:
	  if (c == '\n') skip_state = after_endl;
	  break;
	case after_endl:
	  switch (c) {
	  case '@':
	    is.putback(c);
	    return;
	  case '\n':
	    continue;
	  default:
	    if (!isspace(c)) {
	      skip_state = skip;
	    }
	    break;
	  }
	}	
      }
    }

  public:

    /// \brief Abstract base class for reading a section.
    ///
    /// This class has an \c header() member function what get a 
    /// header line string and decides if it want to process the next 
    /// section. Several SectionReaders can be attached to an LemonReader 
    /// and the first attached what can process the section will be used. 
    /// Its \c read() member will called with a stream contains the section. 
    /// From this stream the empty lines and comments are filtered out.
    class SectionReader {
      friend class LemonReader;
    protected:
      /// \brief Constructor for SectionReader.
      ///
      /// Constructor for SectionReader. It attach this reader to
      /// the given LemonReader.
      SectionReader(LemonReader& reader) {
	reader.attach(*this);
      }

      virtual ~SectionReader() {}

      /// \brief Gives back true when the SectionReader can process 
      /// the section with the given header line.
      ///
      /// It gives back true when the SectionReader can process
      /// the section with the given header line.
      virtual bool header(const std::string& line) = 0;

      /// \brief Reader function of the section.
      ///
      /// It reads the content of the section.
      virtual void read(std::istream& is) = 0;

      /// \brief The given section missing in the file.
      ///
      /// The given section missing in the file.
      virtual void missing() {};
    };

    /// \brief Constructor for LemonReader.
    ///
    /// Constructor for LemonReader which reads from the given stream.
    LemonReader(std::istream& _is) 
      : is(&_is), own_is(false) {}

    /// \brief Constructor for LemonReader.
    ///
    /// Constructor for LemonReader which reads from the given file.
    LemonReader(const std::string& filename) 
      : is(0), own_is(true) {
      is = new std::ifstream(filename.c_str());
      if (is->fail()) {
	throw FileOpenError(filename);
      }
    }

    /// \brief Desctructor for LemonReader.
    ///
    /// Desctructor for LemonReader.
    ~LemonReader() {
      if (own_is) {
	delete is;
      }
    }

  private:
    LemonReader(const LemonReader&);
    void operator=(const LemonReader&);

    void attach(SectionReader& reader) {
      readers.push_back(std::make_pair(&reader, false));
    }

  public:
    /// \brief Executes the LemonReader.
    /// 
    /// It executes the LemonReader.
    void run() {
      int line_num = 0;
      std::string line;
      
      SectionReaders::iterator it;
      skipPreSection(*is, line_num);
      while ((++line_num, getline(*is, line)) && line.find("@end") != 0) {
	for (it = readers.begin(); it != readers.end(); ++it) {
	  if (it->first->header(line)) {
	    it->second = true;
	    char buf[2048];
	    FilterStreamBuf buffer(*is, line_num);
	    try {
	      buffer.pubsetbuf(buf, sizeof(buf));
	      std::istream ss(&buffer);
	      it->first->read(ss);
	      skipPreSection(*is, line_num);
	      break;
	    } catch (DataFormatError& error) {
	      error.line(buffer.line_num());
	      throw;
	    }	
	  }
	}
      }
      for (it = readers.begin(); it != readers.end(); ++it) {
	if (!it->second) {
	  try {
	    it->first->missing();
	  } catch (DataFormatError& error) {
	    error.line(line_num);
	    throw;
	  }	
	}
      }
    }


  private:

    std::istream* is;
    bool own_is;

    typedef std::vector<std::pair<SectionReader*, bool> > SectionReaders;
    SectionReaders readers;

  };

  /// \ingroup section_io
  /// \brief SectionReader for reading a graph's nodeset.
  ///
  /// The lemon format can store multiple graph nodesets with several
  /// maps.  The nodeset section's header line is \c \@nodeset \c
  /// nodeset_name, but the \c nodeset_name may be empty.
  ///
  /// The first line of the section contains the names of the maps separated
  /// with white spaces. Each next lines describes a node in the nodeset, and
  /// contains the mapped values for each map.
  ///
  /// If the nodeset contains an \c "label" named map then it will be regarded
  /// as id map. This map should contain only unique values and when the 
  /// \c readLabel() member will read a value from the given stream it will
  /// give back that node which is mapped to this value.
  ///
  /// \relates LemonReader
  template <typename _Graph, typename _Traits = DefaultReaderTraits>
  class NodeSetReader : public LemonReader::SectionReader {
    typedef LemonReader::SectionReader Parent;
  public:

    typedef _Graph Graph;
    typedef _Traits Traits;
    typedef typename Graph::Node Node;
    typedef typename Traits::Skipper DefaultSkipper;

    /// \brief Constructor.
    ///
    /// Constructor for NodeSetReader. It creates the NodeSetReader and
    /// attach it into the given LemonReader. The nodeset reader will
    /// add the read nodes to the given Graph. The reader will read
    /// the section when the \c section_name and the \c _name are the same. 
    NodeSetReader(LemonReader& _reader, 
		  Graph& _graph, 
		  const std::string& _name = std::string(),
		  const DefaultSkipper& _skipper = DefaultSkipper()) 
      : Parent(_reader), graph(_graph), name(_name), skipper(_skipper) {} 


    /// \brief Destructor.
    ///
    /// Destructor for NodeSetReader.
    virtual ~NodeSetReader() {
      for (typename MapReaders::iterator it = readers.begin(); 
	   it != readers.end(); ++it) {
	delete it->second;
      }
    }

  private:
    NodeSetReader(const NodeSetReader&);
    void operator=(const NodeSetReader&);
  
  public:

    /// \brief Add a new node map reader command for the reader.
    ///
    /// Add a new node map reader command for the reader.
    template <typename Map>
    NodeSetReader& readNodeMap(std::string label, Map& map) {
      return _readMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    template <typename Map>
    NodeSetReader& readNodeMap(std::string label, const Map& map) {
      return _readMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    /// \brief Add a new node map reader command for the reader.
    ///
    /// Add a new node map reader command for the reader.
    template <typename ItemReader, typename Map>
    NodeSetReader& readNodeMap(std::string label, Map& map, 
			       const ItemReader& ir = ItemReader()) {
      return _readMap<ItemReader, Map, typename _reader_bits::Arg<Map>::Type>
	(label, map, ir);
    }

    template <typename ItemReader, typename Map>
    NodeSetReader& readNodeMap(std::string label, const Map& map, 
			       const ItemReader& ir = ItemReader()) {
      return _readMap<ItemReader, Map, typename _reader_bits::Arg<Map>::Type>
	(label, map, ir);
    }

  private:

    template <typename ItemReader, typename Map, typename MapParameter>
    NodeSetReader& _readMap(std::string label, MapParameter map, 
			    const ItemReader& ir = ItemReader()) {
      checkConcept<concepts::WriteMap<Node, typename Map::Value>, Map>();
      checkConcept<_reader_bits::ItemReader<typename Map::Value>, ItemReader>();
      if (readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for node map: " << label;
	throw IoParameterError(msg.message());
      }      
      readers.insert(
        make_pair(label, new _reader_bits::
		  MapReader<Node, Map, ItemReader>(map, ir)));
      return *this;
    }

  public:

    /// \brief Add a new node map skipper command for the reader.
    ///
    /// Add a new node map skipper command for the reader.
    template <typename ItemReader>
    NodeSetReader& skipNodeMap(std::string label, 
                               const ItemReader& ir = ItemReader()) {
      if (readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for node map: " << label;
	throw IoParameterError(msg.message());
      }
      readers.insert(make_pair(label, new _reader_bits::
			       SkipReader<Node, ItemReader>(ir)));
      return *this;
    }

  protected:

    /// \brief Gives back true when the SectionReader can process 
    /// the section with the given header line.
    ///
    /// It gives back true when the header line starts with \c \@nodeset,
    /// and the header line's name and the nodeset's name are the same.
    virtual bool header(const std::string& line) {
      std::istringstream ls(line);
      std::string command;
      std::string id;
      ls >> command >> id;
      return command == "@nodeset" && name == id;
    }

    /// \brief Reader function of the section.
    ///
    /// It reads the content of the section.
    virtual void read(std::istream& is) {
      std::vector<_reader_bits::MapReaderBase<Node>* > index;
      std::string line;

      {
        getline(is, line);
        std::istringstream ls(line);
        std::string id;
        while (ls >> id) {
          typename MapReaders::iterator it = readers.find(id);
          if (it != readers.end()) {
            it->second->touch();
            index.push_back(it->second);
          } else {
            index.push_back(&skipper);
          }
          if (id == "label") {
            inverter.reset(index.back()->getInverter());
            index.back() = inverter.get();
          }
        }
      }
      for (typename MapReaders::iterator it = readers.begin();
	   it != readers.end(); ++it) {
	if (!it->second->touched()) {
	  ErrorMessage msg;
	  msg << "Map not found in file: " << it->first;
	  throw IoParameterError(msg.message());
	}
      }
      while (getline(is, line)) {	
	Node node = graph.addNode();
	std::istringstream ls(line);
	for (int i = 0; i < int(index.size()); ++i) {
	  index[i]->read(ls, node);
	}
      }
    }

    virtual void missing() {
      if (readers.empty()) return;
      ErrorMessage msg;
      msg << "NodeSet section not found in file: @nodeset " << name;
      throw IoParameterError(msg.message());
    }

  public:

    /// \brief Returns true if the nodeset can give back the node by its label.
    ///
    /// Returns true if the nodeset can give back the node by its label.
    /// It is possible only if an "label" named map was read.
    bool isLabelReader() const {
      return inverter.get() != 0;
    }

    /// \brief Gives back the node by its label.
    ///
    /// It reads an id from the stream and gives back which node belongs to
    /// it. It is possible only if there was read an "label" named map.
    void readLabel(std::istream& is, Node& node) const {
      node = inverter->read(is);
    } 

  private:

    typedef std::map<std::string, _reader_bits::MapReaderBase<Node>*> MapReaders;
    MapReaders readers;
   
    Graph& graph;   
    std::string name;
    _reader_bits::SkipReader<Node, DefaultSkipper> skipper;

    std::auto_ptr<_reader_bits::MapInverterBase<Node> > inverter;
  };

  /// \ingroup section_io
  /// \brief SectionReader for reading a bipartite graph's nodeset.
  ///
  /// The lemon format can store multiple bipartite graph nodesets
  /// with several maps. The bipartite graph nodeset section's header
  /// line is \c \@bpnodeset \c bpnodeset_name, but the \c bpnodeset_name
  /// may be empty.
  ///
  /// The first line of the section contains \c "&anodeset" and the
  /// the names of the A-node maps and regular maps separated with
  /// white spaces. Each next lines describes an A-node in the anodeset,
  /// and contains the mapped values for each map. If one of the line
  /// starts with \c "&bnodeset" then this line contains the names of
  /// the B-node maps and the regular node maps. And the remaining lines
  /// contains the mapped values to the B-nodes.
  ///
  /// If there is "label" named map then it should be defined in both
  /// nodeset, and it will be regarded as id map. This map should
  /// contain only unique values and when the \c readLabel() member
  /// will read a value from the given stream it will give back that
  /// node which is mapped to this value.
  ///
  /// \relates LemonReader
  template <typename _Graph, typename _Traits = DefaultReaderTraits>
  class BpNodeSetReader : public LemonReader::SectionReader {
    typedef LemonReader::SectionReader Parent;
  public:

    typedef _Graph Graph;
    typedef _Traits Traits;
    typedef typename Graph::Node Node;
    typedef typename Traits::Skipper DefaultSkipper;

    /// \brief Constructor.
    ///
    /// Constructor for BpNodeSetReader. It creates the BpNodeSetReader and
    /// attach it into the given LemonReader. The nodeset reader will
    /// add the read nodes to the given Graph. The reader will read
    /// the section when the \c section_name and the \c _name are the same. 
    BpNodeSetReader(LemonReader& _reader, 
		  Graph& _graph, 
		  const std::string& _name = std::string(),
		  const DefaultSkipper& _skipper = DefaultSkipper()) 
      : Parent(_reader), graph(_graph), name(_name), skipper(_skipper) {} 


    /// \brief Destructor.
    ///
    /// Destructor for BpNodeSetReader.
    virtual ~BpNodeSetReader() {
      for (typename MapReaders::iterator it = readers.begin(); 
	   it != readers.end(); ++it) {
	delete it->second;
      }
    }

  private:
    BpNodeSetReader(const BpNodeSetReader&);
    void operator=(const BpNodeSetReader&);
  
  public:

    /// \brief Add a new node map reader command for the reader.
    ///
    /// Add a new node map reader command for the reader.
    template <typename Map>
    BpNodeSetReader& readNodeMap(std::string label, Map& map) {
      return _readMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    template <typename Map>
    BpNodeSetReader& readNodeMap(std::string label, const Map& map) {
      return _readMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    /// \brief Add a new node map reader command for the reader.
    ///
    /// Add a new node map reader command for the reader.
    template <typename ItemReader, typename Map>
    BpNodeSetReader& readNodeMap(std::string label, Map& map, 
			       const ItemReader& ir = ItemReader()) {
      return _readMap<ItemReader, Map, typename _reader_bits::Arg<Map>::Type>
	(label, map, ir);
    }

    template <typename ItemReader, typename Map>
    BpNodeSetReader& readNodeMap(std::string label, const Map& map, 
			       const ItemReader& ir = ItemReader()) {
      return _readMap<ItemReader, Map, typename _reader_bits::Arg<Map>::Type>
	(label, map, ir);
    }

  private:

    template <typename ItemReader, typename Map, typename MapParameter>
    BpNodeSetReader& _readMap(std::string label, MapParameter map, 
			    const ItemReader& ir = ItemReader()) {
      checkConcept<concepts::WriteMap<Node, typename Map::Value>, Map>();
      checkConcept<_reader_bits::ItemReader<typename Map::Value>, ItemReader>();
      if (areaders.find(label) != areaders.end() ||
	  breaders.find(label) != breaders.end() ||
	  readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for node map: " << label;
	throw IoParameterError(msg.message());
      }      
      readers.insert(make_pair(label, new _reader_bits::
		  MapReader<Node, Map, ItemReader>(map, ir)));
      return *this;
    }

  public:

    /// \brief Add a new A-node map reader command for the reader.
    ///
    /// Add a new A-node map reader command for the reader.
    template <typename Map>
    BpNodeSetReader& readANodeMap(std::string label, Map& map) {
      return _readAMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    template <typename Map>
    BpNodeSetReader& readANodeMap(std::string label, const Map& map) {
      return _readAMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    /// \brief Add a new A-node map reader command for the reader.
    ///
    /// Add a new A-node map reader command for the reader.
    template <typename ItemReader, typename Map>
    BpNodeSetReader& readANodeMap(std::string label, Map& map, 
			       const ItemReader& ir = ItemReader()) {
      return _readAMap<ItemReader, Map, typename _reader_bits::Arg<Map>::Type>
	(label, map, ir);
    }

    template <typename ItemReader, typename Map>
    BpNodeSetReader& readANodeMap(std::string label, const Map& map, 
			       const ItemReader& ir = ItemReader()) {
      return _readAMap<ItemReader, Map, typename _reader_bits::Arg<Map>::Type>
	(label, map, ir);
    }

  private:

    template <typename ItemReader, typename Map, typename MapParameter>
    BpNodeSetReader& _readAMap(std::string label, MapParameter map, 
			    const ItemReader& ir = ItemReader()) {
      checkConcept<concepts::WriteMap<Node, typename Map::Value>, Map>();
      checkConcept<_reader_bits::ItemReader<typename Map::Value>, ItemReader>();
      if (label == "label") {
	throw IoParameterError("Label cannot be A-node map");
      }
      if (areaders.find(label) != areaders.end() ||
	  readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for A-node map: " << label;
	throw IoParameterError(msg.message());
      }
      areaders.insert(make_pair(label, new _reader_bits::
				MapReader<Node, Map, ItemReader>(map, ir)));
      return *this;
    }

  public:

    /// \brief Add a new B-node map reader command for the reader.
    ///
    /// Add a new B-node map reader command for the reader.
    template <typename Map>
    BpNodeSetReader& readBNodeMap(std::string label, Map& map) {
      return _readBMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    template <typename Map>
    BpNodeSetReader& readBNodeMap(std::string label, const Map& map) {
      return _readBMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    /// \brief Add a new B-node map reader command for the reader.
    ///
    /// Add a new B-node map reader command for the reader.
    template <typename ItemReader, typename Map>
    BpNodeSetReader& readBNodeMap(std::string label, Map& map, 
			       const ItemReader& ir = ItemReader()) {
      return _readBMap<ItemReader, Map, typename _reader_bits::Arg<Map>::Type>
	(label, map, ir);
    }

    template <typename ItemReader, typename Map>
    BpNodeSetReader& readBNodeMap(std::string label, const Map& map, 
			       const ItemReader& ir = ItemReader()) {
      return _readBMap<ItemReader, Map, typename _reader_bits::Arg<Map>::Type>
	(label, map, ir);
    }

  private:

    template <typename ItemReader, typename Map, typename MapParameter>
    BpNodeSetReader& _readBMap(std::string label, MapParameter map, 
			    const ItemReader& ir = ItemReader()) {
      checkConcept<concepts::WriteMap<Node, typename Map::Value>, Map>();
      checkConcept<_reader_bits::ItemReader<typename Map::Value>, ItemReader>();
      if (label == "label") {
	throw IoParameterError("Label cannot be B-node map");
      }
      if (breaders.find(label) != breaders.end() ||
	  readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for B-node map: " << label;
	throw IoParameterError(msg.message());
      }
      breaders.insert(make_pair(label, new _reader_bits::
				MapReader<Node, Map, ItemReader>(map, ir)));
      return *this;
    }

  public:

    /// \brief Add a new node map skipper command for the reader.
    ///
    /// Add a new node map skipper command for the reader.
    template <typename ItemReader>
    BpNodeSetReader& skipNodeMap(std::string label, 
				 const ItemReader& ir = ItemReader()) {
      if (areaders.find(label) != areaders.end() ||
	  breaders.find(label) != breaders.end() ||
	  readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for node map: " << label;
	throw IoParameterError(msg.message());
      }
      readers.insert(make_pair(label, new _reader_bits::
			       SkipReader<Node, ItemReader>(ir)));
      return *this;
    }

    /// \brief Add a new A-node map skipper command for the reader.
    ///
    /// Add a new A-node map skipper command for the reader.
    template <typename ItemReader>
    BpNodeSetReader& skipANodeMap(std::string label, 
				  const ItemReader& ir = ItemReader()) {
      if (label == "label") {
	throw IoParameterError("Label cannot be A-node map");
      }
      if (areaders.find(label) != areaders.end() ||
	  readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for A-node map: " << label;
	throw IoParameterError(msg.message());
      }
      areaders.insert(make_pair(label, new _reader_bits::
				SkipReader<Node, ItemReader>(ir)));
      return *this;
    }

    /// \brief Add a new B-node map skipper command for the reader.
    ///
    /// Add a new B-node map skipper command for the reader.
    template <typename ItemReader>
    BpNodeSetReader& skipBNodeMap(std::string label, 
				  const ItemReader& ir = ItemReader()) {
      if (label == "label") {
	throw IoParameterError("Label cannot be B-node map");
      }
      if (breaders.find(label) != breaders.end() ||
	  readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for B-node map: " << label;
	throw IoParameterError(msg.message());
      }
      breaders.insert(make_pair(label, new _reader_bits::
				SkipReader<Node, ItemReader>(ir)));
      return *this;
    }


  protected:

    /// \brief Gives back true when the SectionReader can process 
    /// the section with the given header line.
    ///
    /// It gives back true when the header line starts with \c \@nodeset,
    /// and the header line's name and the nodeset's name are the same.
    virtual bool header(const std::string& line) {
      std::istringstream ls(line);
      std::string command;
      std::string id;
      ls >> command >> id;
      return command == "@bpnodeset" && name == id;
    }

    /// \brief Reader function of the section.
    ///
    /// It reads the content of the section.
    virtual void read(std::istream& is) {
      std::string line;
      {
	std::vector<_reader_bits::MapReaderBase<Node>* > index;
	{
	  getline(is, line);
	  std::istringstream ls(line);
	  std::string id;
	  ls >> id;
	  if (id != "&anodeset") {
	    throw IoParameterError("Cannot find &anodeset subsection");
	  }
	  while (ls >> id) {
	    typename MapReaders::iterator it = readers.find(id);
	    typename MapReaders::iterator ait = areaders.find(id);
	    if (it != readers.end()) {
	      it->second->touch();
	      index.push_back(it->second);
	    } else if (ait != areaders.end()) {
	      ait->second->touch();
	      index.push_back(ait->second);
	    } else {
	      index.push_back(&skipper);
	    }
	    if (id == "label") {
	      inverter.reset(index.back()->getInverter());
	      index.back() = inverter.get();
	    }
	  }
	}
	for (typename MapReaders::iterator it = areaders.begin();
	     it != areaders.end(); ++it) {
	  if (!it->second->touched()) {
	    ErrorMessage msg;
	    msg << "Map not found in file: " << it->first;
	    throw IoParameterError(msg.message());
	  }
	}
	for (typename MapReaders::iterator it = readers.begin();
	     it != readers.end(); ++it) {
	  if (!it->second->touched()) {
	    ErrorMessage msg;
	    msg << "Map not found in file: " << it->first;
	    throw IoParameterError(msg.message());
	  }
	  it->second->touch(false);
	}

	while (getline(is, line)) {
	  if (line[0] == '&') {
	    std::istringstream ls(line);
	    std::string id;
	    ls >> id;
	    if (id == "&bnodeset") break;
	  }
	  Node node = graph.addANode();
	  std::istringstream ls(line);
	  for (int i = 0; i < int(index.size()); ++i) {
	    index[i]->read(ls, node);
	  }
	}
      }

      {
	std::vector<_reader_bits::MapReaderBase<Node>* > index;
	{
	  std::istringstream ls(line);
	  std::string id;
	  ls >> id;
	  if (id != "&bnodeset") {
	    throw IoParameterError("Cannot find &bnodeset subsection");
	  }
	  while (ls >> id) {
	    typename MapReaders::iterator it = readers.find(id);
	    typename MapReaders::iterator bit = breaders.find(id);
	    if (it != readers.end()) {
	      it->second->touch();
	      index.push_back(it->second);
	    } else if (bit != breaders.end()) {
	      bit->second->touch();
	      index.push_back(bit->second);
	    } else {
	      index.push_back(&skipper);
	    }
	    if (id == "label" && inverter.get() != 0) {
	      index.back() = inverter.get();
	    }
	  }
	}
	for (typename MapReaders::iterator it = breaders.begin();
	     it != breaders.end(); ++it) {
	  if (!it->second->touched()) {
	    ErrorMessage msg;
	    msg << "Map not found in file: " << it->first;
	    throw IoParameterError(msg.message());
	  }
	}
	for (typename MapReaders::iterator it = readers.begin();
	     it != readers.end(); ++it) {
	  if (!it->second->touched()) {
	    ErrorMessage msg;
	    msg << "Map not found in file: " << it->first;
	    throw IoParameterError(msg.message());
	  }
	}
	while (getline(is, line)) {	
	  Node node = graph.addBNode();
	  std::istringstream ls(line);
	  for (int i = 0; i < int(index.size()); ++i) {
	    index[i]->read(ls, node);
	  }
	}
      }
    }

    virtual void missing() {
      if (readers.empty()) return;
      ErrorMessage msg;
      msg << "BpNodeSet section not found in file: @bpnodeset " << name;
      throw IoParameterError(msg.message());
    }

  public:

    /// \brief Returns true if the nodeset can give back the node by its label.
    ///
    /// Returns true if the nodeset can give back the node by its label.
    /// It is possible only if an "label" named map was read.
    bool isLabelReader() const {
      return inverter.get() != 0;
    }

    /// \brief Gives back the node by its label.
    ///
    /// It reads an id from the stream and gives back which node belongs to
    /// it. It is possible only if there was read an "label" named map.
    void readLabel(std::istream& is, Node& node) const {
      node = inverter->read(is);
    } 

  private:

    typedef std::map<std::string, _reader_bits::MapReaderBase<Node>*> 
    MapReaders;
    
    MapReaders areaders, breaders, readers;
   
    Graph& graph;
    std::string name;
    _reader_bits::SkipReader<Node, DefaultSkipper> skipper;

    std::auto_ptr<_reader_bits::MapInverterBase<Node> > inverter;
  };


  /// \ingroup section_io
  /// \brief SectionReader for reading a graph's edgeset.
  ///
  /// The lemon format can store multiple graph edgesets with several maps.
  /// The edgeset section's header line is \c \@edgeset \c edgeset_name, but the
  /// \c edgeset_name may be empty.
  ///
  /// The first line of the section contains the names of the maps separated
  /// with white spaces. Each next lines describes an edge in the edgeset. The
  /// line contains the source and the target nodes' id and the mapped 
  /// values for each map.
  ///
  /// If the edgeset contains an \c "label" named map then it will be regarded
  /// as id map. This map should contain only unique values and when the 
  /// \c readLabel() member will read a value from the given stream it will
  /// give back that edge which is mapped to this value.
  ///
  /// The edgeset reader needs a node id reader to identify which nodes
  /// have to be connected. If a NodeSetReader reads an "label" named map,
  /// it will be able to resolve the nodes by ids.
  ///
  /// \relates LemonReader
  template <typename _Graph, typename _Traits = DefaultReaderTraits>
  class EdgeSetReader : public LemonReader::SectionReader {
    typedef LemonReader::SectionReader Parent;
  public:

    typedef _Graph Graph;
    typedef _Traits Traits;
    typedef typename Graph::Node Node;
    typedef typename Graph::Edge Edge;
    typedef typename Traits::Skipper DefaultSkipper;

    /// \brief Constructor.
    ///
    /// Constructor for EdgeSetReader. It creates the EdgeSetReader and
    /// attach it into the given LemonReader. The edgeset reader will
    /// add the read edges to the given Graph. It will use the given
    /// node id reader to read the source and target nodes of the edges.
    /// The reader will read the section only if the \c _name and the 
    /// \c edgset_name are the same. 
    template <typename NodeLabelReader>
    EdgeSetReader(LemonReader& _reader, 
		  Graph& _graph, 
		  const NodeLabelReader& _nodeLabelReader, 
		  const std::string& _name = std::string(),
		  const DefaultSkipper& _skipper = DefaultSkipper()) 
      : Parent(_reader), graph(_graph), name(_name), skipper(_skipper) {
      checkConcept<_reader_bits::ItemLabelReader<Node>, NodeLabelReader>();
      nodeLabelReader.reset(new _reader_bits::
			 LabelReader<Node, NodeLabelReader>(_nodeLabelReader));
    }
    /// \brief Destructor.
    ///
    /// Destructor for EdgeSetReader.
    virtual ~EdgeSetReader() {
      for (typename MapReaders::iterator it = readers.begin(); 
	   it != readers.end(); ++it) {
	delete it->second;
      }
    }

  private:
    EdgeSetReader(const EdgeSetReader&);
    void operator=(const EdgeSetReader&);

  public:

    /// \brief Add a new edge map reader command for the reader.
    ///
    /// Add a new edge map reader command for the reader.
    template <typename Map>
    EdgeSetReader& readEdgeMap(std::string label, Map& map) {
      return _readMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    template <typename Map>
    EdgeSetReader& readEdgeMap(std::string label, const Map& map) {
      return _readMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    /// \brief Add a new edge map reader command for the reader.
    ///
    /// Add a new edge map reader command for the reader.
    template <typename ItemReader, typename Map>
    EdgeSetReader& readEdgeMap(std::string label, Map& map, 
                               const ItemReader& ir = ItemReader()) {
      return _readMap<ItemReader, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map, ir);
    }

    template <typename ItemReader, typename Map>
    EdgeSetReader& readEdgeMap(std::string label, const Map& map, 
			       const ItemReader& ir = ItemReader()) {
      return _readMap<ItemReader, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map, ir);
    }

  private:

    template <typename ItemReader, typename Map, typename MapParameter>
    EdgeSetReader& _readMap(std::string label, MapParameter map, 
			    const ItemReader& ir = ItemReader()) {
      checkConcept<concepts::WriteMap<Edge, typename Map::Value>, Map>();
      checkConcept<_reader_bits::ItemReader<typename Map::Value>, ItemReader>();
      if (readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for edge map: " << label;
	throw IoParameterError(msg.message());
      }
      readers.insert(
	make_pair(label, new _reader_bits::
		  MapReader<Edge, Map, ItemReader>(map, ir)));
      return *this;
    }

  public:

    /// \brief Add a new edge map skipper command for the reader.
    ///
    /// Add a new edge map skipper command for the reader.
    template <typename ItemReader>
    EdgeSetReader& skipEdgeMap(std::string label, 
			       const ItemReader& ir = ItemReader()) {
      if (readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for edge map: " << label;
	throw IoParameterError(msg.message());
      }
      readers.insert(make_pair(label, new _reader_bits::
			       SkipReader<Edge, ItemReader>(ir)));
      return *this;
    }

  protected:

    /// \brief Gives back true when the SectionReader can process 
    /// the section with the given header line.
    ///
    /// It gives back true when the header line starts with \c \@edgeset,
    /// and the header line's name and the edgeset's name are the same.
    /// The sections with \@uedgeset head line could be read with this
    /// section reader too.
    virtual bool header(const std::string& line) {
      std::istringstream ls(line);
      std::string command;
      std::string id;
      ls >> command >> id;
      return (command == "@edgeset" || command == "@uedgeset") && name == id;
    }

    /// \brief Reader function of the section.
    ///
    /// It reads the content of the section.
    virtual void read(std::istream& is) {
      if (!nodeLabelReader->isLabelReader()) {
	throw DataFormatError("Cannot find nodeset or label map");
      }
      std::vector<_reader_bits::MapReaderBase<Edge>* > index;
      std::string line;

      {
        getline(is, line);
        std::istringstream ls(line);	
        std::string id;
        while (ls >> id) {
          typename MapReaders::iterator it = readers.find(id);
          if (it != readers.end()) {
            index.push_back(it->second);
            it->second->touch();
          } else {
            index.push_back(&skipper);
          }
          if (id == "label") {
            inverter.reset(index.back()->getInverter());
            index.back() = inverter.get();
          }
        }
      }
      for (typename MapReaders::iterator it = readers.begin();
	   it != readers.end(); ++it) {
	if (!it->second->touched()) {
	  ErrorMessage msg;
	  msg << "Map not found in file: " << it->first;
	  throw IoParameterError(msg.message());
	}
      }
      while (getline(is, line)) {	
	std::istringstream ls(line);
	Node from = nodeLabelReader->read(ls);
	Node to = nodeLabelReader->read(ls);
	Edge edge = graph.addEdge(from, to);
	for (int i = 0; i < int(index.size()); ++i) {
	  index[i]->read(ls, edge);
	}
      }
    }

    virtual void missing() {
      if (readers.empty()) return;
      ErrorMessage msg;
      msg << "EdgeSet section not found in file: @edgeset " << name;
      throw IoParameterError(msg.message());
    }

  public:

    /// \brief Returns true if the edgeset can give back the edge by its label.
    ///
    /// Returns true if the edgeset can give back the edge by its label.
    /// It is possible only if an "label" named map was read.
    bool isLabelReader() const {
      return inverter.get() != 0;
    }

    /// \brief Gives back the edge by its label.
    ///
    /// It reads an id from the stream and gives back which edge belongs to
    /// it. It is possible only if there was read an "label" named map.
    void readLabel(std::istream& is, Edge& edge) const {
      edge = inverter->read(is);
    } 

  private:

    typedef std::map<std::string, _reader_bits::MapReaderBase<Edge>*> 
    MapReaders;
    
    MapReaders readers;
   
    Graph& graph;   
    std::string name;
    _reader_bits::SkipReader<Edge, DefaultSkipper> skipper;

    std::auto_ptr<_reader_bits::MapInverterBase<Edge> > inverter;
    std::auto_ptr<_reader_bits::LabelReaderBase<Node> > nodeLabelReader;
  };

  /// \ingroup section_io
  /// \brief SectionReader for reading a undirected graph's edgeset.
  ///
  /// The lemon format can store multiple undirected edgesets with several 
  /// maps. The undirected edgeset section's header line is \c \@uedgeset 
  /// \c uedgeset_name, but the \c uedgeset_name may be empty.
  ///
  /// The first line of the section contains the names of the maps separated
  /// with white spaces. Each next lines describes an edge in the edgeset. The
  /// line contains the connected nodes' id and the mapped values for each map.
  ///
  /// The section can handle the directed as a syntactical sugar. Two
  /// undirected edge map describes one directed edge map. This two maps
  /// are the forward map and the backward map and the names of this map
  /// is near the same just with a prefix \c '+' or \c '-' character 
  /// difference.
  ///
  /// If the edgeset contains an \c "label" named map then it will be regarded
  /// as id map. This map should contain only unique values and when the 
  /// \c readLabel() member will read a value from the given stream it will
  /// give back that uicted edge which is mapped to this value.
  ///
  /// The undirected edgeset reader needs a node id reader to identify which 
  /// nodes have to be connected. If a NodeSetReader reads an "label" named 
  /// map, it will be able to resolve the nodes by ids.
  ///
  /// \relates LemonReader
  template <typename _Graph, typename _Traits = DefaultReaderTraits>
  class UEdgeSetReader : public LemonReader::SectionReader {
    typedef LemonReader::SectionReader Parent;
  public:

    typedef _Graph Graph;
    typedef _Traits Traits;
    typedef typename Graph::Node Node;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::UEdge UEdge;
    typedef typename Traits::Skipper DefaultSkipper;

    /// \brief Constructor.
    ///
    /// Constructor for UEdgeSetReader. It creates the UEdgeSetReader 
    /// and attach it into the given LemonReader. The undirected edgeset 
    /// reader will add the read undirected edges to the given Graph. It 
    /// will use the given node id reader to read the source and target 
    /// nodes of the edges. The reader will read the section only if the 
    /// \c _name and the \c uedgset_name are the same. 
    template <typename NodeLabelReader>
    UEdgeSetReader(LemonReader& _reader, 
		       Graph& _graph, 
		       const NodeLabelReader& _nodeLabelReader, 
		       const std::string& _name = std::string(),
		       const DefaultSkipper& _skipper = DefaultSkipper()) 
      : Parent(_reader), graph(_graph), name(_name), skipper(_skipper) {
      checkConcept<_reader_bits::ItemLabelReader<Node>, NodeLabelReader>();
      nodeLabelReader.reset(new _reader_bits::
			 LabelReader<Node, NodeLabelReader>(_nodeLabelReader));
    }
    /// \brief Destructor.
    ///
    /// Destructor for UEdgeSetReader.
    virtual ~UEdgeSetReader() {
      for (typename MapReaders::iterator it = readers.begin(); 
	   it != readers.end(); ++it) {
	delete it->second;
      }
    }

  private:
    UEdgeSetReader(const UEdgeSetReader&);
    void operator=(const UEdgeSetReader&);

  public:

    /// \brief Add a new undirected edge map reader command for the reader.
    ///
    /// Add a new edge undirected map reader command for the reader.
    template <typename Map>
    UEdgeSetReader& readUEdgeMap(std::string label, Map& map) {
      return _readMap<
	typename Traits::template Reader<typename Map::Value>, Map, 
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    template <typename Map>
    UEdgeSetReader& readUEdgeMap(std::string label, const Map& map) {
      return _readMap<
	typename Traits::template Reader<typename Map::Value>, Map, 
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    /// \brief Add a new undirected edge map reader command for the reader.
    ///
    /// Add a new edge undirected map reader command for the reader.
    template <typename ItemReader, typename Map>
    UEdgeSetReader& readUEdgeMap(std::string label, Map& map, 
                                 const ItemReader& ir = ItemReader()) {
      return _readMap<ItemReader, Map, typename _reader_bits::Arg<Map>::Type>
	(label, map, ir);
    }

    template <typename ItemReader, typename Map>
    UEdgeSetReader& readUEdgeMap(std::string label, const Map& map, 
                                 const ItemReader& ir = ItemReader()) {
      return _readMap<ItemReader, Map, typename _reader_bits::Arg<Map>::Type >
	(label, map, ir);
    }

  private:

    template <typename ItemReader, typename Map, typename MapParameter>
    UEdgeSetReader& _readMap(std::string label, MapParameter map,
                             const ItemReader& ir = ItemReader()) {
      checkConcept<concepts::WriteMap<UEdge, typename Map::Value>, Map>();
      checkConcept<_reader_bits::ItemReader<typename Map::Value>, ItemReader>();
      if (readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for edge map: " << label;
	throw IoParameterError(msg.message());
      }
      readers.insert(
	make_pair(label, new _reader_bits::
		  MapReader<UEdge, Map, ItemReader>(map, ir)));
      return *this;
    }

  public:

    /// \brief Add a new undirected edge map skipper command for the reader.
    ///
    /// Add a new undirected edge map skipper command for the reader.
    template <typename ItemReader>
    UEdgeSetReader& skipUEdgeMap(std::string label, 
                                 const ItemReader& ir = ItemReader()) {
      if (readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for node map: " << label;
	throw IoParameterError(msg.message());
      }
      readers.insert(make_pair(label, new _reader_bits::
			       SkipReader<UEdge, ItemReader>(ir)));
      return *this;
    }

    /// \brief Add a new directed edge map reader command for the reader.
    ///
    /// Add a new directed edge map reader command for the reader.
    template <typename Map>
    UEdgeSetReader& readEdgeMap(std::string label, Map& map) {
      return _readDirMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    template <typename Map>
    UEdgeSetReader& readEdgeMap(std::string label, const Map& map) {
      return _readDirMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    /// \brief Add a new directed edge map reader command for the reader.
    ///
    /// Add a new directed edge map reader command for the reader.
    template <typename ItemReader, typename Map>
    UEdgeSetReader& readEdgeMap(std::string label, Map& map, 
				    const ItemReader& ir = ItemReader()) {
      return _readDirMap<ItemReader, Map, 
        typename _reader_bits::Arg<Map>::Type>(label, map, ir);
    }

    template <typename ItemReader, typename Map>
    UEdgeSetReader& readEdgeMap(std::string label, const Map& map, 
				    const ItemReader& ir = ItemReader()) {
      return _readDirMap<ItemReader, Map, 
        typename _reader_bits::Arg<Map>::Type>(label, map, ir);
    }

  private:

    template <typename ItemReader, typename Map, typename MapParameter>
    UEdgeSetReader& _readDirMap(std::string label, MapParameter map,
				    const ItemReader& ir = ItemReader()) { 
      checkConcept<_reader_bits::ItemReader<typename Map::Value>, ItemReader>();
      checkConcept<concepts::WriteMap<Edge, typename Map::Value>, Map>();
      readUEdgeMap("+" + label, 
                   _reader_bits::forwardComposeMap(graph, map), ir);
      readUEdgeMap("-" + label, 
                   _reader_bits::backwardComposeMap(graph, map), ir);
      return *this;      
    }

  public:

    /// \brief Add a new directed edge map skipper command for the reader.
    ///
    /// Add a new directed edge map skipper command for the reader.
    template <typename ItemReader>
    UEdgeSetReader& skipEdgeMap(std::string label, 
                                const ItemReader& ir = ItemReader()) {
      skipUEdgeMap("+" + label, ir);
      skipUEdgeMap("-" + label, ir);
      return *this;
    }

  protected:

    /// \brief Gives back true when the SectionReader can process 
    /// the section with the given header line.
    ///
    /// It gives back true when the header line starts with \c \@uedgeset,
    /// and the header line's name and the edgeset's name are the same.
    /// The sections with \@edgeset head line could be read with this
    /// section reader too.
    virtual bool header(const std::string& line) {
      std::istringstream ls(line);
      std::string command;
      std::string id;
      ls >> command >> id;
      return (command == "@edgeset" || command == "@uedgeset") && name == id;
    }

    /// \brief Reader function of the section.
    ///
    /// It reads the content of the section.
    virtual void read(std::istream& is) {
      if (!nodeLabelReader->isLabelReader()) {
	throw DataFormatError("Cannot find nodeset or label map");
      }
      std::vector<_reader_bits::MapReaderBase<UEdge>* > index;
      std::string line;

      {
        getline(is, line);
        std::istringstream ls(line);	
        std::string id;
        while (ls >> id) {
          typename MapReaders::iterator it = readers.find(id);
          if (it != readers.end()) {
            index.push_back(it->second);
            it->second->touch();
          } else {
            index.push_back(&skipper);
          }
          if (id == "label") {
            inverter.reset(index.back()->getInverter());
            index.back() = inverter.get();
          }
        }
        for (typename MapReaders::iterator it = readers.begin();
             it != readers.end(); ++it) {
          if (!it->second->touched()) {
            ErrorMessage msg;
            msg << "Map not found in file: " << it->first;
            throw IoParameterError(msg.message());
          }
        }
      }
      while (getline(is, line)) {	
	std::istringstream ls(line);
	Node from = nodeLabelReader->read(ls);
	Node to = nodeLabelReader->read(ls);
	UEdge edge = graph.addEdge(from, to);
	for (int i = 0; i < int(index.size()); ++i) {
	  index[i]->read(ls, edge);
	}
      }
    }

    virtual void missing() {
      if (readers.empty()) return;
      ErrorMessage msg;
      msg << "UEdgeSet section not found in file: @uedgeset " << name;
      throw IoParameterError(msg.message());
    }

  public:

    /// \brief Returns true if the edgeset can give back the edge by its label.
    ///
    /// Returns true if the edgeset can give back the undirected edge by its 
    /// id. It is possible only if an "label" named map was read.
    bool isLabelReader() const {
      return inverter.get() != 0;
    }

    /// \brief Gives back the undirected edge by its label.
    ///
    /// It reads an id from the stream and gives back which undirected edge 
    /// belongs to it. It is possible only if there was read an "label" named map.
    void readLabel(std::istream& is, UEdge& uedge) const {
      uedge = inverter->read(is);
    } 

    /// \brief Gives back the directed edge by its label.
    ///
    /// It reads an id from the stream and gives back which directed edge 
    /// belongs to it. The directed edge id is the \c '+' or \c '-' character
    /// and the undirected edge id. It is possible only if there was read 
    /// an "label" named map.
    void readLabel(std::istream& is, Edge& edge) const {
      char c;
      is >> c;
      UEdge uedge = inverter->read(is);
      if (c == '+') {
	edge = graph.direct(uedge, true);
      } else if (c == '-') {
        edge = graph.direct(uedge, false);
      } else {
	throw DataFormatError("Wrong id format for edge "
			      "in undirected edgeset");
      }
    } 

  private:

    typedef std::map<std::string, 
		     _reader_bits::MapReaderBase<UEdge>*> MapReaders;
    MapReaders readers;
   
    Graph& graph;   
    std::string name;
    _reader_bits::SkipReader<UEdge, DefaultSkipper> skipper;

    std::auto_ptr<_reader_bits::MapInverterBase<UEdge> > inverter;
    std::auto_ptr<_reader_bits::LabelReaderBase<Node> > nodeLabelReader;
  };

  /// \ingroup section_io
  /// \brief SectionReader for reading labeled nodes.
  ///
  /// The nodes section's header line is \c \@nodes \c nodes_name, but the
  /// \c nodes_name may be empty.
  ///
  /// Each line in the section contains the name of the node 
  /// and then the node id. 
  ///
  /// \relates LemonReader
  template <typename _Graph>
  class NodeReader : public LemonReader::SectionReader {
    typedef LemonReader::SectionReader Parent;
    typedef _Graph Graph;
    typedef typename Graph::Node Node;
  public:
    
    /// \brief Constructor.
    ///
    /// Constructor for NodeReader. It creates the NodeReader and
    /// attach it into the given LemonReader. It will use the given
    /// node id reader to give back the nodes. The reader will read the 
    /// section only if the \c _name and the \c nodes_name are the same. 
    template <typename _LabelReader>
    NodeReader(LemonReader& _reader, const _LabelReader& _labelReader, 
	       const std::string& _name = std::string()) 
      : Parent(_reader), name(_name) {
      checkConcept<_reader_bits::ItemLabelReader<Node>, _LabelReader>();
      nodeLabelReader.reset(new _reader_bits::
			 LabelReader<Node, _LabelReader>(_labelReader));
    }

    /// \brief Destructor.
    ///
    /// Destructor for NodeReader.
    virtual ~NodeReader() {}

  private:
    NodeReader(const NodeReader&);
    void operator=(const NodeReader&);

  public:

    /// \brief Add a node reader command for the NodeReader.
    ///
    /// Add a node reader command for the NodeReader.
    void readNode(std::string label, Node& item) {
      if (readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for node: " << label;
	throw IoParameterError(msg.message());
      }
      readers.insert(make_pair(label, _reader_bits::ItemStore<Node>(item)));
    }

  protected:

    /// \brief Gives back true when the SectionReader can process 
    /// the section with the given header line.
    ///
    /// It gives back true when the header line start with \c \@nodes,
    /// and the header line's name and the reader's name are the same.
    virtual bool header(const std::string& line) {
      std::istringstream ls(line);
      std::string command;
      std::string id;
      ls >> command >> id;
      return command == "@nodes" && name == id;
    }

    /// \brief Reader function of the section.
    ///
    /// It reads the content of the section.
    virtual void read(std::istream& is) {
      if (!nodeLabelReader->isLabelReader()) {
	throw DataFormatError("Cannot find nodeset or label map");
      }
      std::string line;
      while (getline(is, line)) {
	std::istringstream ls(line);
	std::string id;
	ls >> id;
	typename NodeReaders::iterator it = readers.find(id);
	if (it != readers.end()) {
	  it->second.read(nodeLabelReader->read(ls));
	  it->second.touch();
	}	
      }
      for (typename NodeReaders::iterator it = readers.begin();
	   it != readers.end(); ++it) {
	if (!it->second.touched()) {
	  ErrorMessage msg;
	  msg << "Node not found in file: " << it->first;
	  throw IoParameterError(msg.message());
	}
      }
    }

    virtual void missing() {
      if (readers.empty()) return;
      ErrorMessage msg;
      msg << "Nodes section not found in file: @nodes " << name;
      throw IoParameterError(msg.message());
    }
    
  private:

    std::string name;

    typedef std::map<std::string, _reader_bits::ItemStore<Node> > NodeReaders;
    NodeReaders readers;
    std::auto_ptr<_reader_bits::LabelReaderBase<Node> > nodeLabelReader;
  };

  /// \ingroup section_io
  /// \brief SectionReader for reading labeled edges.
  ///
  /// The edges section's header line is \c \@edges \c edges_name, but the
  /// \c edges_name may be empty.
  ///
  /// Each line in the section contains the name of the edge 
  /// and then the edge id. 
  ///
  /// \relates LemonReader
  template <typename _Graph>
  class EdgeReader : public LemonReader::SectionReader {
    typedef LemonReader::SectionReader Parent;
    typedef _Graph Graph;
    typedef typename Graph::Edge Edge;
  public:
    
    /// \brief Constructor.
    ///
    /// Constructor for EdgeReader. It creates the EdgeReader and
    /// attach it into the given LemonReader. It will use the given
    /// edge id reader to give back the edges. The reader will read the 
    /// section only if the \c _name and the \c edges_name are the same. 
    template <typename _LabelReader>
    EdgeReader(LemonReader& _reader, const _LabelReader& _labelReader, 
	       const std::string& _name = std::string()) 
      : Parent(_reader), name(_name) {
      checkConcept<_reader_bits::ItemLabelReader<Edge>, _LabelReader>();
      edgeLabelReader.reset(new _reader_bits::
			 LabelReader<Edge, _LabelReader>(_labelReader));
    }

    /// \brief Destructor.
    ///
    /// Destructor for EdgeReader.
    virtual ~EdgeReader() {}
  private:
    EdgeReader(const EdgeReader&);
    void operator=(const EdgeReader&);

  public:

    /// \brief Add an edge reader command for the EdgeReader.
    ///
    /// Add an edge reader command for the EdgeReader.
    void readEdge(std::string label, Edge& item) {
      if (readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for edge: " << label;
	throw IoParameterError(msg.message());
      }
      readers.insert(make_pair(label, _reader_bits::ItemStore<Edge>(item)));
    }

  protected:

    /// \brief Gives back true when the SectionReader can process 
    /// the section with the given header line.
    ///
    /// It gives back true when the header line start with \c \@edges,
    /// and the header line's name and the reader's name are the same.
    virtual bool header(const std::string& line) {
      std::istringstream ls(line);
      std::string command;
      std::string id;
      ls >> command >> id;
      return command == "@edges" && name == id;
    }

    /// \brief Reader function of the section.
    ///
    /// It reads the content of the section.
    virtual void read(std::istream& is) {
      if (!edgeLabelReader->isLabelReader()) {
	throw DataFormatError("Cannot find edgeset or label map");
      }
      std::string line;
      while (getline(is, line)) {
	std::istringstream ls(line);
	std::string id;
	ls >> id;
	typename EdgeReaders::iterator it = readers.find(id);
	if (it != readers.end()) {
	  it->second.read(edgeLabelReader->read(ls));
	  it->second.touch();
	}	
      }
      for (typename EdgeReaders::iterator it = readers.begin();
	   it != readers.end(); ++it) {
	if (!it->second.touched()) {
	  ErrorMessage msg;
	  msg << "Edge not found in file: " << it->first;
	  throw IoParameterError(msg.message());
	}
      }
    }

    virtual void missing() {
      if (readers.empty()) return;
      ErrorMessage msg;
      msg << "Edges section not found in file: @edges " << name;
      throw IoParameterError(msg.message());
    }
    
  private:

    std::string name;

    typedef std::map<std::string, _reader_bits::ItemStore<Edge> > EdgeReaders;
    EdgeReaders readers;
    std::auto_ptr<_reader_bits::LabelReaderBase<Edge> > edgeLabelReader;
  };

  /// \ingroup section_io
  /// \brief SectionReader for reading labeled undirected edges.
  ///
  /// The undirected edges section's header line is \c \@uedges 
  /// \c uedges_name, but the \c uedges_name may be empty.
  ///
  /// Each line in the section contains the name of the undirected edge 
  /// and then the undirected edge id. 
  ///
  /// \relates LemonReader
  template <typename _Graph>
  class UEdgeReader : public LemonReader::SectionReader {
    typedef LemonReader::SectionReader Parent;
    typedef _Graph Graph;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::UEdge UEdge;
  public:
    
    /// \brief Constructor.
    ///
    /// Constructor for UEdgeReader. It creates the UEdgeReader and
    /// attach it into the given LemonReader. It will use the given
    /// undirected edge id reader to give back the edges. The reader will 
    /// read the section only if the \c _name and the \c uedges_name are 
    /// the same. 
    template <typename _LabelReader>
    UEdgeReader(LemonReader& _reader, const _LabelReader& _labelReader, 
	       const std::string& _name = std::string()) 
      : Parent(_reader), name(_name) {
      checkConcept<_reader_bits::ItemLabelReader<UEdge>, _LabelReader>();
      checkConcept<_reader_bits::ItemLabelReader<Edge>, _LabelReader>();
      uedgeLabelReader.reset(new _reader_bits::
			     LabelReader<UEdge, _LabelReader>(_labelReader));
      edgeLabelReader.reset(new _reader_bits::
			    LabelReader<Edge, _LabelReader>(_labelReader));
    }

    /// \brief Destructor.
    ///
    /// Destructor for UEdgeReader.
    virtual ~UEdgeReader() {}
  private:
    UEdgeReader(const UEdgeReader&);
    void operator=(const UEdgeReader&);

  public:

    /// \brief Add an undirected edge reader command for the UEdgeReader.
    ///
    /// Add an undirected edge reader command for the UEdgeReader.
    void readUEdge(std::string label, UEdge& item) {
      if (uedgeReaders.find(label) != uedgeReaders.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for undirected edge: " << label;
	throw IoParameterError(msg.message());
      }
      uedgeReaders.insert(make_pair(label, _reader_bits::
					ItemStore<UEdge>(item)));
    }

    /// \brief Add an edge reader command for the UEdgeReader.
    ///
    /// Add an edge reader command for the UEdgeReader.
    void readEdge(std::string label, Edge& item) {
      if (edgeReaders.find(label) != edgeReaders.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for edge: " << label;
	throw IoParameterError(msg.message());
      }
      edgeReaders.insert(make_pair(label, _reader_bits::ItemStore<Edge>(item)));
    }

  protected:

    /// \brief Gives back true when the SectionReader can process 
    /// the section with the given header line.
    ///
    /// It gives back true when the header line start with \c \@edges,
    /// and the header line's name and the reader's name are the same.
    virtual bool header(const std::string& line) {
      std::istringstream ls(line);
      std::string command;
      std::string id;
      ls >> command >> id;
      return command == "@uedges" && name == id;
    }

    /// \brief Reader function of the section.
    ///
    /// It reads the content of the section.
    virtual void read(std::istream& is) {
      if (!edgeLabelReader->isLabelReader()) {
	throw DataFormatError("Cannot find undirected edgeset or label map");
      }
      if (!uedgeLabelReader->isLabelReader()) {
	throw DataFormatError("Cannot find undirected edgeset or label map");
      }
      std::string line;
      while (getline(is, line)) {
	std::istringstream ls(line);
	std::string id;
	ls >> id;
	{
	  typename UEdgeReaders::iterator it = uedgeReaders.find(id);
	  if (it != uedgeReaders.end()) {
	    it->second.read(uedgeLabelReader->read(ls));
	    it->second.touch();
	    continue;
	  }	
	} {
	  typename EdgeReaders::iterator it = edgeReaders.find(id);
	  if (it != edgeReaders.end()) {
	    it->second.read(edgeLabelReader->read(ls));
	    it->second.touch();
	    continue;
	  }	
	}
      }
      for (typename EdgeReaders::iterator it = edgeReaders.begin();
	   it != edgeReaders.end(); ++it) {
	if (!it->second.touched()) {
	  ErrorMessage msg;
	  msg << "Edge not found in file: " << it->first;
	  throw IoParameterError(msg.message());
	}
      }
      for (typename UEdgeReaders::iterator it = uedgeReaders.begin();
	   it != uedgeReaders.end(); ++it) {
	if (!it->second.touched()) {
	  ErrorMessage msg;
	  msg << "UEdge not found in file: " << it->first;
	  throw IoParameterError(msg.message());
	}
      }
    }

    virtual void missing() {
      if (edgeReaders.empty() && uedgeReaders.empty()) return;
      ErrorMessage msg;
      msg << "UEdges section not found in file: @uedges " << name;
      throw IoParameterError(msg.message());
    }
    
  private:

    std::string name;

    typedef std::map<std::string, 
		     _reader_bits::ItemStore<UEdge> > UEdgeReaders;
    UEdgeReaders uedgeReaders;
    std::auto_ptr<_reader_bits::LabelReaderBase<UEdge> > uedgeLabelReader;

    typedef std::map<std::string, _reader_bits::ItemStore<Edge> > EdgeReaders;
    EdgeReaders edgeReaders;
    std::auto_ptr<_reader_bits::LabelReaderBase<Edge> > edgeLabelReader;
  };

  /// \ingroup section_io
  /// \brief SectionReader for attributes.
  ///
  /// The lemon format can store multiple attribute set. Each set has
  /// the header line \c \@attributes \c attributeset_name, but the 
  /// attributeset_name may be empty.
  ///
  /// The attributeset section contains several lines. Each of them starts
  /// with an attribute and then a the value for the id.
  ///
  /// \relates LemonReader
  template <typename _Traits = DefaultReaderTraits>
  class AttributeReader : public LemonReader::SectionReader {
    typedef LemonReader::SectionReader Parent;
    typedef _Traits Traits; 
  public:
    /// \brief Constructor.
    ///
    /// Constructor for AttributeReader. It creates the AttributeReader and
    /// attach it into the given LemonReader. The reader process a section
    /// only if the \c section_name and the \c _name are the same.
    AttributeReader(LemonReader& _reader, 
		    const std::string& _name = std::string()) 
      : Parent(_reader), name(_name) {}

    /// \brief Destructor.
    ///
    /// Destructor for AttributeReader.
    virtual ~AttributeReader() {
      for (typename Readers::iterator it = readers.begin(); 
	   it != readers.end(); ++it) {
	delete it->second;
      }
    }

  private:
    AttributeReader(const AttributeReader&);
    void operator=(AttributeReader&);

  public:
    /// \brief Add an attribute reader command for the reader.
    ///
    /// Add an attribute reader command for the reader.
    template <typename Value>
    AttributeReader& readAttribute(const std::string& label, Value& value) {
      return readAttribute<typename Traits::template Reader<Value> >
	(label, value);
    }

    /// \brief Add an attribute reader command for the reader.
    ///
    /// Add an attribute reader command for the reader.
    template <typename ItemReader, typename Value>
    AttributeReader& readAttribute(const std::string& label, Value& value,
				   const ItemReader& ir = ItemReader()) {
      checkConcept<_reader_bits::ItemReader<Value>, ItemReader>();
      if (readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for attribute: " << label;
	throw IoParameterError(msg.message());
      }
      readers.insert(make_pair(label, new _reader_bits::
			       ValueReader<Value, ItemReader>(value, ir)));
      return *this;
    }

  protected:

    /// \brief Gives back true when the SectionReader can process 
    /// the section with the given header line.
    ///
    /// It gives back true when the header line start with \c \@attributes,
    /// and the header line's id and the attributeset's id are the same.
    bool header(const std::string& line) {
      std::istringstream ls(line);
      std::string command;
      std::string id;
      ls >> command >> id;
      return command == "@attributes" && name == id;
    }

    /// \brief Reader function of the section.
    ///
    /// It reads the content of the section.
    void read(std::istream& is) {
      std::string line;
      while (getline(is, line)) {
	std::istringstream ls(line);
	std::string id;
	ls >> id;
	typename Readers::iterator it = readers.find(id);
	if (it != readers.end()) {
	  it->second->read(ls);
 	  it->second->touch();
	}
      }
      for (typename Readers::iterator it = readers.begin();
	   it != readers.end(); ++it) {
	if (!it->second->touched()) {
	  ErrorMessage msg;
	  msg << "Attribute not found in file: " << it->first;
	  throw IoParameterError(msg.message());
	}	
      }
    }    
    
    virtual void missing() {
      if (readers.empty()) return;
      ErrorMessage msg;
      msg << "Attribute section not found in file: @attributes " << name;
      throw IoParameterError(msg.message());
    }

  private:
    std::string name;

    typedef std::map<std::string, _reader_bits::ValueReaderBase*> Readers;
    Readers readers;  
  };

  /// \ingroup section_io
  /// \brief SectionReader for reading extra node maps.
  ///
  /// The lemon format can store maps in the nodeset sections. This
  /// class let you make distinict section to store maps.  The main
  /// purpose of this class is a logical separation of some maps. The
  /// other useful application could be to store paths in node maps.
  ///
  /// The first line of the section contains the names of the maps
  /// separated with white spaces. Each next line describes an item
  /// in the itemset, and contains in the first column the label of
  /// the item and then the mapped values for each map.
  ///
  /// \relates LemonReader
  template <typename _Graph, typename _Traits = DefaultReaderTraits>
  class NodeMapReader : public LemonReader::SectionReader {
    typedef LemonReader::SectionReader Parent;
  public:

    typedef _Graph Graph;
    typedef typename Graph::Node Node;
    typedef _Traits Traits;
    typedef typename Traits::Skipper DefaultSkipper;

    /// \brief Constructor.
    ///
    /// Constructor for NodeMapReader. It creates the NodeMapReader and
    /// attach it into the given LemonReader. The reader will read
    /// the section when the \c section_name and the \c _name are the same.
    template <typename _LabelReader>
    NodeMapReader(LemonReader& _reader, 
		  const Graph& _graph, 
		  const _LabelReader& _labelReader,
		  const std::string& _name = std::string(),
		  const DefaultSkipper& _skipper = DefaultSkipper()) 
      : Parent(_reader), graph(_graph), name(_name), skipper(_skipper) {
      labelReader.reset(new _reader_bits::
			LabelReader<Node, _LabelReader>(_labelReader));
    } 


    /// \brief Destructor.
    ///
    /// Destructor for NodeMapReader.
    virtual ~NodeMapReader() {
      for (typename MapReaders::iterator it = readers.begin(); 
	   it != readers.end(); ++it) {
	delete it->second;
      }
    }

  private:
    NodeMapReader(const NodeMapReader&);
    void operator=(const NodeMapReader&);
  
  public:

    /// \brief Add a new node map reader command for the reader.
    ///
    /// Add a new node map reader command for the reader.
    template <typename Map>
    NodeMapReader& readNodeMap(std::string label, Map& map) {
      return _readMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    template <typename Map>
    NodeMapReader& readNodeMap(std::string label, const Map& map) {
      return _readMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    /// \brief Add a new node map reader command for the reader.
    ///
    /// Add a new node map reader command for the reader.
    template <typename ItemReader, typename Map>
    NodeMapReader& readNodeMap(std::string label, Map& map, 
			       const ItemReader& ir = ItemReader()) {
      return _readMap<ItemReader, Map, typename _reader_bits::Arg<Map>::Type>
	(label, map, ir);
    }

    template <typename ItemReader, typename Map>
    NodeMapReader& readNodeMap(std::string label, const Map& map, 
			       const ItemReader& ir = ItemReader()) {
      return _readMap<ItemReader, Map, typename _reader_bits::Arg<Map>::Type>
	(label, map, ir);
    }

  private:

    template <typename ItemReader, typename Map, typename MapParameter>
    NodeMapReader& _readMap(std::string label, MapParameter map, 
			   const ItemReader& ir = ItemReader()) {
      checkConcept<concepts::WriteMap<Node, typename Map::Value>, Map>();
      checkConcept<_reader_bits::ItemReader<typename Map::Value>, ItemReader>();
      if (readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for map: " << label;
	throw IoParameterError(msg.message());
      }      
      readers.insert(
	make_pair(label, new _reader_bits::
		  MapReader<Node, Map, ItemReader>(map, ir)));
      return *this;
    }

  public:

    /// \brief Add a new node map skipper command for the reader.
    ///
    /// Add a new node map skipper command for the reader.
    template <typename ItemReader>
    NodeMapReader& skipNodeMap(std::string label, 
			       const ItemReader& ir = ItemReader()) {
      if (readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for map: " << label;
	throw IoParameterError(msg.message());
      }
      readers.insert(make_pair(label, new _reader_bits::
			       SkipReader<Node, ItemReader>(ir)));
      return *this;
    }

  protected:

    /// \brief Gives back true when the SectionReader can process 
    /// the section with the given header line.
    ///
    /// It gives back true when the header line starts with \c \@mapset,
    /// and the header line's name and the mapset's name are the same.
    virtual bool header(const std::string& line) {
      std::istringstream ls(line);
      std::string command;
      std::string id;
      ls >> command >> id;
      return command == "@nodemaps" && name == id;
    }

    /// \brief Reader function of the section.
    ///
    /// It reads the content of the section.
    virtual void read(std::istream& is) {
      std::vector<_reader_bits::MapReaderBase<Node>* > index;
      std::string line;

      {
        getline(is, line);
        std::istringstream ls(line);
        std::string id;
        while (ls >> id) {
          typename MapReaders::iterator it = readers.find(id);
          if (it != readers.end()) {
            it->second->touch();
            index.push_back(it->second);
          } else {
            index.push_back(&skipper);
          }
        }
      }
      for (typename MapReaders::iterator it = readers.begin();
	   it != readers.end(); ++it) {
	if (!it->second->touched()) {
	  ErrorMessage msg;
	  msg << "Map not found in file: " << it->first;
	  throw IoParameterError(msg.message());
	}
      }
      while (getline(is, line)) {	
	std::istringstream ls(line);
	Node node = labelReader->read(ls);
	for (int i = 0; i < int(index.size()); ++i) {
	  index[i]->read(ls, node);
	}
      }
    }

    virtual void missing() {
      if (readers.empty()) return;
      ErrorMessage msg;
      msg << "NodeMap section not found in file: @nodemaps " << name;
      throw IoParameterError(msg.message());
    }

  private:

    typedef std::map<std::string, _reader_bits::MapReaderBase<Node>*> MapReaders;
    MapReaders readers;
   
    const Graph& graph;   
    std::string name;
    _reader_bits::SkipReader<Node, DefaultSkipper> skipper;
    std::auto_ptr<_reader_bits::LabelReaderBase<Node> > labelReader;

  };

  /// \ingroup section_io
  /// \brief SectionReader for reading extra edge maps.
  ///
  /// The lemon format can store maps in the edgeset sections. This
  /// class let you make distinict section to store maps.  The main
  /// purpose of this class is a logical separation of some maps. The
  /// other useful application could be to store paths in edge maps.
  ///
  /// The first line of the section contains the names of the maps
  /// separated with white spaces. Each next line describes an item
  /// in the itemset, and contains in the first column the label of
  /// the item and then the mapped values for each map.
  ///
  /// \relates LemonReader
  template <typename _Graph, typename _Traits = DefaultReaderTraits>
  class EdgeMapReader : public LemonReader::SectionReader {
    typedef LemonReader::SectionReader Parent;
  public:

    typedef _Graph Graph;
    typedef typename Graph::Edge Edge;
    typedef _Traits Traits;
    typedef typename Traits::Skipper DefaultSkipper;

    /// \brief Constructor.
    ///
    /// Constructor for EdgeMapReader. It creates the EdgeMapReader and
    /// attach it into the given LemonReader. The reader will read
    /// the section when the \c section_name and the \c _name are the same.
    template <typename _LabelReader>
    EdgeMapReader(LemonReader& _reader, 
		   const Graph& _graph, 
		   const _LabelReader& _labelReader,
		   const std::string& _name = std::string(),
		   const DefaultSkipper& _skipper = DefaultSkipper()) 
      : Parent(_reader), graph(_graph), name(_name), skipper(_skipper) {
      labelReader.reset(new _reader_bits::
			LabelReader<Edge, _LabelReader>(_labelReader));
    } 


    /// \brief Destructor.
    ///
    /// Destructor for EdgeMapReader.
    virtual ~EdgeMapReader() {
      for (typename MapReaders::iterator it = readers.begin(); 
	   it != readers.end(); ++it) {
	delete it->second;
      }
    }

  private:
    EdgeMapReader(const EdgeMapReader&);
    void operator=(const EdgeMapReader&);
  
  public:

    /// \brief Add a new edge map reader command for the reader.
    ///
    /// Add a new edge map reader command for the reader.
    template <typename Map>
    EdgeMapReader& readEdgeMap(std::string label, Map& map) {
      return _readMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    template <typename Map>
    EdgeMapReader& readEdgeMap(std::string label, const Map& map) {
      return _readMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    /// \brief Add a new edge map reader command for the reader.
    ///
    /// Add a new edge map reader command for the reader.
    template <typename ItemReader, typename Map>
    EdgeMapReader& readEdgeMap(std::string label, Map& map, 
			  const ItemReader& ir = ItemReader()) {
      return _readMap<ItemReader, Map, typename _reader_bits::Arg<Map>::Type>
	(label, map, ir);
    }

    template <typename ItemReader, typename Map>
    EdgeMapReader& readEdgeMap(std::string label, const Map& map, 
			  const ItemReader& ir = ItemReader()) {
      return _readMap<ItemReader, Map, typename _reader_bits::Arg<Map>::Type>
	(label, map, ir);
    }

  private:

    template <typename ItemReader, typename Map, typename MapParameter>
    EdgeMapReader& _readMap(std::string label, MapParameter map, 
				const ItemReader& ir = ItemReader()) {
      checkConcept<concepts::WriteMap<Edge, typename Map::Value>, Map>();
      checkConcept<_reader_bits::ItemReader<typename Map::Value>, ItemReader>();
      if (readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for map: " << label;
	throw IoParameterError(msg.message());
      }      
      readers.insert(
	make_pair(label, new _reader_bits::
		  MapReader<Edge, Map, ItemReader>(map, ir)));
      return *this;
    }

  public:

    /// \brief Add a new edge map skipper command for the reader.
    ///
    /// Add a new edge map skipper command for the reader.
    template <typename ItemReader>
    EdgeMapReader& skipEdgeMap(std::string label, 
			  const ItemReader& ir = ItemReader()) {
      if (readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for map: " << label;
	throw IoParameterError(msg.message());
      }
      readers.insert(make_pair(label, new _reader_bits::
			       SkipReader<Edge, ItemReader>(ir)));
      return *this;
    }

  protected:

    /// \brief Gives back true when the SectionReader can process 
    /// the section with the given header line.
    ///
    /// It gives back true when the header line starts with \c \@mapset,
    /// and the header line's name and the mapset's name are the same.
    virtual bool header(const std::string& line) {
      std::istringstream ls(line);
      std::string command;
      std::string id;
      ls >> command >> id;
      return (command == "@edgemaps" || command == "@uedgemaps") && name == id;
    }

    /// \brief Reader function of the section.
    ///
    /// It reads the content of the section.
    virtual void read(std::istream& is) {
      std::vector<_reader_bits::MapReaderBase<Edge>* > index;
      std::string line;

      {
        getline(is, line);
        std::istringstream ls(line);
        std::string id;
        while (ls >> id) {
          typename MapReaders::iterator it = readers.find(id);
          if (it != readers.end()) {
            it->second->touch();
            index.push_back(it->second);
          } else {
            index.push_back(&skipper);
          }
        }
      }
      for (typename MapReaders::iterator it = readers.begin();
	   it != readers.end(); ++it) {
	if (!it->second->touched()) {
	  ErrorMessage msg;
	  msg << "Map not found in file: " << it->first;
	  throw IoParameterError(msg.message());
	}
      }
      while (getline(is, line)) {	
	std::istringstream ls(line);
	Edge edge = labelReader->read(ls);
	for (int i = 0; i < int(index.size()); ++i) {
	  index[i]->read(ls, edge);
	}
      }
    }

    virtual void missing() {
      if (readers.empty()) return;
      ErrorMessage msg;
      msg << "EdgeMap section not found in file: @edgemaps " << name;
      throw IoParameterError(msg.message());
    }

  private:

    typedef std::map<std::string, _reader_bits::MapReaderBase<Edge>*> MapReaders;
    MapReaders readers;
   
    const Graph& graph;   
    std::string name;
    _reader_bits::SkipReader<Edge, DefaultSkipper> skipper;
    std::auto_ptr<_reader_bits::LabelReaderBase<Edge> > labelReader;

  };

  /// \ingroup section_io
  /// \brief SectionReader for reading extra undirected edge maps.
  ///
  /// The lemon format can store maps in the uedgeset sections. This
  /// class let you make distinict section to store maps.  The main
  /// purpose of this class is a logical separation of some maps. The
  /// other useful application could be to store paths in undirected
  /// edge maps.
  ///
  /// The first line of the section contains the names of the maps
  /// separated with white spaces. Each next line describes an item
  /// in the itemset, and contains in the first column the label of
  /// the item and then the mapped values for each map.
  ///
  /// \relates LemonReader
  template <typename _Graph, typename _Traits = DefaultReaderTraits>
  class UEdgeMapReader : public LemonReader::SectionReader {
    typedef LemonReader::SectionReader Parent;
  public:

    typedef _Graph Graph;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::UEdge UEdge;
    typedef _Traits Traits;
    typedef typename Traits::Skipper DefaultSkipper;

    /// \brief Constructor.
    ///
    /// Constructor for UEdgeMapReader. It creates the UEdgeMapReader and
    /// attach it into the given LemonReader. The reader will read
    /// the section when the \c section_name and the \c _name are the same.
    template <typename _LabelReader>
    UEdgeMapReader(LemonReader& _reader, const Graph& _graph, 
		   const _LabelReader& _labelReader,
		   const std::string& _name = std::string(),
		   const DefaultSkipper& _skipper = DefaultSkipper()) 
      : Parent(_reader), graph(_graph), name(_name), skipper(_skipper) {
      labelReader.reset(new _reader_bits::
			LabelReader<UEdge, _LabelReader>(_labelReader));
    } 


    /// \brief Destructor.
    ///
    /// Destructor for UEdgeMapReader.
    virtual ~UEdgeMapReader() {
      for (typename MapReaders::iterator it = readers.begin(); 
	   it != readers.end(); ++it) {
	delete it->second;
      }
    }

  private:
    UEdgeMapReader(const UEdgeMapReader&);
    void operator=(const UEdgeMapReader&);
  
  public:

    /// \brief Add a new undirected edge map reader command for the
    /// reader.
    ///
    /// Add a new undirected edge map reader command for the reader.
    template <typename Map>
    UEdgeMapReader& readUEdgeMap(std::string label, Map& map) {
      return _readMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    template <typename Map>
    UEdgeMapReader& readUEdgeMap(std::string label, const Map& map) {
      return _readMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    /// \brief Add a new undirected edge map reader command for the
    /// reader.
    ///
    /// Add a new undirected edge map reader command for the reader.
    template <typename ItemReader, typename Map>
    UEdgeMapReader& readUEdgeMap(std::string label, Map& map, 
			  const ItemReader& ir = ItemReader()) {
      return _readMap<ItemReader, Map, typename _reader_bits::Arg<Map>::Type>
	(label, map, ir);
    }

    template <typename ItemReader, typename Map>
    UEdgeMapReader& readUEdgeMap(std::string label, const Map& map, 
			  const ItemReader& ir = ItemReader()) {
      return _readMap<ItemReader, Map, typename _reader_bits::Arg<Map>::Type>
	(label, map, ir);
    }

  private:

    template <typename ItemReader, typename Map, typename MapParameter>
    UEdgeMapReader& _readMap(std::string label, MapParameter map, 
				const ItemReader& ir = ItemReader()) {
      checkConcept<concepts::WriteMap<Edge, typename Map::Value>, Map>();
      checkConcept<_reader_bits::ItemReader<typename Map::Value>, ItemReader>();
      if (readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for map: " << label;
	throw IoParameterError(msg.message());
      }      
      readers.insert(
	make_pair(label, new _reader_bits::
		  MapReader<UEdge, Map, ItemReader>(map, ir)));
      return *this;
    }

  public:

    /// \brief Add a new undirected edge map skipper command for the
    /// reader.
    ///
    /// Add a new undirected edge map skipper command for the reader.
    template <typename ItemReader>
    UEdgeMapReader& skipUEdgeMap(std::string label, 
			  const ItemReader& ir = ItemReader()) {
      if (readers.find(label) != readers.end()) {
	ErrorMessage msg;
	msg << "Multiple read rule for map: " << label;
	throw IoParameterError(msg.message());
      }
      readers.insert(make_pair(label, new _reader_bits::
			       SkipReader<Edge, ItemReader>(ir)));
      return *this;
    }

    /// \brief Add a new directed edge map reader command for the reader.
    ///
    /// Add a new directed edge map reader command for the reader.
    template <typename Map>
    UEdgeMapReader& readEdgeMap(std::string label, Map& map) {
      return _readDirMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    template <typename Map>
    UEdgeMapReader& readEdgeMap(std::string label, const Map& map) {
      return _readDirMap<
	typename Traits::template Reader<typename Map::Value>, Map,
	typename _reader_bits::Arg<Map>::Type>(label, map);
    }

    /// \brief Add a new directed edge map reader command for the reader.
    ///
    /// Add a new directed edge map reader command for the reader.
    template <typename ItemReader, typename Map>
    UEdgeMapReader& readEdgeMap(std::string label, Map& map, 
				    const ItemReader& ir = ItemReader()) {
      return _readDirMap<ItemReader, Map, 
        typename _reader_bits::Arg<Map>::Type>(label, map, ir);
    }

    template <typename ItemReader, typename Map>
    UEdgeMapReader& readEdgeMap(std::string label, const Map& map, 
				    const ItemReader& ir = ItemReader()) {
      return _readDirMap<ItemReader, Map, 
        typename _reader_bits::Arg<Map>::Type>(label, map, ir);
    }

  private:

    template <typename ItemReader, typename Map, typename MapParameter>
    UEdgeMapReader& _readDirMap(std::string label, MapParameter map,
				    const ItemReader& ir = ItemReader()) { 
      checkConcept<_reader_bits::ItemReader<typename Map::Value>, ItemReader>();
      checkConcept<concepts::WriteMap<Edge, typename Map::Value>, Map>();
      readUEdgeMap("+" + label, 
                   _reader_bits::forwardComposeMap(graph, map), ir);
      readUEdgeMap("-" + label, 
                   _reader_bits::backwardComposeMap(graph, map), ir);
      return *this;      
    }

  public:

    /// \brief Add a new directed edge map skipper command for the reader.
    ///
    /// Add a new directed edge map skipper command for the reader.
    template <typename ItemReader>
    UEdgeMapReader& skipEdgeMap(std::string label, 
                                const ItemReader& ir = ItemReader()) {
      skipUEdgeMap("+" + label, ir);
      skipUEdgeMap("-" + label, ir);
      return *this;
    }

  protected:

    /// \brief Gives back true when the SectionReader can process 
    /// the section with the given header line.
    ///
    /// It gives back true when the header line starts with \c \@mapset,
    /// and the header line's name and the mapset's name are the same.
    virtual bool header(const std::string& line) {
      std::istringstream ls(line);
      std::string command;
      std::string id;
      ls >> command >> id;
      return (command == "@edgemaps" || command == "@uedgemaps") && name == id;
    }

    /// \brief Reader function of the section.
    ///
    /// It reads the content of the section.
    virtual void read(std::istream& is) {
      std::vector<_reader_bits::MapReaderBase<UEdge>* > index;
      std::string line;

      {
        getline(is, line);
        std::istringstream ls(line);
        std::string id;
        while (ls >> id) {
          typename MapReaders::iterator it = readers.find(id);
          if (it != readers.end()) {
            it->second->touch();
            index.push_back(it->second);
          } else {
            index.push_back(&skipper);
          }
        }
      }
      for (typename MapReaders::iterator it = readers.begin();
	   it != readers.end(); ++it) {
	if (!it->second->touched()) {
	  ErrorMessage msg;
	  msg << "Map not found in file: " << it->first;
	  throw IoParameterError(msg.message());
	}
      }
      while (getline(is, line)) {	
	std::istringstream ls(line);
	UEdge uedge = labelReader->read(ls);
	for (int i = 0; i < int(index.size()); ++i) {
	  index[i]->read(ls, uedge);
	}
      }
    }

    virtual void missing() {
      if (readers.empty()) return;
      ErrorMessage msg;
      msg << "UEdgeMap section not found in file: @uedgemaps " << name;
      throw IoParameterError(msg.message());
    }

  private:

    const Graph& graph;   
    std::string name;

    typedef std::map<std::string, 
		     _reader_bits::MapReaderBase<UEdge>*> MapReaders;
   
    MapReaders readers;
    _reader_bits::SkipReader<UEdge, DefaultSkipper> skipper;

    std::auto_ptr<_reader_bits::LabelReaderBase<UEdge> > labelReader;

  };

  /// \ingroup section_io
  /// \brief SectionReader for retrieve what is in the file.
  ///
  /// SectionReader for retrieve what is in the file. If you want
  /// to know which sections, maps and items are in the file
  /// use the next code:
  ///\code
  /// LemonReader reader("input.lgf");
  /// ContentReader content(reader);
  /// reader.run();
  ///\endcode
  class ContentReader : public LemonReader::SectionReader {
    typedef LemonReader::SectionReader Parent;
  public:
    /// \brief Constructor.
    ///
    /// Constructor for
    ContentReader(LemonReader& _reader) : Parent(_reader) {}

    /// \brief Desctructor.
    ///
    /// Desctructor.
    virtual ~ContentReader() {}

    /// \brief Gives back how many nodesets are in the file.
    ///
    /// Gives back how many nodesets are in the file.
    int nodeSetNum() const {
      return nodesets.size();
    }

    /// \brief Gives back the name of nodeset on the indiced position.
    ///
    /// Gives back the name of nodeset on the indiced position.
    std::string nodeSetName(int index) const {
      return nodesets[index].name;
    }

    /// \brief Gives back the map names of nodeset on the indiced position.
    ///
    /// Gives back the map names of nodeset on the indiced position.
    const std::vector<std::string>& nodeSetMaps(int index) const {
      return nodesets[index].items;
    }

    /// \brief Gives back how many edgesets are in the file.
    ///
    /// Gives back how many edgesets are in the file.
    int edgeSetNum() const {
      return edgesets.size();
    }

    /// \brief Gives back the name of edgeset on the indiced position.
    ///
    /// Gives back the name of edgeset on the indiced position.
    std::string edgeSetName(int index) const {
      return edgesets[index].name;
    }

    /// \brief Gives back the map names of edgeset on the indiced position.
    ///
    /// Gives back the map names of edgeset on the indiced position.
    const std::vector<std::string>& edgeSetMaps(int index) const {
      return edgesets[index].items;
    }

    /// \brief Gives back how many undirected edgesets are in the file.
    ///
    /// Gives back how many undirected edgesets are in the file.
    int uEdgeSetNum() const {
      return uedgesets.size();
    }

    /// \brief Gives back the name of undirected edgeset on the indiced 
    /// position.
    ///
    /// Gives back the name of undirected edgeset on the indiced position.
    std::string uEdgeSetName(int index) const {
      return uedgesets[index].name;
    }

    /// \brief Gives back the map names of undirected edgeset on the indiced 
    /// position.
    ///
    /// Gives back the map names of undirected edgeset on the indiced position.
    const std::vector<std::string>& uEdgeSetMaps(int index) const {
      return uedgesets[index].items;
    }

    /// \brief Gives back how many labeled nodes section are in the file.
    ///
    /// Gives back how many labeled nodes section are in the file.
    int nodesNum() const {
      return nodes.size();
    }

    /// \brief Gives back the name of labeled nodes section on the indiced 
    /// position.
    ///
    /// Gives back the name of labeled nodes section on the indiced position.
    std::string nodesName(int index) const {
      return nodes[index].name;
    }

    /// \brief Gives back the names of the labeled nodes in the indiced 
    /// section.
    ///
    /// Gives back the names of the labeled nodes in the indiced section.
    const std::vector<std::string>& nodesItems(int index) const {
      return nodes[index].items;
    }

    /// \brief Gives back how many labeled edges section are in the file.
    ///
    /// Gives back how many labeled edges section are in the file.
    int edgesNum() const {
      return edges.size();
    }

    /// \brief Gives back the name of labeled edges section on the indiced 
    /// position.
    ///
    /// Gives back the name of labeled edges section on the indiced position.
    std::string edgesName(int index) const {
      return edges[index].name;
    }

    /// \brief Gives back the names of the labeled edges in the indiced 
    /// section.
    ///
    /// Gives back the names of the labeled edges in the indiced section.
    const std::vector<std::string>& edgesItems(int index) const {
      return edges[index].items;
    }
 
    /// \brief Gives back how many labeled undirected edges section are 
    /// in the file.
    ///
    /// Gives back how many labeled undirected edges section are in the file.
    int uEdgesNum() const {
      return uedges.size();
    }

    /// \brief Gives back the name of labeled undirected edges section 
    /// on the indiced position.
    ///
    /// Gives back the name of labeled undirected edges section on the 
    /// indiced position.
    std::string uEdgesName(int index) const {
      return uedges[index].name;
    }

    /// \brief Gives back the names of the labeled undirected edges in 
    /// the indiced section.
    ///
    /// Gives back the names of the labeled undirected edges in the 
    /// indiced section.
    const std::vector<std::string>& uEdgesItems(int index) const {
      return uedges[index].items;
    }

 
    /// \brief Gives back how many attributes section are in the file.
    ///
    /// Gives back how many attributes section are in the file.
    int attributesNum() const {
      return attributes.size();
    }

    /// \brief Gives back the name of attributes section on the indiced 
    /// position.
    ///
    /// Gives back the name of attributes section on the indiced position.
    std::string attributesName(int index) const {
      return attributes[index].name;
    }

    /// \brief Gives back the names of the attributes in the indiced section.
    ///
    /// Gives back the names of the attributes in the indiced section.
    const std::vector<std::string>& attributesItems(int index) const {
      return attributes[index].items;
    }

    const std::vector<std::string>& otherSections() const {
      return sections;
    }

  protected:
    
    /// \brief Gives back true when the SectionReader can process 
    /// the section with the given header line.
    ///
    /// It gives back true when the section is common section.
    bool header(const std::string& line) {
      std::istringstream ls(line);
      std::string command, name;
      ls >> command >> name;
      if (command == "@nodeset") {
	current = command;
	nodesets.push_back(SectionInfo(name));
      } else if (command == "@edgeset") {
	current = command;
	edgesets.push_back(SectionInfo(name));
      } else if (command == "@uedgeset") {
	current = command;
	uedgesets.push_back(SectionInfo(name));
      } else if (command == "@nodes") {
	current = command;
	nodes.push_back(SectionInfo(name));
      } else if (command == "@edges") {
	current = command;
	edges.push_back(SectionInfo(name));
      } else if (command == "@uedges") {
	current = command;
	uedges.push_back(SectionInfo(name));
      } else if (command == "@attributes") {
	current = command;
	attributes.push_back(SectionInfo(name));
      } else {
	sections.push_back(line);
	return false;
      }
      return true;
    }

    /// \brief Retrieve the items from various sections.
    ///
    /// Retrieve the items from various sections.
    void read(std::istream& is) {
      if (current == "@nodeset") {
	readMapNames(is, nodesets.back().items);
      } else if (current == "@edgeset") {
	readMapNames(is, edgesets.back().items);
      } else if (current == "@uedgeset") {
	readMapNames(is, uedgesets.back().items);
      } else if (current == "@nodes") {
	readItemNames(is, nodes.back().items);
      } else if (current == "@edges") {
	readItemNames(is, edges.back().items);
      } else if (current == "@uedges") {
	readItemNames(is, uedges.back().items);
      } else if (current == "@attributes") {
	readItemNames(is, attributes.back().items);
      }
    }    

  private:

    void readMapNames(std::istream& is, std::vector<std::string>& maps) {
      std::string line, name;
      std::getline(is, line);
      std::istringstream ls(line);
      while (ls >> name) {
	maps.push_back(name);
      }
      while (getline(is, line));
    }

    void readItemNames(std::istream& is, std::vector<std::string>& maps) {
      std::string line, name;
      while (std::getline(is, line)) {
	std::istringstream ls(line);
	ls >> name;
	maps.push_back(name);
      }
    }

    struct SectionInfo {
      std::string name;
      std::vector<std::string> items;

      SectionInfo(const std::string& _name) : name(_name) {}
    };

    std::vector<SectionInfo> nodesets;
    std::vector<SectionInfo> edgesets;
    std::vector<SectionInfo> uedgesets;

    std::vector<SectionInfo> nodes;
    std::vector<SectionInfo> edges;
    std::vector<SectionInfo> uedges;

    std::vector<SectionInfo> attributes;

    std::vector<std::string> sections;

    std::string current;

  };

}
#endif
