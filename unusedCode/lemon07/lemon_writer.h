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
///\brief Lemon Format writer.

#ifndef LEMON_LEMON_WRITER_H
#define LEMON_LEMON_WRITER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>

#include <lemon/error.h>
#include <lemon/bits/invalid.h>
#include <lemon/graph_utils.h>
#include <lemon/bits/item_writer.h>
#include <lemon/bits/utility.h>
#include <lemon/maps.h>
#include <lemon/dim2.h>

#include <lemon/concept_check.h>
#include <lemon/concepts/maps.h>


namespace lemon {

  namespace _writer_bits {
    
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

    template <typename Map>
    struct ComposeLess {
      ComposeLess(const Map& _map) : map(_map), less() {}

      bool operator()(const typename Map::Key& p, 
                      const typename Map::Key& q) const {
	return less(map[p], map[q]);
      }
      const Map& map;
      Less<typename Map::Value> less;
    };

    template <typename UGraph, typename Map>
    struct UEdgeComposeLess {
      UEdgeComposeLess(const UGraph& _ugraph, const Map& _map) 
	: ugraph(_ugraph), map(_map), less() {}

      bool operator()(const typename UGraph::Edge& p, 
                      const typename UGraph::Edge& q) const {
	return p != q ? less(map[p], map[q]) : 
	  (!ugraph.direction(p) && ugraph.direction(q));
      }

      const UGraph& ugraph;
      const Map& map;
      Less<typename Map::Value> less;
    };

    template <typename Item>
    class ItemLabelWriter {
    public:

      bool isLabelWriter() { return true; }

      void writeLabel(std::ostream&, const Item&) {}
      
      template <class _ItemLabelWriter>
      struct Constraints {
	void constraints() {
	  bool b = writer.isLabelWriter();
	  ignore_unused_variable_warning(b);
	  writer.writeLabel(os, item);
	}
	_ItemLabelWriter& writer;
	std::ostream& os;
	const Item& item;
      };

    };

    template <typename Item>
    class ItemWriter {
    public:

      void write(std::ostream&, const Item&) {}
      
      template <class _ItemWriter>
      struct Constraints {
	void constraints() {
	  writer.write(os, item);
	}
	_ItemWriter& writer;
	std::ostream& os;
	const Item& item;
      };

    };

    template <typename Map>
    struct Ref { typedef const Map& Type; };

    template <typename Graph, typename Map>
    class ForwardComposeMap {
    public:
      typedef typename Graph::UEdge Key;
      typedef typename Map::Value Value;

      ForwardComposeMap(const Graph& _graph, const Map& _map) 
	: graph(_graph), map(_map) {}
      
      Value operator[](const Key& key) const {
	return map[graph.direct(key, true)];
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
    class BackwardComposeMap {
    public:
      typedef typename Graph::UEdge Key;
      typedef typename Map::Value Value;

      BackwardComposeMap(const Graph& _graph, const Map& _map) 
	: graph(_graph), map(_map) {}
      
      Value operator[](const Key& key) const {
	return map[graph.direct(key, false)];
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
    struct Ref<ForwardComposeMap<Graph, Map> > { 
      typedef ForwardComposeMap<Graph, Map> Type;
    };

    template <typename Graph, typename Map>
    struct Ref<BackwardComposeMap<Graph, Map> > { 
      typedef BackwardComposeMap<Graph, Map> Type; 
    };

    template <typename Map>
    struct Ref<dim2::XMap<Map> > { 
      typedef dim2::XMap<Map> Type;
    };
    template <typename Map>
    struct Ref<dim2::ConstXMap<Map> > { 
      typedef dim2::ConstXMap<Map> Type;
    };

    template <typename Map>
    struct Ref<dim2::YMap<Map> > { 
      typedef dim2::YMap<Map> Type;
    };
    template <typename Map>
    struct Ref<dim2::ConstYMap<Map> > { 
      typedef dim2::ConstYMap<Map> Type;
    };


    template <typename _Item>    
    class MapWriterBase {
    public:
      typedef _Item Item;

      virtual ~MapWriterBase() {}

      virtual void write(std::ostream& os, const Item& item) const = 0;
      virtual void sort(std::vector<Item>&) const = 0;
    };


    template <typename _Item, typename _Map, typename _Writer>
    class MapWriter : public MapWriterBase<_Item> {
    public:
      typedef _Map Map;
      typedef _Writer Writer;
      typedef typename Writer::Value Value;
      typedef _Item Item;
      
      typename _writer_bits::Ref<Map>::Type map;
      Writer writer;

      MapWriter(const Map& _map, const Writer& _writer) 
	: map(_map), writer(_writer) {}

      virtual ~MapWriter() {}

      virtual void write(std::ostream& os, const Item& item) const {
	Value value = map[item];
	writer.write(os, value);
      }

      virtual void sort(std::vector<Item>& items) const {
        ComposeLess<Map> less(map);
        std::sort(items.begin(), items.end(), less);
      }

    };

    template <typename _UGraph>    
    class UEdgeMapWriterBase {
    public:
      typedef typename _UGraph::Edge Edge;
      typedef typename _UGraph::UEdge UEdge;

      typedef UEdge Item;

      virtual ~UEdgeMapWriterBase() {}

      virtual void write(std::ostream& os, const Item& item) const = 0;
      virtual void sort(const _UGraph&, std::vector<Edge>&) const = 0;
      virtual void sort(std::vector<UEdge>&) const = 0;
    };


    template <typename _UGraph, typename _Map, typename _Writer>
    class UEdgeMapWriter : public UEdgeMapWriterBase<_UGraph> {
    public:
      typedef _Map Map;
      typedef _Writer Writer;
      typedef typename Writer::Value Value;

      typedef typename _UGraph::Edge Edge;
      typedef typename _UGraph::UEdge UEdge;
      typedef UEdge Item;
      
      typename _writer_bits::Ref<Map>::Type map;
      Writer writer;

      UEdgeMapWriter(const Map& _map, const Writer& _writer) 
	: map(_map), writer(_writer) {}

      virtual ~UEdgeMapWriter() {}

      virtual void write(std::ostream& os, const Item& item) const {
	Value value = map[item];
	writer.write(os, value);
      }

      virtual void sort(const _UGraph& ugraph, std::vector<Edge>& items) const {
        UEdgeComposeLess<_UGraph, Map> less(ugraph, map);
        std::sort(items.begin(), items.end(), less);
      }

      virtual void sort(std::vector<UEdge>& items) const {
        ComposeLess<Map> less(map);
        std::sort(items.begin(), items.end(), less);
      }

    };


    class ValueWriterBase {
    public:
      virtual ~ValueWriterBase() {}
      virtual void write(std::ostream&) = 0;
    };

    template <typename _Value, typename _Writer>
    class ValueWriter : public ValueWriterBase {
    public:
      typedef _Value Value;
      typedef _Writer Writer;

      ValueWriter(const Value& _value, const Writer& _writer)
 	: value(_value), writer(_writer) {}

      virtual void write(std::ostream& os) {
	writer.write(os, value);
      }
    private:
      const Value& value;
      Writer writer;
    };
    

    template <typename _Item>
    class LabelWriterBase {
    public:
      typedef _Item Item;
      virtual ~LabelWriterBase() {}
      virtual void write(std::ostream&, const Item&) const = 0;
      virtual void sort(std::vector<Item>&) const = 0;
      virtual bool isLabelWriter() const = 0;
      virtual LabelWriterBase* clone() const = 0;
    };

    template <typename _Item, typename _BoxedLabelWriter>
    class LabelWriter : public LabelWriterBase<_Item> {
    public:
      typedef _Item Item;
      typedef _BoxedLabelWriter BoxedLabelWriter;

      const BoxedLabelWriter& labelWriter;

      LabelWriter(const BoxedLabelWriter& _labelWriter) 
	: labelWriter(_labelWriter) {}

      virtual void write(std::ostream& os, const Item& item) const {
	labelWriter.writeLabel(os, item);
      }
      virtual void sort(std::vector<Item>& items) const {
	labelWriter.sortByLabel(items);
      }

      virtual bool isLabelWriter() const {
	return labelWriter.isLabelWriter();
      }

      virtual LabelWriter* clone() const {
	return new LabelWriter(labelWriter);
      }
    };

  }

  /// \ingroup lemon_io
  /// \brief Lemon Format writer class.
  /// 
  /// The Lemon Format contains several sections. We do not want to
  /// determine what sections are in a lemon file we give only a framework
  /// to write a section oriented format.
  ///
  /// In the Lemon Format each section starts with a line contains a \c \@
  /// character on the first not white space position. This line is the
  /// header line of the section. Each next lines belong to this section
  /// while it does not starts with \c \@ character. This line can start a 
  /// new section or if it can close the file with the \c \@end line.
  /// The file format ignore the empty lines and it may contain comments
  /// started with a \c # character to the end of the line. 
  ///
  /// The framework provides an abstract LemonWriter::SectionWriter class
  /// what defines the interface of a SectionWriter. The SectionWriter
  /// has the \c header() member function what gives back the header of the
  /// section. After that it will be called the \c write() member which
  /// should write the content of the section.
  ///
  /// \relates GraphWriter
  /// \relates NodeSetWriter
  /// \relates EdgeSetWriter
  /// \relates NodesWriter
  /// \relates EdgesWriter
  /// \relates AttributeWriter
  class LemonWriter {
  public:

    /// \brief Abstract base class for writing a section.
    ///
    /// This class has an \c header() member function what gives back
    /// the header line of the section. The \c write() member should
    /// write the content of the section to the stream.
    class SectionWriter {
      friend class LemonWriter;
    protected:
      /// \brief Constructor for SectionWriter.
      ///
      /// Constructor for SectionWriter. It attach this writer to
      /// the given LemonWriter.
      SectionWriter(LemonWriter& writer) {
	writer.attach(*this);
      }
      
      virtual ~SectionWriter() {}

      /// \brief The header of section.
      ///
      /// It gives back the header of the section.
      virtual std::string header() = 0;

      /// \brief Writer function of the section.
      ///
      /// Write the content of the section.
      virtual void write(std::ostream& os) = 0;
      
      /// \brief Gives back true when the section should be written.
      ///
      /// Gives back true when the section should be written.
      virtual bool valid() { return true; }
    };

    /// \brief Constructor for LemonWriter.
    ///
    /// Constructor for LemonWriter which writes to the given stream.
    LemonWriter(std::ostream& _os) 
      : os(&_os), own_os(false) {}

    /// \brief Constructor for LemonWriter.
    ///
    /// Constructor for LemonWriter which writes to the given file.
    LemonWriter(const std::string& filename) 
      : os(0), own_os(true) {
      os = new std::ofstream(filename.c_str());
    }

    /// \brief Desctructor for LemonWriter.
    ///
    /// Desctructor for LemonWriter.
    ~LemonWriter() {
      if (own_os) {
	delete os;
      }
    }

  private:
    LemonWriter(const LemonWriter&);
    void operator=(const LemonWriter&);

    void attach(SectionWriter& writer) {
      writers.push_back(&writer);
    }

  public:

    /// \brief Executes the LemonWriter.
    /// 
    /// It executes the LemonWriter.
    void run() {
      SectionWriters::iterator it;
      for (it = writers.begin(); it != writers.end(); ++it) {
        if ((*it)->valid()) {
          *os << (*it)->header() << std::endl;
          (*it)->write(*os);
        }
      }
      *os << "@end" << std::endl;
    }


  private:

    std::ostream* os;
    bool own_os;

    typedef std::vector<SectionWriter*> SectionWriters;
    SectionWriters writers;

  };

  /// \ingroup section_io
  /// \brief SectionWriter for writing a graph's nodeset.
  ///
  /// The lemon format can store multiple graph nodesets with several maps.
  /// The nodeset section's header line is \c \@nodeset \c nodeset_name, but 
  /// the \c nodeset_name may be empty.
  ///
  /// The first line of the section contains the names of the maps separated
  /// with white spaces. Each next lines describes a node in the nodeset, and
  /// contains the mapped values for each map.
  ///
  /// If the nodeset contains an \c "label" named map then it will be regarded
  /// as label map. This map should contain only unique values and when the 
  /// \c writeLabel() member will be called with a node it will write it's 
  /// label. Otherwise if the \c _forceLabelMap constructor parameter is true 
  /// then the label map will be the id in the graph. In addition if the
  /// the \c _forceSort is true then the writer will write the nodes
  /// sorted by the labels.
  ///
  /// \relates LemonWriter
  template <typename _Graph, typename _Traits = DefaultWriterTraits>
  class NodeSetWriter : public LemonWriter::SectionWriter {
    typedef LemonWriter::SectionWriter Parent;
  public:

    typedef _Graph Graph;
    typedef _Traits Traits;
    typedef typename Graph::Node Node;

    /// \brief Constructor.
    ///
    /// Constructor for NodeSetWriter. It creates the NodeSetWriter and
    /// attach it into the given LemonWriter. If the \c _forceLabelMap
    /// parameter is true then the writer will write own label map when
    /// the user does not give "label" named map. In addition if the
    /// the \c _forceSort is true then the writer will write the edges
    /// sorted by the labels.
    NodeSetWriter(LemonWriter& _writer, const Graph& _graph, 
		  const std::string& _name = std::string(), 
		  bool _forceLabelMap = true, bool _forceSort = true) 
      : Parent(_writer), labelMap(0), forceLabelMap(_forceLabelMap), 
	forceSort(_forceSort), graph(_graph), name(_name) {}

    /// \brief Destructor.
    ///
    /// Destructor for NodeSetWriter.
    virtual ~NodeSetWriter() {
      typename MapWriters::iterator it;
      for (it = writers.begin(); it != writers.end(); ++it) {
	delete it->second;
      }
    }

  private:
    NodeSetWriter(const NodeSetWriter&);
    void operator=(const NodeSetWriter&);
  
  public:

    /// \brief Add a new node map writer command for the writer.
    ///
    /// Add a new node map writer command for the writer.
    template <typename Map>
    NodeSetWriter& writeNodeMap(std::string label, const Map& map) {
      return writeNodeMap<typename Traits::
	template Writer<typename Map::Value>, Map>(label, map);
    }

    /// \brief Add a new node map writer command for the writer.
    ///
    /// Add a new node map writer command for the writer.
    template <typename ItemWriter, typename Map>
    NodeSetWriter& writeNodeMap(std::string label, const Map& map, 
			    const ItemWriter& iw = ItemWriter()) {
      checkConcept<concepts::ReadMap<Node, typename Map::Value>, Map>();
      checkConcept<_writer_bits::ItemWriter<typename Map::Value>,ItemWriter>();
      writers.push_back(
	make_pair(label, new _writer_bits::
		  MapWriter<Node, Map, ItemWriter>(map, iw)));
      return *this;
    }

  protected:

    /// \brief The header of the section.
    ///
    /// It gives back the header of the section.
    virtual std::string header() {
      return "@nodeset " + name;
    }

    /// \brief  Writer function of the section.
    ///
    /// Write the content of the section.
    virtual void write(std::ostream& os) {
      for (int i = 0; i < int(writers.size()); ++i) {
	if (writers[i].first == "label") {
	  labelMap = writers[i].second;
	  forceLabelMap = false;
	  break;
	}
      }
      std::vector<Node> items;
      for (typename Graph::NodeIt it(graph); it != INVALID; ++it) {
        items.push_back(it);
      }
      if (forceSort) {
        if (labelMap) {
          labelMap->sort(items);
        } else {
          typedef IdMap<Graph, Node> Map;
          Map map(graph);
          _writer_bits::ComposeLess<Map> less(map);
          std::sort(items.begin(), items.end(), less);
        }
      }
      if (forceLabelMap) {
	os << "label\t";
      }
      for (int i = 0; i < int(writers.size()); ++i) {
	os << writers[i].first << '\t';
      }
      os << std::endl;
      for (typename std::vector<Node>::iterator it = items.begin();
           it != items.end(); ++it) {
	if (forceLabelMap) {
	  os << graph.id(*it) << '\t';
	}
	for (int i = 0; i < int(writers.size()); ++i) {
	  writers[i].second->write(os, *it);
	  os << '\t';
	}
	os << std::endl;
      }
    }

  public:

    /// \brief Returns true if the nodeset can write the labels of the nodes.
    ///
    /// Returns true if the nodeset can write the labels of the nodes.
    /// It is possible only if a "label" named map was written or the 
    /// \c _forceLabelMap constructor parameter was true.
    bool isLabelWriter() const {
      return labelMap != 0 || forceLabelMap;
    }

    /// \brief Write the label of the given node.
    ///
    /// It writes the label of the given node. If there was written a "label"
    /// named map then it will write the map value belongs to the node.
    /// Otherwise if the \c forceLabel parameter was true it will write
    /// its label in the graph. 
    void writeLabel(std::ostream& os, const Node& item) const {
      if (forceLabelMap) {
	os << graph.id(item);
      } else {
	labelMap->write(os, item);
      }
    }

    /// \brief Sorts the given node vector by label.
    ///
    /// Sorts the given node vector by label. If there was written an
    /// "label" named map then the vector will be sorted by the values
    /// of this map. Otherwise if the \c forceLabel parameter was true
    /// it will be sorted by its id in the graph.
    void sortByLabel(std::vector<Node>& nodes) const {
      if (labelMap) {
	labelMap->sort(nodes);
      } else {
	typedef IdMap<Graph, Node> Map;
	Map map(graph);
	_writer_bits::ComposeLess<Map> less(map);
	std::sort(nodes.begin(), nodes.end(), less);
      }
    }

  private:

    typedef std::vector<std::pair<std::string, _writer_bits::
				  MapWriterBase<Node>*> > MapWriters;
    MapWriters writers;

    _writer_bits::MapWriterBase<Node>* labelMap;
    bool forceLabelMap;
    bool forceSort;
   
    const Graph& graph;   
    std::string name;

  };

  /// \ingroup section_io
  /// \brief SectionWriter for writing a bipartite graph's nodeset.
  ///
  /// The lemon format can store multiple bipartite graph nodesets
  /// with several maps.  The nodeset section's header line is \c
  /// \@bpnodeset \c bpnodeset_name, but the \c bpnodeset_name may be empty.
  ///
  /// The first line of the section contains the names of the maps separated
  /// with white spaces. Each next lines describes a node in the nodeset, and
  /// contains the mapped values for each map.
  ///
  /// If the nodeset contains an \c "label" named map then it will be regarded
  /// as label map. This map should contain only unique values and when the 
  /// \c writeLabel() member will be called with a node it will write it's 
  /// label. Otherwise if the \c _forceLabelMap constructor parameter is true 
  /// then the label map will be the id in the graph. In addition if the
  /// the \c _forceSort is true then the writer will write the edges
  /// sorted by the labels.
  ///
  /// \relates LemonWriter
  template <typename _Graph, typename _Traits = DefaultWriterTraits>
  class BpNodeSetWriter : public LemonWriter::SectionWriter {
    typedef LemonWriter::SectionWriter Parent;
  public:

    typedef _Graph Graph;
    typedef _Traits Traits;
    typedef typename Graph::Node Node;

    /// \brief Constructor.
    ///
    /// Constructor for BpNodeSetWriter. It creates the BpNodeSetWriter and
    /// attach it into the given LemonWriter. If the \c _forceLabelMap
    /// parameter is true then the writer will write own label map when
    /// the user does not give "label" named map. In addition if the
    /// the \c _forceSort is true then the writer will write the nodes
    /// sorted by the labels.
    BpNodeSetWriter(LemonWriter& _writer, const Graph& _graph, 
		  const std::string& _name = std::string(), 
		  bool _forceLabelMap = true, bool _forceSort = true) 
      : Parent(_writer), labelMap(0), forceLabelMap(_forceLabelMap), 
	forceSort(_forceSort), graph(_graph), name(_name) {}

    /// \brief Destructor.
    ///
    /// Destructor for BpNodeSetWriter.
    virtual ~BpNodeSetWriter() {
      typename MapWriters::iterator it;
      for (it = writers.begin(); it != writers.end(); ++it) {
	delete it->second;
      }
    }

  private:
    BpNodeSetWriter(const BpNodeSetWriter&);
    void operator=(const BpNodeSetWriter&);
  
  public:

    /// \brief Add a new A-node map writer command for the writer.
    ///
    /// Add a new A-node map writer command for the writer.
    template <typename Map>
    BpNodeSetWriter& writeANodeMap(std::string label, const Map& map) {
      return writeANodeMap<typename Traits::
	template Writer<typename Map::Value>, Map>(label, map);
    }

    /// \brief Add a new A-node map writer command for the writer.
    ///
    /// Add a new A-node map writer command for the writer.
    template <typename ItemWriter, typename Map>
    BpNodeSetWriter& writeANodeMap(std::string label, const Map& map, 
				   const ItemWriter& iw = ItemWriter()) {
      checkConcept<concepts::ReadMap<Node, typename Map::Value>, Map>();
      checkConcept<_writer_bits::ItemWriter<typename Map::Value>,ItemWriter>();
      if (label == "label") {
	throw IoParameterError("Label cannot be A-node map");
      }
      awriters.push_back(make_pair(label, new _writer_bits::
				   MapWriter<Node, Map, ItemWriter>(map, iw)));
      return *this;
    }

    /// \brief Add a new B-node map writer command for the writer.
    ///
    /// Add a new B-node map writer command for the writer.
    template <typename Map>
    BpNodeSetWriter& writeBNodeMap(std::string label, const Map& map) {
      return writeBNodeMap<typename Traits::
	template Writer<typename Map::Value>, Map>(label, map);
    }

    /// \brief Add a new B-node map writer command for the writer.
    ///
    /// Add a new B-node map writer command for the writer.
    template <typename ItemWriter, typename Map>
    BpNodeSetWriter& writeBNodeMap(std::string label, const Map& map, 
				   const ItemWriter& iw = ItemWriter()) {
      checkConcept<concepts::ReadMap<Node, typename Map::Value>, Map>();
      checkConcept<_writer_bits::ItemWriter<typename Map::Value>,ItemWriter>();
      if (label == "label") {
	throw IoParameterError("Label cannot be B-node map");
      }
      bwriters.push_back(make_pair(label, new _writer_bits::
				   MapWriter<Node, Map, ItemWriter>(map, iw)));
      return *this;
    }

    /// \brief Add a new node map writer command for the writer.
    ///
    /// Add a new node map writer command for the writer.
    template <typename Map>
    BpNodeSetWriter& writeNodeMap(std::string label, const Map& map) {
      return writeNodeMap<typename Traits::
	template Writer<typename Map::Value>, Map>(label, map);
    }

    /// \brief Add a new node map writer command for the writer.
    ///
    /// Add a new node map writer command for the writer.
    template <typename ItemWriter, typename Map>
    BpNodeSetWriter& writeNodeMap(std::string label, const Map& map, 
				  const ItemWriter& iw = ItemWriter()) {
      checkConcept<concepts::ReadMap<Node, typename Map::Value>, Map>();
      checkConcept<_writer_bits::ItemWriter<typename Map::Value>,ItemWriter>();
      writers.push_back(make_pair(label, new _writer_bits::
				  MapWriter<Node, Map, ItemWriter>(map, iw)));
      return *this;
    }

  protected:

    /// \brief The header of the section.
    ///
    /// It gives back the header of the section.
    virtual std::string header() {
      return "@bpnodeset " + name;
    }

    /// \brief Writer function of the section.
    ///
    /// Write the content of the section.
    virtual void write(std::ostream& os) {
      for (int i = 0; i < int(writers.size()); ++i) {
	if (writers[i].first == "label") {
	  labelMap = writers[i].second;
	  forceLabelMap = false;
	  break;
	}
      }
      {
	os << "&anodeset ";
	std::vector<Node> items;
	for (typename Graph::ANodeIt it(graph); it != INVALID; ++it) {
	  items.push_back(it);
	}
	if (forceSort) {
	  if (labelMap) {
	    labelMap->sort(items);
	  } else {
	    typedef IdMap<Graph, Node> Map;
	    Map map(graph);
	    _writer_bits::ComposeLess<Map> less(map);
	    std::sort(items.begin(), items.end(), less);
	  }
	}
	if (forceLabelMap) {
	  os << "label\t";
	}
	for (int i = 0; i < int(writers.size()); ++i) {
	  os << writers[i].first << '\t';
	}
	for (int i = 0; i < int(awriters.size()); ++i) {
	  os << awriters[i].first << '\t';
	}
	os << std::endl;
	for (typename std::vector<Node>::iterator it = items.begin();
	     it != items.end(); ++it) {
	  if (forceLabelMap) {
	    os << graph.id(*it) << '\t';
	  }
	  for (int i = 0; i < int(writers.size()); ++i) {
	    writers[i].second->write(os, *it);
	    os << '\t';
	  }
	  for (int i = 0; i < int(awriters.size()); ++i) {
	    awriters[i].second->write(os, *it);
	    os << '\t';
	  }
	  os << std::endl;
	}
      }
      {
	os << "&bnodeset ";
	std::vector<Node> items;
	for (typename Graph::BNodeIt it(graph); it != INVALID; ++it) {
	  items.push_back(it);
	}
	if (forceSort) {
	  if (labelMap) {
	    labelMap->sort(items);
	  } else {
	    typedef IdMap<Graph, Node> Map;
	    Map map(graph);
	    _writer_bits::ComposeLess<Map> less(map);
	    std::sort(items.begin(), items.end(), less);
	  }
	}
	if (forceLabelMap) {
	  os << "label\t";
	}
	for (int i = 0; i < int(writers.size()); ++i) {
	  os << writers[i].first << '\t';
	}
	for (int i = 0; i < int(bwriters.size()); ++i) {
	  os << bwriters[i].first << '\t';
	}
	os << std::endl;
	for (typename std::vector<Node>::iterator it = items.begin();
	     it != items.end(); ++it) {
	  if (forceLabelMap) {
	    os << graph.id(*it) << '\t';
	  }
	  for (int i = 0; i < int(writers.size()); ++i) {
	    writers[i].second->write(os, *it);
	    os << '\t';
	  }
	  for (int i = 0; i < int(bwriters.size()); ++i) {
	    bwriters[i].second->write(os, *it);
	    os << '\t';
	  }
	  os << std::endl;
	}
      }
    }

  public:

    /// \brief Returns true if the nodeset can write the labels of the nodes.
    ///
    /// Returns true if the nodeset can write the labels of the nodes.
    /// It is possible only if a "label" named map was written or the 
    /// \c _forceLabelMap constructor parameter was true.
    bool isLabelWriter() const {
      return labelMap != 0 || forceLabelMap;
    }

    /// \brief Write the label of the given node.
    ///
    /// It writes the label of the given node. If there was written a "label"
    /// named map then it will write the map value belongs to the node.
    /// Otherwise if the \c forceLabel parameter was true it will write
    /// its label in the graph. 
    void writeLabel(std::ostream& os, const Node& item) const {
      if (forceLabelMap) {
	os << graph.id(item);
      } else {
	labelMap->write(os, item);
      }
    }

    /// \brief Sorts the given node vector by label.
    ///
    /// Sorts the given node vector by label. If there was written an
    /// "label" named map then the vector will be sorted by the values
    /// of this map. Otherwise if the \c forceLabel parameter was true
    /// it will be sorted by its id in the graph.
    void sortByLabel(std::vector<Node>& nodes) const {
      if (labelMap) {
	labelMap->sort(nodes);
      } else {
	typedef IdMap<Graph, Node> Map;
	Map map(graph);
	_writer_bits::ComposeLess<Map> less(map);
	std::sort(nodes.begin(), nodes.end(), less);
      }
    }

  private:

    typedef std::vector<std::pair<std::string, _writer_bits::
				  MapWriterBase<Node>*> > MapWriters;
    MapWriters awriters, bwriters, writers;
    
    _writer_bits::MapWriterBase<Node>* labelMap;
    bool forceLabelMap;
    bool forceSort;
   
    const Graph& graph;   
    std::string name;

  };

  /// \ingroup section_io
  /// \brief SectionWriter for writing a graph's edgesets.
  ///
  /// The lemon format can store multiple graph edgesets with several maps. 
  /// The edgeset section's header line is \c \@edgeset \c edgeset_name, but 
  /// the \c edgeset_name may be empty.
  ///
  /// The first line of the section contains the names of the maps separated
  /// with white spaces. Each next lines describes a edge in the edgeset. The
  /// line contains the source and the target nodes' label and the mapped 
  /// values for each map.
  ///
  /// If the edgeset contains an \c "label" named map then it will be regarded
  /// as label map. This map should contain only unique values and when the 
  /// \c writeLabel() member will be called with an edge it will write it's 
  /// label. Otherwise if the \c _forceLabelMap constructor parameter is true 
  /// then the label map will be the id in the graph. In addition if the
  /// the \c _forceSort is true then the writer will write the edges
  /// sorted by the labels.
  ///
  /// The edgeset writer needs a node label writer to identify which nodes
  /// have to be connected. If a NodeSetWriter can write the nodes' label,
  /// it will be able to use with this class.
  ///
  /// \relates LemonWriter
  template <typename _Graph, typename _Traits = DefaultWriterTraits>
  class EdgeSetWriter : public LemonWriter::SectionWriter {
    typedef LemonWriter::SectionWriter Parent;
  public:

    typedef _Graph Graph;
    typedef _Traits Traits;
    typedef typename Graph::Node Node;
    typedef typename Graph::Edge Edge;

    /// \brief Constructor.
    ///
    /// Constructor for EdgeSetWriter. It creates the EdgeSetWriter
    /// and attach it into the given LemonWriter. It will write node
    /// labels by the \c _nodeLabelWriter. If the \c _forceLabelMap
    /// parameter is true then the writer will write own label map if
    /// the user does not give "label" named map. In addition if the
    /// the \c _forceSort is true then the writer will write the
    /// edges sorted by the labels.
    template <typename NodeLabelWriter>
    EdgeSetWriter(LemonWriter& _writer, const Graph& _graph, 
		  const NodeLabelWriter& _nodeLabelWriter, 
		  const std::string& _name = std::string(),
		  bool _forceLabelMap = true, bool _forceSort = true)
      : Parent(_writer), labelMap(0), forceLabelMap(_forceLabelMap),
	forceSort(_forceSort), graph(_graph), name(_name) {
      checkConcept<_writer_bits::ItemLabelWriter<Node>, NodeLabelWriter>();
      nodeLabelWriter.reset(new _writer_bits::
			 LabelWriter<Node, NodeLabelWriter>(_nodeLabelWriter));
    } 

    /// \brief Destructor.
    ///
    /// Destructor for EdgeSetWriter.
    virtual ~EdgeSetWriter() {
      typename MapWriters::iterator it;
      for (it = writers.begin(); it != writers.end(); ++it) {
	delete it->second;
      }
    }

  private:
    EdgeSetWriter(const EdgeSetWriter&);
    void operator=(const EdgeSetWriter&);

  public:

    /// \brief Add a new edge map writer command for the writer.
    ///
    /// Add a new edge map writer command for the writer.
    template <typename Map>
    EdgeSetWriter& writeEdgeMap(std::string label, const Map& map) {
      return writeEdgeMap<typename Traits::
	template Writer<typename Map::Value>, Map>(label, map);
    }

    /// \brief Add a new edge map writer command for the writer.
    ///
    /// Add a new edge map writer command for the writer.
    template <typename ItemWriter, typename Map>
    EdgeSetWriter& writeEdgeMap(std::string label, const Map& map, 
			    const ItemWriter& iw = ItemWriter()) {
      checkConcept<concepts::ReadMap<Edge, typename Map::Value>, Map>();
      checkConcept<_writer_bits::ItemWriter<typename Map::Value>, ItemWriter>();
      writers.push_back(
	make_pair(label, new _writer_bits::
		  MapWriter<Edge, Map, ItemWriter>(map, iw)));
      return *this;
    }

  protected:

    /// \brief The header of the section.
    ///
    /// It gives back the header of the section.
    virtual std::string header() {
      return "@edgeset " + name;
    }

    /// \brief  Writer function of the section.
    ///
    /// Write the content of the section.
    virtual void write(std::ostream& os) {
      if (!nodeLabelWriter->isLabelWriter()) {
	throw DataFormatError("Cannot find nodeset or label map");
      }
      for (int i = 0; i < int(writers.size()); ++i) {
	if (writers[i].first == "label") {
	  labelMap = writers[i].second;
	  forceLabelMap = false;
	  break;
	}
      }
      std::vector<Edge> items;
      for (typename Graph::EdgeIt it(graph); it != INVALID; ++it) {
        items.push_back(it);
      }
      if (forceSort) {
        if (labelMap) {
          labelMap->sort(items);
        } else {
          typedef IdMap<Graph, Edge> Map;
          Map map(graph);
          _writer_bits::ComposeLess<Map> less(map);
          std::sort(items.begin(), items.end(), less);
        }
      }
      os << "\t\t";
      if (forceLabelMap) {
	os << "label\t";
      }
      for (int i = 0; i < int(writers.size()); ++i) {
	os << writers[i].first << '\t';
      }
      os << std::endl;
      for (typename std::vector<Edge>::iterator it = items.begin();
           it != items.end(); ++it) {
	nodeLabelWriter->write(os, graph.source(*it));
	os << '\t';
	nodeLabelWriter->write(os, graph.target(*it));
	os << '\t';
	if (forceLabelMap) {
	  os << graph.id(*it) << '\t';
	}
	for (int i = 0; i < int(writers.size()); ++i) {
	  writers[i].second->write(os, *it);
	  os << '\t';
	}
	os << std::endl;
      }
    }

  public:

    /// \brief Returns true if the edgeset can write the labels of the edges.
    ///
    /// Returns true if the edgeset can write the labels of the edges.
    /// It is possible only if a "label" named map was written or the 
    /// \c _forceLabelMap constructor parameter was true.
    bool isLabelWriter() const {
      return forceLabelMap || labelMap != 0;
    }

    /// \brief Write the label of the given edge.
    ///
    /// It writes the label of the given edge. If there was written a "label"
    /// named map then it will write the map value belongs to the edge.
    /// Otherwise if the \c forceLabel parameter was true it will write
    /// its label in the graph. 
    void writeLabel(std::ostream& os, const Edge& item) const {
      if (forceLabelMap) {
	os << graph.id(item);
      } else {
	labelMap->write(os, item);
      }
    } 

    /// \brief Sorts the given edge vector by label.
    ///
    /// Sorts the given edge vector by label. If there was written an
    /// "label" named map then the vector will be sorted by the values
    /// of this map. Otherwise if the \c forceLabel parameter was true
    /// it will be sorted by its id in the graph.
    void sortByLabel(std::vector<Edge>& edges) const {
      if (labelMap) {
	labelMap->sort(edges);
      } else {
	typedef IdMap<Graph, Edge> Map;
	Map map(graph);
	_writer_bits::ComposeLess<Map> less(map);
	std::sort(edges.begin(), edges.end(), less);
      }
    }

  private:

    typedef std::vector<std::pair<std::string, _writer_bits::
				  MapWriterBase<Edge>*> > MapWriters;
    MapWriters writers;

    _writer_bits::MapWriterBase<Edge>* labelMap;
    bool forceLabelMap;
    bool forceSort;
   
    const Graph& graph;   
    std::string name;

    std::auto_ptr<_writer_bits::LabelWriterBase<Node> > nodeLabelWriter;
  };

  /// \ingroup section_io
  /// \brief SectionWriter for writing a undirected edgeset.
  ///
  /// The lemon format can store multiple undirected edgesets with several 
  /// maps. The undirected edgeset section's header line is \c \@uedgeset 
  /// \c uedgeset_name, but the \c uedgeset_name may be empty.
  ///
  /// The first line of the section contains the names of the maps separated
  /// with white spaces. Each next lines describes an undirected edge in the 
  /// edgeset. The line contains the two connected nodes' label and the mapped 
  /// values for each undirected map.
  ///
  /// The section can handle the directed as a syntactical sugar. Two
  /// undirected edge map describes one directed edge map. This two maps
  /// are the forward map and the backward map and the names of this map
  /// is near the same just with a prefix \c '+' or \c '-' character 
  /// difference.
  ///
  /// If the edgeset contains an \c "label" named map then it will be
  /// regarded as label map. This map should contain only unique
  /// values and when the \c writeLabel() member will be called with
  /// an undirected edge it will write it's label. Otherwise if the \c
  /// _forceLabelMap constructor parameter is true then the label map
  /// will be the id in the graph.  In addition if the the \c
  /// _forceSort is true then the writer will write the edges sorted
  /// by the labels.
  ///
  /// The undirected edgeset writer needs a node label writer to identify 
  /// which nodes have to be connected. If a NodeSetWriter can write the 
  /// nodes' label, it will be able to use with this class.
  ///
  /// \relates LemonWriter
  template <typename _Graph, typename _Traits = DefaultWriterTraits>
  class UEdgeSetWriter : public LemonWriter::SectionWriter {
    typedef LemonWriter::SectionWriter Parent;
  public:

    typedef _Graph Graph;
    typedef _Traits Traits;
    typedef typename Graph::Node Node;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::UEdge UEdge;

    /// \brief Constructor.
    ///
    /// Constructor for UEdgeSetWriter. It creates the UEdgeSetWriter
    /// and attach it into the given LemonWriter. It will write node
    /// labels by the \c _nodeLabelWriter. If the \c _forceLabelMap
    /// parameter is true then the writer will write own label map if
    /// the user does not give "label" named map. In addition if the
    /// the \c _forceSort is true then the writer will write the
    /// edges sorted by the labels.
    template <typename NodeLabelWriter>
    UEdgeSetWriter(LemonWriter& _writer, const Graph& _graph, 
		       const NodeLabelWriter& _nodeLabelWriter, 
		       const std::string& _name = std::string(),
		       bool _forceLabelMap = true, bool _forceSort = true)
      : Parent(_writer), labelMap(0), forceLabelMap(_forceLabelMap),
	forceSort(_forceSort), graph(_graph), name(_name) {
      checkConcept<_writer_bits::ItemLabelWriter<Node>, NodeLabelWriter>();
      nodeLabelWriter.reset(new _writer_bits::
	LabelWriter<Node, NodeLabelWriter>(_nodeLabelWriter));
    } 

    /// \brief Destructor.
    ///
    /// Destructor for UEdgeSetWriter.
    virtual ~UEdgeSetWriter() {
      typename MapWriters::iterator it;
      for (it = writers.begin(); it != writers.end(); ++it) {
	delete it->second;
      }
    }

  private:
    UEdgeSetWriter(const UEdgeSetWriter&);
    void operator=(const UEdgeSetWriter&);

  public:

    /// \brief Add a new undirected edge map writer command for the writer.
    ///
    /// Add a new undirected map writer command for the writer.
    template <typename Map>
    UEdgeSetWriter& writeUEdgeMap(std::string label, const Map& map) {
      return writeUEdgeMap<typename Traits::
	template Writer<typename Map::Value>, Map>(label, map);
    }

    /// \brief Add a new undirected map writer command for the writer.
    ///
    /// Add a new undirected map writer command for the writer.
    template <typename ItemWriter, typename Map>
    UEdgeSetWriter& writeUEdgeMap(std::string label, const Map& map, 
                                  const ItemWriter& iw = ItemWriter()) {
      checkConcept<concepts::ReadMap<UEdge, typename Map::Value>, Map>();
      checkConcept<_writer_bits::ItemWriter<typename Map::Value>, ItemWriter>();
      writers.push_back(
	make_pair(label, new _writer_bits::
		  UEdgeMapWriter<Graph, Map, ItemWriter>(map, iw)));
      return *this;
    }

    /// \brief Add a new directed edge map writer command for the writer.
    ///
    /// Add a new directed map writer command for the writer.
    template <typename Map>
    UEdgeSetWriter& writeEdgeMap(std::string label, const Map& map) {
      return writeEdgeMap<typename Traits::
	template Writer<typename Map::Value>, Map>(label, map);
    }

    /// \brief Add a new directed map writer command for the writer.
    ///
    /// Add a new directed map writer command for the writer.
    template <typename ItemWriter, typename Map>
    UEdgeSetWriter& writeEdgeMap(std::string label, const Map& map, 
                                 const ItemWriter& iw = ItemWriter()) {
      checkConcept<concepts::ReadMap<Edge, typename Map::Value>, Map>();
      checkConcept<_writer_bits::ItemWriter<typename Map::Value>, ItemWriter>();
      writeUEdgeMap("+" + label, 
                    _writer_bits::forwardComposeMap(graph, map), iw);
      writeUEdgeMap("-" + label, 
                    _writer_bits::backwardComposeMap(graph, map), iw);
      return *this;
    }

  protected:

    /// \brief The header of the section.
    ///
    /// It gives back the header of the section.
    virtual std::string header() {
      return "@uedgeset " + name;
    }

    /// \brief  Writer function of the section.
    ///
    /// Write the content of the section.
    virtual void write(std::ostream& os) {
      if (!nodeLabelWriter->isLabelWriter()) {
	throw DataFormatError("Cannot find nodeset or label map");
      }
      for (int i = 0; i < int(writers.size()); ++i) {
	if (writers[i].first == "label") {
	  labelMap = writers[i].second;
	  forceLabelMap = false;
	  break;
	}
      }
      std::vector<UEdge> items;
      for (typename Graph::UEdgeIt it(graph); it != INVALID; ++it) {
        items.push_back(it);
      }
      if (forceSort) {
        if (labelMap) {
          labelMap->sort(items);
        } else {
          typedef IdMap<Graph, UEdge> Map;
          Map map(graph);
          _writer_bits::ComposeLess<Map> less(map);
          std::sort(items.begin(), items.end(), less);
        }
      }
      os << "\t\t";
      if (forceLabelMap) {
	os << "label\t";
      }
      for (int i = 0; i < int(writers.size()); ++i) {
	os << writers[i].first << '\t';
      }
      os << std::endl;
      for (typename std::vector<UEdge>::iterator it = items.begin();
           it != items.end(); ++it) {
	nodeLabelWriter->write(os, graph.source(*it));
	os << '\t';
	nodeLabelWriter->write(os, graph.target(*it));
	os << '\t';
	if (forceLabelMap) {
	  os << graph.id(*it) << '\t';
	}
	for (int i = 0; i < int(writers.size()); ++i) {
	  writers[i].second->write(os, *it);
	  os << '\t';
	}
	os << std::endl;
      }
    }

  public:

    /// \brief Returns true if the undirected edgeset can write the labels of 
    /// the edges.
    ///
    /// Returns true if the undirected edgeset can write the labels of the 
    /// undirected edges. It is possible only if a "label" named map was 
    /// written or the \c _forceLabelMap constructor parameter was true.
    bool isLabelWriter() const {
      return forceLabelMap || labelMap != 0;
    }

    /// \brief Write the label of the given undirected edge.
    ///
    /// It writes the label of the given undirected edge. If there was written 
    /// a "label" named map then it will write the map value belongs to the 
    /// undirected edge. Otherwise if the \c forceLabel parameter was true it 
    /// will write its id in the graph. 
    void writeLabel(std::ostream& os, const UEdge& item) const {
      if (forceLabelMap) {
	os << graph.id(item);
      } else {
	labelMap->write(os, item);
      }
    } 

    /// \brief Write the label of the given edge.
    ///
    /// It writes the label of the given edge. If there was written 
    /// a "label" named map then it will write the map value belongs to the 
    /// edge. Otherwise if the \c forceLabel parameter was true it 
    /// will write its id in the graph. If the edge is forward map
    /// then its prefix character is \c '+' elsewhere \c '-'.
    void writeLabel(std::ostream& os, const Edge& item) const {
      if (graph.direction(item)) {
	os << "+";
      } else {
	os << "-";
      }
      if (forceLabelMap) {
	os << graph.id(static_cast<const UEdge&>(item));
      } else {
	labelMap->write(os, item);
      }
    } 

    /// \brief Sorts the given undirected edge vector by label.
    ///
    /// Sorts the given undirected edge vector by label. If there was
    /// written a "label" named map then the vector will be sorted by
    /// the values of this map. Otherwise if the \c forceLabel
    /// parameter was true it will be sorted by its id in the graph.
    void sortByLabel(std::vector<UEdge>& uedges) const {
      if (labelMap) {
	labelMap->sort(uedges);
      } else {
	typedef IdMap<Graph, UEdge> Map;
	Map map(graph);
	_writer_bits::ComposeLess<Map> less(map);
	std::sort(uedges.begin(), uedges.end(), less);
      }
    }

    /// \brief Sorts the given edge vector by label.
    ///
    /// Sorts the given edge vector by label. If there was written a
    /// "label" named map then the vector will be sorted by the values
    /// of this map. Otherwise if the \c forceLabel parameter was true
    /// it will be sorted by its id in the graph.
    void sortByLabel(std::vector<Edge>& edges) const {
      if (labelMap) {
	labelMap->sort(graph, edges);
      } else {
	typedef IdMap<Graph, Edge> Map;
	Map map(graph);
	_writer_bits::ComposeLess<Map> less(map);
	std::sort(edges.begin(), edges.end(), less);
      }
    }

  private:

    typedef std::vector<std::pair<std::string, _writer_bits::
				  UEdgeMapWriterBase<Graph>*> > MapWriters;
    MapWriters writers;

    _writer_bits::UEdgeMapWriterBase<Graph>* labelMap;
    bool forceLabelMap;
    bool forceSort;
   
    const Graph& graph;   
    std::string name;

    std::auto_ptr<_writer_bits::LabelWriterBase<Node> > nodeLabelWriter;
  };

  /// \ingroup section_io
  /// \brief SectionWriter for writing named nodes.
  ///
  /// The nodes section's header line is \c \@nodes \c nodes_name, but the
  /// \c nodes_name may be empty.
  ///
  /// Each line in the section contains the name of the node and 
  /// then the node label. 
  ///
  /// \relates LemonWriter
  template <typename _Graph>
  class NodeWriter : public LemonWriter::SectionWriter {
    typedef LemonWriter::SectionWriter Parent;
    typedef _Graph Graph;
    typedef typename Graph::Node Node;
  public:
    
    /// \brief Constructor.
    ///
    /// Constructor for NodeWriter. It creates the NodeWriter and
    /// attach it into the given LemonWriter. The given \c _LabelWriter
    /// will write the nodes' label what can be a nodeset writer.
    template <typename _LabelWriter>
    NodeWriter(LemonWriter& _writer, const _LabelWriter& _labelWriter, 
	       const std::string& _name = std::string()) 
      : Parent(_writer), name(_name) {
      checkConcept<_writer_bits::ItemLabelWriter<Node>, _LabelWriter>();
      labelWriter.reset(new _writer_bits::LabelWriter<Node, _LabelWriter>
                        (_labelWriter));
    }


    /// \brief Destructor.
    ///
    /// Destructor for NodeWriter.
    virtual ~NodeWriter() {}

  private:
    NodeWriter(const NodeWriter&);
    void operator=(const NodeWriter&);

  public:

    /// \brief Add a node writer command for the NodeWriter.
    ///
    /// Add a node writer command for the NodeWriter.
    void writeNode(std::string label, const Node& item) {
      writers.push_back(make_pair(label, &item));
    }

  protected:

    /// \brief The header of the section.
    ///
    /// It gives back the header of the section.
    virtual std::string header() {
      return "@nodes " + name;
    }

    /// \brief  Writer function of the section.
    ///
    /// Write the content of the section.
    virtual void write(std::ostream& os) {
      if (!labelWriter->isLabelWriter()) {
	throw DataFormatError("Cannot find nodeset or label map");
      }
      for (int i = 0; i < int(writers.size()); ++i) {
	os << writers[i].first << ' ';
	labelWriter->write(os, *(writers[i].second));
	os << std::endl;
      }
    }

    /// \brief Gives back true when the section should be written.
    ///
    /// Gives back true when the section should be written.
    virtual bool valid() { return !writers.empty(); }
    
  private:

    std::string name;

    typedef std::vector<std::pair<std::string, const Node*> > NodeWriters;
    NodeWriters writers;
    std::auto_ptr<_writer_bits::LabelWriterBase<Node> > labelWriter;
  };

  /// \ingroup section_io
  /// \brief SectionWriter for writing named edges.
  ///
  /// The edges section's header line is \c \@edges \c edges_name, but the
  /// \c edges_name may be empty.
  ///
  /// Each line in the section contains the name of the edge and 
  /// then the edge label. 
  ///
  /// \relates LemonWriter
  template <typename _Graph>
  class EdgeWriter : public LemonWriter::SectionWriter {
    typedef LemonWriter::SectionWriter Parent;
    typedef _Graph Graph;
    typedef typename Graph::Edge Edge;
  public:
    
    /// \brief Constructor.
    ///
    /// Constructor for EdgeWriter. It creates the EdgeWriter and
    /// attach it into the given LemonWriter. The given \c _LabelWriter
    /// will write the edges' label what can be a edgeset writer.
    template <typename _LabelWriter>
    EdgeWriter(LemonWriter& _writer, const _LabelWriter& _labelWriter, 
	       const std::string& _name = std::string()) 
      : Parent(_writer), name(_name) {
      checkConcept<_writer_bits::ItemLabelWriter<Edge>, _LabelWriter>();
      labelWriter.reset(new _writer_bits::LabelWriter<Edge, _LabelWriter>(_labelWriter));
    }

    /// \brief Destructor.
    ///
    /// Destructor for EdgeWriter.
    virtual ~EdgeWriter() {}
  private:
    EdgeWriter(const EdgeWriter&);
    void operator=(const EdgeWriter&);

  public:

    /// \brief Add an edge writer command for the EdgeWriter.
    ///
    /// Add an edge writer command for the EdgeWriter.
    void writeEdge(std::string label, const Edge& item) {
      writers.push_back(make_pair(label, &item));
    }

  protected:

    /// \brief The header of the section.
    ///
    /// It gives back the header of the section.
    virtual std::string header() {
      return "@edges " + name;
    }

    /// \brief  Writer function of the section.
    ///
    /// Write the content of the section.
    virtual void write(std::ostream& os) {
      if (!labelWriter->isLabelWriter()) {
	throw DataFormatError("Cannot find edgeset or label map");
      }
      for (int i = 0; i < int(writers.size()); ++i) {
	os << writers[i].first << ' ';
	labelWriter->write(os, *(writers[i].second));
	os << std::endl;
      }
    }

    /// \brief Gives back true when the section should be written.
    ///
    /// Gives back true when the section should be written.
    virtual bool valid() { return !writers.empty(); }
    
  private:

    std::string name;

    typedef std::vector<std::pair<std::string, const Edge*> > EdgeWriters;
    EdgeWriters writers;

    std::auto_ptr<_writer_bits::LabelWriterBase<Edge> > labelWriter;
  };


  /// \ingroup section_io
  /// \brief SectionWriter for writing named undirected edges.
  ///
  /// The undirected edges section's header line is \c \@uedges 
  /// \c uedges_name, but the \c uedges_name may be empty.
  ///
  /// Each line in the section contains the name of the undirected edge and 
  /// then the undirected edge label. 
  ///
  /// \relates LemonWriter
  template <typename _Graph>
  class UEdgeWriter : public LemonWriter::SectionWriter {
    typedef LemonWriter::SectionWriter Parent;
    typedef _Graph Graph;
    typedef typename Graph::Node Node;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::UEdge UEdge;
  public:
    
    /// \brief Constructor.
    ///
    /// Constructor for UEdgeWriter. It creates the UEdgeWriter and
    /// attach it into the given LemonWriter. The given \c _LabelWriter
    /// will write the undirected edges' label what can be an undirected 
    /// edgeset writer.
    template <typename _LabelWriter>
    UEdgeWriter(LemonWriter& _writer, const _LabelWriter& _labelWriter, 
	       const std::string& _name = std::string()) 
      : Parent(_writer), name(_name) {
      checkConcept<_writer_bits::ItemLabelWriter<Edge>, _LabelWriter>();
      checkConcept<_writer_bits::ItemLabelWriter<UEdge>, _LabelWriter>();
      uEdgeLabelWriter.reset(new _writer_bits::
			      LabelWriter<UEdge, _LabelWriter>(_labelWriter));
      edgeLabelWriter.reset(new _writer_bits::
			 LabelWriter<Edge, _LabelWriter>(_labelWriter));
    }

    /// \brief Destructor.
    ///
    /// Destructor for UEdgeWriter.
    virtual ~UEdgeWriter() {}
  private:
    UEdgeWriter(const UEdgeWriter&);
    void operator=(const UEdgeWriter&);

  public:

    /// \brief Add an edge writer command for the UEdgeWriter.
    ///
    /// Add an edge writer command for the UEdgeWriter.
    void writeEdge(std::string label, const Edge& item) {
      edgeWriters.push_back(make_pair(label, &item));
    }

    /// \brief Add an undirected edge writer command for the UEdgeWriter.
    ///
    /// Add an undirected edge writer command for the UEdgeWriter.
    void writeUEdge(std::string label, const UEdge& item) {
      uEdgeWriters.push_back(make_pair(label, &item));
    }

  protected:

    /// \brief The header of the section.
    ///
    /// It gives back the header of the section.
    virtual std::string header() {
      return "@uedges " + name;
    }

    /// \brief  Writer function of the section.
    ///
    /// Write the content of the section.
    virtual void write(std::ostream& os) {
      if (!edgeLabelWriter->isLabelWriter()) {
	throw DataFormatError("Cannot find undirected edgeset or label map");
      }
      if (!uEdgeLabelWriter->isLabelWriter()) {
	throw DataFormatError("Cannot find undirected edgeset or label map");
      }
      for (int i = 0; i < int(uEdgeWriters.size()); ++i) {
	os << uEdgeWriters[i].first << ' ';
	uEdgeLabelWriter->write(os, *(uEdgeWriters[i].second));
	os << std::endl;
      }
      for (int i = 0; i < int(edgeWriters.size()); ++i) {
	os << edgeWriters[i].first << ' ';
	edgeLabelWriter->write(os, *(edgeWriters[i].second));
	os << std::endl;
      }
    }

    /// \brief Gives back true when the section should be written.
    ///
    /// Gives back true when the section should be written.
    virtual bool valid() { 
      return !uEdgeWriters.empty() || !edgeWriters.empty(); 
    }
    
  private:

    std::string name;

    typedef std::vector<std::pair<std::string, 
				  const UEdge*> > UEdgeWriters;
    UEdgeWriters uEdgeWriters;
    std::auto_ptr<_writer_bits::LabelWriterBase<UEdge> > uEdgeLabelWriter;

    typedef std::vector<std::pair<std::string, const Edge*> > EdgeWriters;
    EdgeWriters edgeWriters;
    std::auto_ptr<_writer_bits::LabelWriterBase<Edge> > edgeLabelWriter;

  };

  /// \ingroup section_io
  /// \brief SectionWriter for writing extra node maps.
  ///
  /// The lemon format can store maps in the nodeset. This class let
  /// you make distinict section to store maps. The main purpose of
  /// this class is a logical separation of some maps. The other
  /// useful application could be to store paths in node maps.
  ///
  /// The first line of the section contains the names of the maps
  /// separated with white spaces. Each next line describes an item
  /// in the itemset, and contains in the first column the label of
  /// the item and then the mapped values for each map.
  ///
  /// \relates LemonWriter
  template <typename _Graph, typename _Traits = DefaultWriterTraits>
  class NodeMapWriter : public LemonWriter::SectionWriter {
    typedef LemonWriter::SectionWriter Parent;
  public:

    typedef _Graph Graph;
    typedef _Traits Traits;
    typedef typename Graph::Node Node;

    /// \brief Constructor.
    ///
    /// Constructor for NodeMapWriter. It creates the NodeMapWriter and
    /// attach it into the given LemonWriter. If the the
    /// \c _forceSort is true then the writer will write the edges
    /// sorted by the labels.
    template <typename _LabelWriter>
    NodeMapWriter(LemonWriter& _writer, const Graph& _graph,
		 const _LabelWriter& _labelWriter,
		 const std::string& _name = std::string(),
		 bool _forceSort = true) 
      : Parent(_writer), graph(_graph), name(_name), forceSort(_forceSort) {
      checkConcept<_writer_bits::ItemLabelWriter<Node>, _LabelWriter>();
      labelWriter.reset(new _writer_bits::LabelWriter<Node, 
			_LabelWriter>(_labelWriter));
    }

    /// \brief Destructor.
    ///
    /// Destructor for NodeMapWriter.
    virtual ~NodeMapWriter() {
      typename MapWriters::iterator it;
      for (it = writers.begin(); it != writers.end(); ++it) {
	delete it->second;
      }
    }

  private:
    NodeMapWriter(const NodeMapWriter&);
    void operator=(const NodeMapWriter&);
  
  public:

    /// \brief Add a new node map writer command for the writer.
    ///
    /// Add a new node map writer command for the writer.
    template <typename Map>
    NodeMapWriter& writeNodeMap(std::string label, const Map& map) {
      return writeNodeMap<typename Traits::
	template Writer<typename Map::Value>, Map>(label, map);
    }

    /// \brief Add a new node map writer command for the writer.
    ///
    /// Add a new node map writer command for the writer.
    template <typename ItemWriter, typename Map>
    NodeMapWriter& writeNodeMap(std::string label, const Map& map, 
			   const ItemWriter& iw = ItemWriter()) {
      checkConcept<concepts::ReadMap<Node, typename Map::Value>, Map>();
      checkConcept<_writer_bits::ItemWriter<typename Map::Value>,ItemWriter>();
      writers.push_back(
	make_pair(label, new _writer_bits::
		  MapWriter<Node, Map, ItemWriter>(map, iw)));
      return *this;
    }

  protected:

    /// \brief The header of the section.
    ///
    /// It gives back the header of the section.
    virtual std::string header() {
      return "@nodemaps " + name;
    }

    /// \brief  Writer function of the section.
    ///
    /// Write the content of the section.
    virtual void write(std::ostream& os) {
      std::vector<Node> nodes;
      for (typename Graph::NodeIt it(graph); it != INVALID; ++it) {
        nodes.push_back(it);
      }
      if (forceSort) {
	labelWriter->sort(nodes);
      }
      os << '\t';
      for (int i = 0; i < int(writers.size()); ++i) {
	os << writers[i].first << '\t';
      }
      os << std::endl;
      for (typename std::vector<Node>::iterator it = nodes.begin();
           it != nodes.end(); ++it) {

	labelWriter->write(os, *it); os << '\t';
	for (int i = 0; i < int(writers.size()); ++i) {
	  writers[i].second->write(os, *it);
	  os << '\t';
	}
	os << std::endl;
      }
    }


  private:

    typedef std::vector<std::pair<std::string, _writer_bits::
				  MapWriterBase<Node>*> > MapWriters;
    MapWriters writers;

    _writer_bits::MapWriterBase<Node>* labelMap;

    const Graph& graph;   
    std::string name;
    bool forceSort;

    std::auto_ptr<_writer_bits::LabelWriterBase<Node> > labelWriter;
  };

  /// \ingroup section_io
  /// \brief SectionWriter for writing extra edge maps.
  ///
  /// The lemon format can store maps in the edgeset. This class let
  /// you make distinict section to store maps. The main purpose of
  /// this class is a logical separation of some maps. The other
  /// useful application could be to store paths in edge maps.
  ///
  /// The first line of the section contains the names of the maps
  /// separated with white spaces. Each next line describes an item
  /// in the itemset, and contains in the first column the label of
  /// the item and then the mapped values for each map.
  ///
  /// \relates LemonWriter
  template <typename _Graph, typename _Traits = DefaultWriterTraits>
  class EdgeMapWriter : public LemonWriter::SectionWriter {
    typedef LemonWriter::SectionWriter Parent;
  public:

    typedef _Graph Graph;
    typedef _Traits Traits;
    typedef typename Graph::Edge Edge;

    /// \brief Constructor.
    ///
    /// Constructor for EdgeMapWriter. It creates the EdgeMapWriter and
    /// attach it into the given LemonWriter. If the the
    /// \c _forceSort is true then the writer will write the edges
    /// sorted by the labels.
    template <typename _LabelWriter>
    EdgeMapWriter(LemonWriter& _writer, const Graph& _graph,
		 const _LabelWriter& _labelWriter,
		 const std::string& _name = std::string(),
		 bool _forceSort = true) 
      : Parent(_writer), graph(_graph), name(_name), forceSort(_forceSort) {
      checkConcept<_writer_bits::ItemLabelWriter<Edge>, _LabelWriter>();
      labelWriter.reset(new _writer_bits::LabelWriter<Edge, 
			_LabelWriter>(_labelWriter));
    }

    /// \brief Destructor.
    ///
    /// Destructor for EdgeMapWriter.
    virtual ~EdgeMapWriter() {
      typename MapWriters::iterator it;
      for (it = writers.begin(); it != writers.end(); ++it) {
	delete it->second;
      }
    }

  private:
    EdgeMapWriter(const EdgeMapWriter&);
    void operator=(const EdgeMapWriter&);
  
  public:

    /// \brief Add a new edge map writer command for the writer.
    ///
    /// Add a new edge map writer command for the writer.
    template <typename Map>
    EdgeMapWriter& writeEdgeMap(std::string label, const Map& map) {
      return writeEdgeMap<typename Traits::
	template Writer<typename Map::Value>, Map>(label, map);
    }

    /// \brief Add a new edge map writer command for the writer.
    ///
    /// Add a new edge map writer command for the writer.
    template <typename ItemWriter, typename Map>
    EdgeMapWriter& writeEdgeMap(std::string label, const Map& map, 
				const ItemWriter& iw = ItemWriter()) {
      checkConcept<concepts::ReadMap<Edge, typename Map::Value>, Map>();
      checkConcept<_writer_bits::ItemWriter<typename Map::Value>,ItemWriter>();
      writers.push_back(
	make_pair(label, new _writer_bits::
		  MapWriter<Edge, Map, ItemWriter>(map, iw)));
      return *this;
    }

  protected:

    /// \brief The header of the section.
    ///
    /// It gives back the header of the section.
    virtual std::string header() {
      return "@edgemaps " + name;
    }

    /// \brief  Writer function of the section.
    ///
    /// Write the content of the section.
    virtual void write(std::ostream& os) {
      std::vector<Edge> edges;
      for (typename Graph::EdgeIt it(graph); it != INVALID; ++it) {
        edges.push_back(it);
      }
      if (forceSort) {
	labelWriter->sort(edges);
      }
      os << '\t';
      for (int i = 0; i < int(writers.size()); ++i) {
	os << writers[i].first << '\t';
      }
      os << std::endl;
      for (typename std::vector<Edge>::iterator it = edges.begin();
           it != edges.end(); ++it) {

	labelWriter->write(os, *it); os << '\t';
	for (int i = 0; i < int(writers.size()); ++i) {
	  writers[i].second->write(os, *it);
	  os << '\t';
	}
	os << std::endl;
      }
    }


  private:

    typedef std::vector<std::pair<std::string, _writer_bits::
				  MapWriterBase<Edge>*> > MapWriters;
    MapWriters writers;

    _writer_bits::MapWriterBase<Edge>* labelMap;

    const Graph& graph;   
    std::string name;
    bool forceSort;

    std::auto_ptr<_writer_bits::LabelWriterBase<Edge> > labelWriter;
  };

  /// \ingroup section_io
  /// \brief SectionWriter for writing extra undirected edge maps.
  ///
  /// The lemon format can store maps in the uedgeset. This class let
  /// you make distinict section to store maps. The main purpose of
  /// this class is a logical separation of some maps. The other
  /// useful application could be to store paths in undirected edge
  /// maps.
  ///
  /// The first line of the section contains the names of the maps
  /// separated with white spaces. Each next line describes an item
  /// in the itemset, and contains in the first column the label of
  /// the item and then the mapped values for each map.
  ///
  /// \relates LemonWriter
  template <typename _Graph, typename _Traits = DefaultWriterTraits>
  class UEdgeMapWriter : public LemonWriter::SectionWriter {
    typedef LemonWriter::SectionWriter Parent;
  public:

    typedef _Graph Graph;
    typedef _Traits Traits;
    typedef typename Graph::UEdge UEdge;
    typedef typename Graph::Edge Edge;

    /// \brief Constructor.
    ///
    /// Constructor for UEdgeMapWriter. It creates the UEdgeMapWriter and
    /// attach it into the given LemonWriter. If the the
    /// \c _forceSort is true then the writer will write the uedges
    /// sorted by the labels.
    template <typename _LabelWriter>
    UEdgeMapWriter(LemonWriter& _writer, const Graph& _graph,
		 const _LabelWriter& _labelWriter,
		 const std::string& _name = std::string(),
		 bool _forceSort = true) 
      : Parent(_writer), graph(_graph), name(_name), forceSort(_forceSort) {
      checkConcept<_writer_bits::ItemLabelWriter<UEdge>, _LabelWriter>();
      labelWriter.reset(new _writer_bits::LabelWriter<UEdge, 
			    _LabelWriter>(_labelWriter));
    }

    /// \brief Destructor.
    ///
    /// Destructor for UEdgeMapWriter.
    virtual ~UEdgeMapWriter() {
      typename MapWriters::iterator it;
      for (it = writers.begin(); it != writers.end(); ++it) {
	delete it->second;
      }
    }

  private:
    UEdgeMapWriter(const UEdgeMapWriter&);
    void operator=(const UEdgeMapWriter&);
  
  public:

    /// \brief Add a new undirected edge map writer command for the writer.
    ///
    /// Add a new undirected edge map writer command for the writer.
    template <typename Map>
    UEdgeMapWriter& writeUEdgeMap(std::string label, const Map& map) {
      return writeUEdgeMap<typename Traits::
	template Writer<typename Map::Value>, Map>(label, map);
    }

    /// \brief Add a new undirected edge map writer command for the writer.
    ///
    /// Add a new undirected edge map writer command for the writer.
    template <typename ItemWriter, typename Map>
    UEdgeMapWriter& writeUEdgeMap(std::string label, const Map& map, 
			   const ItemWriter& iw = ItemWriter()) {
      checkConcept<concepts::ReadMap<UEdge, typename Map::Value>, Map>();
      checkConcept<_writer_bits::ItemWriter<typename Map::Value>,ItemWriter>();
      writers.push_back(
	make_pair(label, new _writer_bits::
		  MapWriter<UEdge, Map, ItemWriter>(map, iw)));
      return *this;
    }

    /// \brief Add a new directed edge map writer command for the writer.
    ///
    /// Add a new directed map writer command for the writer.
    template <typename Map>
    UEdgeMapWriter& writeEdgeMap(std::string label, const Map& map) {
      return writeEdgeMap<typename Traits::
	template Writer<typename Map::Value>, Map>(label, map);
    }

    /// \brief Add a new directed map writer command for the writer.
    ///
    /// Add a new directed map writer command for the writer.
    template <typename ItemWriter, typename Map>
    UEdgeMapWriter& writeEdgeMap(std::string label, const Map& map, 
                                 const ItemWriter& iw = ItemWriter()) {
      checkConcept<concepts::ReadMap<Edge, typename Map::Value>, Map>();
      checkConcept<_writer_bits::ItemWriter<typename Map::Value>, ItemWriter>();
      writeUEdgeMap("+" + label, 
                    _writer_bits::forwardComposeMap(graph, map), iw);
      writeUEdgeMap("-" + label, 
                    _writer_bits::backwardComposeMap(graph, map), iw);
      return *this;
    }

  protected:

    /// \brief The header of the section.
    ///
    /// It gives back the header of the section.
    virtual std::string header() {
      return "@uedgemaps " + name;
    }

    /// \brief  Writer function of the section.
    ///
    /// Write the content of the section.
    virtual void write(std::ostream& os) {
      std::vector<UEdge> uedges;
      for (typename Graph::UEdgeIt it(graph); it != INVALID; ++it) {
        uedges.push_back(it);
      }
      if (forceSort) {
	labelWriter->sort(uedges);
      }
      os << '\t';
      for (int i = 0; i < int(writers.size()); ++i) {
	os << writers[i].first << '\t';
      }
      os << std::endl;
      for (typename std::vector<UEdge>::iterator it = uedges.begin();
           it != uedges.end(); ++it) {

	labelWriter->write(os, *it); os << '\t';
	for (int i = 0; i < int(writers.size()); ++i) {
	  writers[i].second->write(os, *it);
	  os << '\t';
	}
	os << std::endl;
      }
    }


  private:

    typedef std::vector<std::pair<std::string, _writer_bits::
				  MapWriterBase<UEdge>*> > MapWriters;
    MapWriters writers;

    _writer_bits::MapWriterBase<UEdge>* labelMap;

    const Graph& graph;   
    std::string name;
    bool forceSort;

    std::auto_ptr<_writer_bits::LabelWriterBase<UEdge> > labelWriter;
  };


  /// \ingroup section_io
  /// \brief SectionWriter for attributes.
  ///
  /// The lemon format can store multiple attribute set. Each set has
  /// the header line \c \@attributes \c attributes_name, but the 
  /// attributeset_name may be empty.
  ///
  /// The attributeset section contains several lines. Each of them starts
  /// with the name of attribute and then the value.
  ///
  /// \relates LemonWriter
  template <typename _Traits = DefaultWriterTraits>
  class AttributeWriter : public LemonWriter::SectionWriter {
    typedef LemonWriter::SectionWriter Parent;
    typedef _Traits Traits; 
  public:
    /// \brief Constructor.
    ///
    /// Constructor for AttributeWriter. It creates the AttributeWriter and
    /// attach it into the given LemonWriter.
    AttributeWriter(LemonWriter& _writer, 
		    const std::string& _name = std::string()) 
      : Parent(_writer), name(_name) {}

    /// \brief Destructor.
    ///
    /// Destructor for AttributeWriter.
    virtual ~AttributeWriter() {
      typename Writers::iterator it;
      for (it = writers.begin(); it != writers.end(); ++it) {
	delete it->second;
      }
    }

  private:
    AttributeWriter(const AttributeWriter&);
    void operator=(AttributeWriter&);

  public:
    /// \brief Add an attribute writer command for the writer.
    ///
    /// Add an attribute writer command for the writer.
    template <typename Value>
    AttributeWriter& writeAttribute(std::string label, 
				    const Value& value) {
      return 
	writeAttribute<typename Traits::template Writer<Value> >(label, value);
    }

    /// \brief Add an attribute writer command for the writer.
    ///
    /// Add an attribute writer command for the writer.
    template <typename ItemWriter, typename Value>
    AttributeWriter& writeAttribute(std::string label, const Value& value,
				    const ItemWriter& iw = ItemWriter()) {
      checkConcept<_writer_bits::ItemWriter<Value>, ItemWriter>();
      writers.push_back(make_pair(label, new _writer_bits::
				  ValueWriter<Value, ItemWriter>(value, iw)));
      return *this;
    }

  protected:

    /// \brief The header of section.
    ///
    /// It gives back the header of the section.
    std::string header() {
      return "@attributes " + name;
    }

    /// \brief  Writer function of the section.
    ///
    /// Write the content of the section.
    void write(std::ostream& os) {
      typename Writers::iterator it;
      for (it = writers.begin(); it != writers.end(); ++it) {
	os << it->first << ' ';
	it->second->write(os);
	os << std::endl;
      }
    }    

    /// \brief Gives back true when the section should be written.
    ///
    /// Gives back true when the section should be written.
    virtual bool valid() { return !writers.empty(); }

  private:
    std::string name;

    typedef std::vector<std::pair<std::string, 
				  _writer_bits::ValueWriterBase*> > Writers;
    Writers writers;  
  };


}
#endif
