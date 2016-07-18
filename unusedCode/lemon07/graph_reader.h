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
///\brief Lemon Graph Format reader.

#ifndef LEMON_GRAPH_READER_H
#define LEMON_GRAPH_READER_H

#include <iostream>

#include <lemon/error.h>
#include <lemon/lemon_reader.h>

namespace lemon {

  /// \addtogroup lemon_io
  /// @{

  /// \brief The graph reader class.
  ///
  /// The \c GraphReader class provides the graph input. 
  /// Before you read this documentation it might be useful to read the general
  /// description of  \ref graph-io-page "Graph Input-Output".
  ///
  /// The file to be read may contain several maps and labeled
  /// (designated) nodes or edges.
  ///
  /// If you read a graph you need not read all the maps and items just those
  /// that you need. The interface of the \c GraphReader is very similar to
  /// the GraphWriter but the reading method does not depend on the order the
  /// given commands (i.e. you don't have to insist on the order in which the
  /// maps are given in the file).
  ///
  /// The reader object assumes that not read values do not contain 
  /// whitespaces, therefore it has some extra possibilities to control how
  /// it should skip the values when the string representation contains spaces.
  ///
  ///\code
  /// GraphReader<ListGraph> reader(std::cin, graph);
  ///\endcode
  ///
  /// The \c readNodeMap() function reads a map from the \c \@nodeset section.
  /// If there is a map that you do not want to read from the file and there is
  /// whitespace in the string represenation of the values then you should
  /// call the \c skipNodeMap() template member function with proper 
  /// parameters.
  ///
  ///\code
  /// reader.readNodeMap("coords", coords);
  ///
  /// reader.skipNodeMap("description", desc);
  ///
  /// reader.readNodeMap("color", colorMap);
  ///\endcode
  ///
  /// With the \c readEdgeMap() member function you can give an edge map
  /// reading command similar to the NodeMaps. 
  ///
  ///\code
  /// reader.readEdgeMap("weight", weightMap);
  /// reader.readEdgeMap("label", labelMap);
  ///\endcode
  ///
  /// With \c readNode() and \c readEdge() functions you can read 
  /// labeled Nodes and Edges.
  ///
  ///\code
  /// reader.readNode("source", sourceNode);
  /// reader.readNode("target", targetNode);
  ///
  /// reader.readEdge("observed", edge);
  ///\endcode
  ///
  /// With the \c readAttribute() functions you can read an attribute
  /// into a variable. You can specify the reader for the attribute as
  /// the nodemaps.
  ///
  /// After you give all read commands you must call the \c run() member
  /// function, which executes all the commands.
  ///
  ///\code
  /// reader.run();
  ///\endcode
  ///
  /// \see DefaultReaderTraits
  /// \see QuotedStringReader
  /// \see \ref GraphWriter
  /// \see \ref graph-io-page
  /// \author Balazs Dezso
  template <typename _Graph, typename _ReaderTraits = DefaultReaderTraits> 
  class GraphReader {
  public:
    
    typedef _Graph Graph;
    typedef typename Graph::Node Node;
    typedef typename Graph::Edge Edge;

    typedef _ReaderTraits ReaderTraits;
    typedef typename ReaderTraits::Skipper DefaultSkipper;

    /// \brief Construct a new GraphReader.
    ///
    /// Construct a new GraphReader. It reads into the given graph
    /// and it uses the given reader as the default skipper.
    GraphReader(std::istream& _is, Graph& _graph, 
		const DefaultSkipper& _skipper = DefaultSkipper()) 
      : reader(new LemonReader(_is)), own_reader(true), skipper(_skipper),
	nodeset_reader(*reader, _graph, std::string(), skipper),
	edgeset_reader(*reader, _graph, nodeset_reader, 
		       std::string(), skipper),
	node_reader(*reader, nodeset_reader, std::string()),
	edge_reader(*reader, edgeset_reader, std::string()),
	attribute_reader(*reader, std::string()) {}

    /// \brief Construct a new GraphReader.
    ///
    /// Construct a new GraphReader. It reads into the given graph
    /// and it uses the given reader as the default skipper.
    GraphReader(const std::string& _filename, Graph& _graph, 
		const DefaultSkipper& _skipper = DefaultSkipper()) 
      : reader(new LemonReader(_filename)), own_reader(true), 
	skipper(_skipper),
	nodeset_reader(*reader, _graph, std::string(), skipper),
	edgeset_reader(*reader, _graph, nodeset_reader, 
		       std::string(), skipper),
	node_reader(*reader, nodeset_reader, std::string()),
	edge_reader(*reader, edgeset_reader, std::string()),
	attribute_reader(*reader, std::string()) {}

    /// \brief Construct a new GraphReader.
    ///
    /// Construct a new GraphReader. It reads into the given graph
    /// and it uses the given reader as the default skipper.
    GraphReader(LemonReader& _reader, Graph& _graph, 
		const DefaultSkipper& _skipper = DefaultSkipper()) 
      : reader(_reader), own_reader(false), skipper(_skipper),
	nodeset_reader(*reader, _graph, std::string(), skipper),
	edgeset_reader(*reader, _graph, nodeset_reader, 
		       std::string(), skipper),
	node_reader(*reader, nodeset_reader, std::string()),
	edge_reader(*reader, edgeset_reader, std::string()),
	attribute_reader(*reader, std::string()) {}

    /// \brief Destruct the graph reader.
    ///
    /// Destruct the graph reader.
    ~GraphReader() {
      if (own_reader) 
	delete reader;
    }

    /// \brief Give a new node map reading command to the reader.
    ///
    /// Give a new node map reading command to the reader.
    template <typename Map>
    GraphReader& readNodeMap(std::string name, Map& map) {
      nodeset_reader.readNodeMap(name, map);
      return *this;
    }

    template <typename Map>
    GraphReader& readNodeMap(std::string name, const Map& map) {
      nodeset_reader.readNodeMap(name, map);
      return *this;
    }

    /// \brief Give a new node map reading command to the reader.
    ///
    /// Give a new node map reading command to the reader.
    template <typename ItemReader, typename Map>
    GraphReader& readNodeMap(std::string name, Map& map, 
			     const ItemReader& ir = ItemReader()) {
      nodeset_reader.readNodeMap(name, map, ir);
      return *this;
    }

    template <typename ItemReader, typename Map>
    GraphReader& readNodeMap(std::string name, const Map& map, 
			     const ItemReader& ir = ItemReader()) {
      nodeset_reader.readNodeMap(name, map, ir);
      return *this;
    }

    /// \brief Give a new node map skipping command to the reader.
    ///
    /// Give a new node map skipping command to the reader.
    template <typename ItemReader>
    GraphReader& skipNodeMap(std::string name, 
			     const ItemReader& ir = ItemReader()) {
      nodeset_reader.skipNodeMap(name, ir);
      return *this;
    }

    /// \brief Give a new edge map reading command to the reader.
    ///
    /// Give a new edge map reading command to the reader.
    template <typename Map>
    GraphReader& readEdgeMap(std::string name, Map& map) { 
      edgeset_reader.readEdgeMap(name, map);
      return *this;
    }

    template <typename Map>
    GraphReader& readEdgeMap(std::string name, const Map& map) { 
      edgeset_reader.readEdgeMap(name, map);
      return *this;
    }


    /// \brief Give a new edge map reading command to the reader.
    ///
    /// Give a new edge map reading command to the reader.
    template <typename ItemReader, typename Map>
    GraphReader& readEdgeMap(std::string name, Map& map,
			     const ItemReader& ir = ItemReader()) {
      edgeset_reader.readEdgeMap(name, map, ir);
      return *this;
    }

    template <typename ItemReader, typename Map>
    GraphReader& readEdgeMap(std::string name, const Map& map,
			     const ItemReader& ir = ItemReader()) {
      edgeset_reader.readEdgeMap(name, map, ir);
      return *this;
    }

    /// \brief Give a new edge map skipping command to the reader.
    ///
    /// Give a new edge map skipping command to the reader.
    template <typename ItemReader>
    GraphReader& skipEdgeMap(std::string name, 
			     const ItemReader& ir = ItemReader()) {
      edgeset_reader.skipEdgeMap(name, ir);
      return *this;
    }

    /// \brief Give a new labeled node reading command to the reader.
    ///
    /// Give a new labeled node reading command to the reader.
    GraphReader& readNode(std::string name, Node& node) {
      node_reader.readNode(name, node);
      return *this;
    }

    /// \brief Give a new labeled edge reading command to the reader.
    ///
    /// Give a new labeled edge reading command to the reader.
    GraphReader& readEdge(std::string name, Edge& edge) {
      edge_reader.readEdge(name, edge);
      return *this;
    }

    /// \brief Give a new attribute reading command.
    ///
    ///  Give a new attribute reading command.
    template <typename Value>
    GraphReader& readAttribute(std::string name, Value& value) {
      attribute_reader.readAttribute(name, value);
      return *this;
    }
    
    /// \brief Give a new attribute reading command.
    ///
    ///  Give a new attribute reading command.
    template <typename ItemReader, typename Value>
    GraphReader& readAttribute(std::string name, Value& value, 
			       const ItemReader& ir = ItemReader()) {
      attribute_reader.readAttribute(name, value, ir);
      return *this;
    }

    /// \brief Conversion operator to LemonReader.
    ///
    /// Conversion operator to LemonReader. It makes possible to access the
    /// encapsulated \e LemonReader, this way you can attach to this reader
    /// new instances of \e LemonReader::SectionReader. For more details see
    /// the \ref rwbackground "Background of Reading and Writing".
    operator LemonReader&() {
      return *reader;
    }

    /// \brief Executes the reading commands.
    ///
    /// Executes the reading commands.
    void run() {
      reader->run();
    }


    /// \brief Returns true if the reader can give back the items by its label.
    ///
    /// \brief Returns true if the reader can give back the items by its label.
    bool isLabelReader() const {
      return nodeset_reader.isLabelReader() && edgeset_reader.isLabelReader();
    }

    /// \brief Gives back the node by its label.
    ///
    /// It reads an label from the stream and gives back which node belongs to
    /// it. It is possible only if there was read a "label" named node map.
    void readLabel(std::istream& is, Node& node) const {
      nodeset_reader.readLabel(is, node);
    } 

    /// \brief Gives back the edge by its label.
    ///
    /// It reads an label from the stream and gives back which edge belongs to
    /// it. It is possible only if there was read a "label" named edge map.
    void readLabel(std::istream& is, Edge& edge) const {
      edgeset_reader.readLabel(is, edge);
    } 

  private:

    LemonReader* reader;
    bool own_reader;

    DefaultSkipper skipper;

    NodeSetReader<Graph, ReaderTraits> nodeset_reader;
    EdgeSetReader<Graph, ReaderTraits> edgeset_reader;

    NodeReader<Graph> node_reader;
    EdgeReader<Graph> edge_reader;
    
    AttributeReader<ReaderTraits> attribute_reader;
  };


  /// \brief The undirected graph reader class.
  ///
  /// The \c UGraphReader class provides the graph input. 
  /// Before you read this documentation it might be useful to read the general
  /// description of  \ref graph-io-page "Graph Input-Output".
  ///
  /// The given file format may contain several maps and labeled nodes or 
  /// edges.
  ///
  /// If you read a graph you need not read all the maps and items just those
  /// that you need. The interface of the \c UGraphReader is very similar
  /// to the UGraphWriter but the reading method does not depend on the
  /// order of the given commands.
  ///
  /// The reader object suppose that each not read value does not contain 
  /// whitespaces, therefore it has some extra possibilities to control how
  /// it should skip the values when the string representation contains spaces.
  ///
  ///\code
  /// UGraphReader<ListUGraph> reader(std::cin, graph);
  ///\endcode
  ///
  /// The \c readNodeMap() function reads a map from the \c \@nodeset section.
  /// If there is a map that you do not want to read from the file and there is
  /// whitespace in the string represenation of the values then you should
  /// call the \c skipNodeMap() template member function with proper 
  /// parameters.
  ///
  ///\code
  /// reader.readNodeMap("coords", coords);
  ///
  /// reader.skipNodeMap("description", desc);
  ///
  /// reader.readNodeMap("color", colorMap);
  ///\endcode
  ///
  /// With the \c readUEdgeMap() member function you can give an 
  /// uedge map reading command similar to the NodeMaps. 
  ///
  ///\code
  /// reader.readUEdgeMap("capacity", capacityMap);
  ///\endcode
  ///
  /// The reading of the directed edge maps is just a syntactical sugar.
  /// It reads two undirected edgemaps into a directed edge map. The 
  /// undirected edge maps' name should be start with the \c '+' and the
  /// \c '-' character and the same.
  ///
  ///\code
  /// reader.readEdgeMap("flow", flowMap);
  ///\endcode 
  ///
  /// With \c readNode() and \c readUEdge() functions you can read 
  /// labeled Nodes and UEdges.
  ///
  ///\code
  /// reader.readNode("source", sourceNode);
  /// reader.readNode("target", targetNode);
  ///
  /// reader.readUEdge("observed", uEdge);
  ///\endcode
  ///
  /// With the \c readAttribute() functions you can read an attribute
  /// in a variable. You can specify the reader for the attribute as
  /// the nodemaps.
  ///
  /// After you give all read commands you must call the \c run() member
  /// function, which execute all the commands.
  ///
  ///\code
  /// reader.run();
  ///\endcode
  ///
  /// \see GraphReader
  /// \see DefaultReaderTraits
  /// \see \ref UGraphWriter
  /// \see \ref graph-io-page
  ///
  /// \author Balazs Dezso
  template <typename _Graph, typename _ReaderTraits = DefaultReaderTraits> 
  class UGraphReader {
  public:
    
    typedef _Graph Graph;
    typedef typename Graph::Node Node;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::UEdge UEdge;

    typedef _ReaderTraits ReaderTraits;
    typedef typename ReaderTraits::Skipper DefaultSkipper;

    /// \brief Construct a new UGraphReader.
    ///
    /// Construct a new UGraphReader. It reads into the given graph
    /// and it use the given reader as the default skipper.
    UGraphReader(std::istream& _is, Graph& _graph, 
		     const DefaultSkipper& _skipper = DefaultSkipper()) 
      : reader(new LemonReader(_is)), own_reader(true), skipper(_skipper),
	nodeset_reader(*reader, _graph, std::string(), skipper),
	uedgeset_reader(*reader, _graph, nodeset_reader, 
			     std::string(), skipper),
	node_reader(*reader, nodeset_reader, std::string()),
	uedge_reader(*reader, uedgeset_reader, std::string()),
	attribute_reader(*reader, std::string()) {}

    /// \brief Construct a new UGraphReader.
    ///
    /// Construct a new UGraphReader. It reads into the given graph
    /// and it use the given reader as the default skipper.
    UGraphReader(const std::string& _filename, Graph& _graph, 
		     const DefaultSkipper& _skipper = DefaultSkipper()) 
      : reader(new LemonReader(_filename)), own_reader(true), 
	skipper(_skipper),
	nodeset_reader(*reader, _graph, std::string(), skipper),
	uedgeset_reader(*reader, _graph, nodeset_reader, 
			     std::string(), skipper),
	node_reader(*reader, nodeset_reader, std::string()),
	uedge_reader(*reader, uedgeset_reader, std::string()),
	attribute_reader(*reader, std::string()) {}

    /// \brief Construct a new UGraphReader.
    ///
    /// Construct a new UGraphReader. It reads into the given graph
    /// and it use the given reader as the default skipper.
    UGraphReader(LemonReader& _reader, Graph& _graph, 
		     const DefaultSkipper& _skipper = DefaultSkipper()) 
      : reader(_reader), own_reader(false), skipper(_skipper),
	nodeset_reader(*reader, _graph, std::string(), skipper),
	uedgeset_reader(*reader, _graph, nodeset_reader, 
			     std::string(), skipper),
	node_reader(*reader, nodeset_reader, std::string()),
	uedge_reader(*reader, uedgeset_reader, std::string()),
	attribute_reader(*reader, std::string()) {}

    /// \brief Destruct the graph reader.
    ///
    /// Destruct the graph reader.
    ~UGraphReader() {
      if (own_reader) 
	delete reader;
    }

    /// \brief Give a new node map reading command to the reader.
    ///
    /// Give a new node map reading command to the reader.
    template <typename Map>
    UGraphReader& readNodeMap(std::string name, Map& map) {
      nodeset_reader.readNodeMap(name, map);
      return *this;
    }

    template <typename Map>
    UGraphReader& readNodeMap(std::string name, const Map& map) {
      nodeset_reader.readNodeMap(name, map);
      return *this;
    }

    /// \brief Give a new node map reading command to the reader.
    ///
    /// Give a new node map reading command to the reader.
    template <typename ItemReader, typename Map>
    UGraphReader& readNodeMap(std::string name, Map& map, 
                              const ItemReader& ir = ItemReader()) {
      nodeset_reader.readNodeMap(name, map, ir);
      return *this;
    }

    template <typename ItemReader, typename Map>
    UGraphReader& readNodeMap(std::string name, const Map& map, 
                              const ItemReader& ir = ItemReader()) {
      nodeset_reader.readNodeMap(name, map, ir);
      return *this;
    }

    /// \brief Give a new node map skipping command to the reader.
    ///
    /// Give a new node map skipping command to the reader.
    template <typename ItemReader>
    UGraphReader& skipNodeMap(std::string name, 
                              const ItemReader& ir = ItemReader()) {
      nodeset_reader.skipNodeMap(name, ir);
      return *this;
    }

    /// \brief Give a new undirected edge map reading command to the reader.
    ///
    /// Give a new undirected edge map reading command to the reader.
    template <typename Map>
    UGraphReader& readUEdgeMap(std::string name, Map& map) { 
      uedgeset_reader.readUEdgeMap(name, map);
      return *this;
    }

    template <typename Map>
    UGraphReader& readUEdgeMap(std::string name, const Map& map) { 
      uedgeset_reader.readUEdgeMap(name, map);
      return *this;
    }


    /// \brief Give a new undirected edge map reading command to the reader.
    ///
    /// Give a new undirected edge map reading command to the reader.
    template <typename ItemReader, typename Map>
    UGraphReader& readUEdgeMap(std::string name, Map& map,
                               const ItemReader& ir = ItemReader()) {
      uedgeset_reader.readUEdgeMap(name, map, ir);
      return *this;
    }

    template <typename ItemReader, typename Map>
    UGraphReader& readUEdgeMap(std::string name, const Map& map,
                               const ItemReader& ir = ItemReader()) {
      uedgeset_reader.readUEdgeMap(name, map, ir);
      return *this;
    }

    /// \brief Give a new undirected edge map skipping command to the reader.
    ///
    /// Give a new undirected edge map skipping command to the reader.
    template <typename ItemReader>
    UGraphReader& skipUEdgeMap(std::string name,
				       const ItemReader& ir = ItemReader()) {
      uedgeset_reader.skipUMap(name, ir);
      return *this;
    }


    /// \brief Give a new edge map reading command to the reader.
    ///
    /// Give a new edge map reading command to the reader.
    template <typename Map>
    UGraphReader& readEdgeMap(std::string name, Map& map) { 
      uedgeset_reader.readEdgeMap(name, map);
      return *this;
    }

    template <typename Map>
    UGraphReader& readEdgeMap(std::string name, const Map& map) { 
      uedgeset_reader.readEdgeMap(name, map);
      return *this;
    }


    /// \brief Give a new edge map reading command to the reader.
    ///
    /// Give a new edge map reading command to the reader.
    template <typename ItemReader, typename Map>
    UGraphReader& readEdgeMap(std::string name, Map& map,
                              const ItemReader& ir = ItemReader()) {
      uedgeset_reader.readEdgeMap(name, map, ir);
      return *this;
    }

    template <typename ItemReader, typename Map>
    UGraphReader& readEdgeMap(std::string name, const Map& map,
                              const ItemReader& ir = ItemReader()) {
      uedgeset_reader.readEdgeMap(name, map, ir);
      return *this;
    }

    /// \brief Give a new edge map skipping command to the reader.
    ///
    /// Give a new edge map skipping command to the reader.
    template <typename ItemReader>
    UGraphReader& skipEdgeMap(std::string name,
                              const ItemReader& ir = ItemReader()) {
      uedgeset_reader.skipEdgeMap(name, ir);
      return *this;
    }

    /// \brief Give a new labeled node reading command to the reader.
    ///
    /// Give a new labeled node reading command to the reader.
    UGraphReader& readNode(std::string name, Node& node) {
      node_reader.readNode(name, node);
      return *this;
    }

    /// \brief Give a new labeled edge reading command to the reader.
    ///
    /// Give a new labeled edge reading command to the reader.
    UGraphReader& readEdge(std::string name, Edge& edge) {
      uedge_reader.readEdge(name, edge);
    }

    /// \brief Give a new labeled undirected edge reading command to the
    /// reader.
    ///
    /// Give a new labeled undirected edge reading command to the reader.
    UGraphReader& readUEdge(std::string name, UEdge& edge) {
      uedge_reader.readUEdge(name, edge);
    }

    /// \brief Give a new attribute reading command.
    ///
    ///  Give a new attribute reading command.
    template <typename Value>
    UGraphReader& readAttribute(std::string name, Value& value) {
      attribute_reader.readAttribute(name, value);
      return *this;
    }
    
    /// \brief Give a new attribute reading command.
    ///
    ///  Give a new attribute reading command.
    template <typename ItemReader, typename Value>
    UGraphReader& readAttribute(std::string name, Value& value, 
			       const ItemReader& ir = ItemReader()) {
      attribute_reader.readAttribute(name, value, ir);
      return *this;
    }

    /// \brief Conversion operator to LemonReader.
    ///
    /// Conversion operator to LemonReader. It make possible
    /// to access the encapsulated \e LemonReader, this way
    /// you can attach to this reader new instances of 
    /// \e LemonReader::SectionReader.
    operator LemonReader&() {
      return *reader;
    }

    /// \brief Executes the reading commands.
    ///
    /// Executes the reading commands.
    void run() {
      reader->run();
    }


    /// \brief Returns true if the reader can give back the items by its label.
    ///
    /// Returns true if the reader can give back the items by its label.
    bool isLabelReader() const {
      return nodeset_reader.isLabelReader() && 
        uedgeset_reader.isLabelReader();
    }

    /// \brief Gives back the node by its label.
    ///
    /// It reads an label from the stream and gives back which node belongs to
    /// it. It is possible only if there was read a "label" named node map.
    void readLabel(std::istream& is, Node& node) const {
      return nodeset_reader.readLabel(is, node);
    } 

    /// \brief Gives back the edge by its label
    ///
    /// It reads an label from the stream and gives back which edge belongs to
    /// it. It is possible only if there was read a "label" named edge map.
    void readLabel(std::istream& is, Edge& edge) const {
      return uedgeset_reader.readLabel(is, edge);
    } 

    /// \brief Gives back the undirected edge by its label.
    ///
    /// It reads an label from the stream and gives back which undirected edge 
    /// belongs to it. It is possible only if there was read a "label" named 
    /// edge map.
    void readLabel(std::istream& is, UEdge& uedge) const {
      return uedgeset_reader.readLabel(is, uedge);
    } 
    

  private:

    LemonReader* reader;
    bool own_reader;

    DefaultSkipper skipper;

    NodeSetReader<Graph, ReaderTraits> nodeset_reader;
    UEdgeSetReader<Graph, ReaderTraits> uedgeset_reader;

    NodeReader<Graph> node_reader;
    UEdgeReader<Graph> uedge_reader;
    
    AttributeReader<ReaderTraits> attribute_reader;
  };

  /// \brief The bipartite graph reader class.
  ///
  /// The \c BpUGraphReader class provides the graph input. 
  /// Before you read this documentation it might be useful to read the general
  /// description of  \ref graph-io-page "Graph Input-Output".
  ///
  /// The given file format may contain several maps and labeled nodes or 
  /// edges.
  ///
  /// If you read a graph you need not read all the maps and items just those
  /// that you need. The interface of the \c BpUGraphReader is very similar
  /// to the BpUGraphWriter but the reading method does not depend on the
  /// order of the given commands.
  ///
  /// The reader object suppose that each not read value does not contain 
  /// whitespaces, therefore it has some extra possibilities to control how
  /// it should skip the values when the string representation contains spaces.
  ///
  ///\code
  /// BpUGraphReader<ListBpUGraph> reader(std::cin, graph);
  ///\endcode
  ///
  /// The \c readANodeMap() function reads a map from the A-part of
  /// the\c \@bpnodeset section, while the \c readBNodeMap() reads
  /// from the B-part of the section.  If you use the \c readNodeMap()
  /// function, then the given map should appear in both part of the
  /// section. If there is a map that you do not want to read from the
  /// file and there is whitespace in the string represenation of the
  /// values then you should call the \c skipANodeMap(), \c
  /// skipBNodeMap() or \c skipNodeMap() template member function with
  /// proper parameters.
  ///
  ///\code
  /// reader.readNodeMap("coords", coords);
  /// reader.readANodeMap("range", range);
  /// reader.readANodeMap("benefit", benefit);
  ///
  /// reader.skipNodeMap("description", desc);
  ///
  /// reader.readNodeMap("color", colorMap);
  ///\endcode
  ///
  /// With the \c readUEdgeMap() member function you can give an 
  /// uedge map reading command similar to the NodeMaps. 
  ///
  ///\code
  /// reader.readUEdgeMap("capacity", capacityMap);
  /// reader.readEdgeMap("flow", flowMap);
  ///\endcode 
  ///
  /// With \c readNode() and \c readUEdge() functions you can read 
  /// labeled Nodes and UEdges.
  ///
  ///\code
  /// reader.readNode("source", sourceNode);
  /// reader.readNode("target", targetNode);
  ///
  /// reader.readUEdge("observed", uEdge);
  ///\endcode
  ///
  /// With the \c readAttribute() functions you can read an attribute
  /// in a variable. You can specify the reader for the attribute as
  /// the nodemaps.
  ///
  /// After you give all read commands you must call the \c run() member
  /// function, which execute all the commands.
  ///
  ///\code
  /// reader.run();
  ///\endcode
  ///
  /// \see GraphReader
  /// \see DefaultReaderTraits
  /// \see \ref UGraphWriter
  /// \see \ref graph-io-page
  ///
  /// \author Balazs Dezso
  template <typename _Graph, typename _ReaderTraits = DefaultReaderTraits> 
  class BpUGraphReader {
  public:
    
    typedef _Graph Graph;
    typedef typename Graph::Node Node;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::UEdge UEdge;

    typedef _ReaderTraits ReaderTraits;
    typedef typename ReaderTraits::Skipper DefaultSkipper;

    /// \brief Construct a new BpUGraphReader.
    ///
    /// Construct a new BpUGraphReader. It reads into the given graph
    /// and it use the given reader as the default skipper.
    BpUGraphReader(std::istream& _is, Graph& _graph, 
		     const DefaultSkipper& _skipper = DefaultSkipper()) 
      : reader(new LemonReader(_is)), own_reader(true), skipper(_skipper),
	nodeset_reader(*reader, _graph, std::string(), skipper),
	uedgeset_reader(*reader, _graph, nodeset_reader, 
			     std::string(), skipper),
	node_reader(*reader, nodeset_reader, std::string()),
	uedge_reader(*reader, uedgeset_reader, std::string()),
	attribute_reader(*reader, std::string()) {}

    /// \brief Construct a new BpUGraphReader.
    ///
    /// Construct a new BpUGraphReader. It reads into the given graph
    /// and it use the given reader as the default skipper.
    BpUGraphReader(const std::string& _filename, Graph& _graph, 
		     const DefaultSkipper& _skipper = DefaultSkipper()) 
      : reader(new LemonReader(_filename)), own_reader(true), 
	skipper(_skipper),
	nodeset_reader(*reader, _graph, std::string(), skipper),
	uedgeset_reader(*reader, _graph, nodeset_reader, 
			     std::string(), skipper),
	node_reader(*reader, nodeset_reader, std::string()),
	uedge_reader(*reader, uedgeset_reader, std::string()),
	attribute_reader(*reader, std::string()) {}

    /// \brief Construct a new BpUGraphReader.
    ///
    /// Construct a new BpUGraphReader. It reads into the given graph
    /// and it use the given reader as the default skipper.
    BpUGraphReader(LemonReader& _reader, Graph& _graph, 
		     const DefaultSkipper& _skipper = DefaultSkipper()) 
      : reader(_reader), own_reader(false), skipper(_skipper),
	nodeset_reader(*reader, _graph, std::string(), skipper),
	uedgeset_reader(*reader, _graph, nodeset_reader, 
			     std::string(), skipper),
	node_reader(*reader, nodeset_reader, std::string()),
	uedge_reader(*reader, uedgeset_reader, std::string()),
	attribute_reader(*reader, std::string()) {}

    /// \brief Destruct the graph reader.
    ///
    /// Destruct the graph reader.
    ~BpUGraphReader() {
      if (own_reader) 
	delete reader;
    }

    /// \brief Give a new node map reading command to the reader.
    ///
    /// Give a new node map reading command to the reader.
    template <typename Map>
    BpUGraphReader& readNodeMap(std::string name, Map& map) {
      nodeset_reader.readNodeMap(name, map);
      return *this;
    }

    template <typename Map>
    BpUGraphReader& readNodeMap(std::string name, const Map& map) {
      nodeset_reader.readNodeMap(name, map);
      return *this;
    }

    /// \brief Give a new node map reading command to the reader.
    ///
    /// Give a new node map reading command to the reader.
    template <typename ItemReader, typename Map>
    BpUGraphReader& readNodeMap(std::string name, Map& map, 
                              const ItemReader& ir = ItemReader()) {
      nodeset_reader.readNodeMap(name, map, ir);
      return *this;
    }

    template <typename ItemReader, typename Map>
    BpUGraphReader& readNodeMap(std::string name, const Map& map, 
                              const ItemReader& ir = ItemReader()) {
      nodeset_reader.readNodeMap(name, map, ir);
      return *this;
    }

    /// \brief Give a new node map skipping command to the reader.
    ///
    /// Give a new node map skipping command to the reader.
    template <typename ItemReader>
    BpUGraphReader& skipNodeMap(std::string name, 
                              const ItemReader& ir = ItemReader()) {
      nodeset_reader.skipNodeMap(name, ir);
      return *this;
    }

    /// \brief Give a new A-node map reading command to the reader.
    ///
    /// Give a new A-node map reading command to the reader.
    template <typename Map>
    BpUGraphReader& readANodeMap(std::string name, Map& map) {
      nodeset_reader.readANodeMap(name, map);
      return *this;
    }

    template <typename Map>
    BpUGraphReader& readANodeMap(std::string name, const Map& map) {
      nodeset_reader.readANodeMap(name, map);
      return *this;
    }

    /// \brief Give a new A-node map reading command to the reader.
    ///
    /// Give a new A-node map reading command to the reader.
    template <typename ItemReader, typename Map>
    BpUGraphReader& readANodeMap(std::string name, Map& map, 
                              const ItemReader& ir = ItemReader()) {
      nodeset_reader.readANodeMap(name, map, ir);
      return *this;
    }

    template <typename ItemReader, typename Map>
    BpUGraphReader& readANodeMap(std::string name, const Map& map, 
                              const ItemReader& ir = ItemReader()) {
      nodeset_reader.readNodeMap(name, map, ir);
      return *this;
    }

    /// \brief Give a new A-node map skipping command to the reader.
    ///
    /// Give a new A-node map skipping command to the reader.
    template <typename ItemReader>
    BpUGraphReader& skipANodeMap(std::string name, 
                              const ItemReader& ir = ItemReader()) {
      nodeset_reader.skipANodeMap(name, ir);
      return *this;
    }

    /// \brief Give a new B-node map reading command to the reader.
    ///
    /// Give a new B-node map reading command to the reader.
    template <typename Map>
    BpUGraphReader& readBNodeMap(std::string name, Map& map) {
      nodeset_reader.readBNodeMap(name, map);
      return *this;
    }

    template <typename Map>
    BpUGraphReader& readBNodeMap(std::string name, const Map& map) {
      nodeset_reader.readBNodeMap(name, map);
      return *this;
    }

    /// \brief Give a new B-node map reading command to the reader.
    ///
    /// Give a new B-node map reading command to the reader.
    template <typename ItemReader, typename Map>
    BpUGraphReader& readBNodeMap(std::string name, Map& map, 
                              const ItemReader& ir = ItemReader()) {
      nodeset_reader.readBNodeMap(name, map, ir);
      return *this;
    }

    template <typename ItemReader, typename Map>
    BpUGraphReader& readBNodeMap(std::string name, const Map& map, 
                              const ItemReader& ir = ItemReader()) {
      nodeset_reader.readNodeMap(name, map, ir);
      return *this;
    }

    /// \brief Give a new B-node map skipping command to the reader.
    ///
    /// Give a new B-node map skipping command to the reader.
    template <typename ItemReader>
    BpUGraphReader& skipBNodeMap(std::string name, 
                              const ItemReader& ir = ItemReader()) {
      nodeset_reader.skipBNodeMap(name, ir);
      return *this;
    }

    /// \brief Give a new undirected edge map reading command to the reader.
    ///
    /// Give a new undirected edge map reading command to the reader.
    template <typename Map>
    BpUGraphReader& readUEdgeMap(std::string name, Map& map) { 
      uedgeset_reader.readUEdgeMap(name, map);
      return *this;
    }

    template <typename Map>
    BpUGraphReader& readUEdgeMap(std::string name, const Map& map) { 
      uedgeset_reader.readUEdgeMap(name, map);
      return *this;
    }


    /// \brief Give a new undirected edge map reading command to the reader.
    ///
    /// Give a new undirected edge map reading command to the reader.
    template <typename ItemReader, typename Map>
    BpUGraphReader& readUEdgeMap(std::string name, Map& map,
                               const ItemReader& ir = ItemReader()) {
      uedgeset_reader.readUEdgeMap(name, map, ir);
      return *this;
    }

    template <typename ItemReader, typename Map>
    BpUGraphReader& readUEdgeMap(std::string name, const Map& map,
                               const ItemReader& ir = ItemReader()) {
      uedgeset_reader.readUEdgeMap(name, map, ir);
      return *this;
    }

    /// \brief Give a new undirected edge map skipping command to the reader.
    ///
    /// Give a new undirected edge map skipping command to the reader.
    template <typename ItemReader>
    BpUGraphReader& skipUEdgeMap(std::string name,
				       const ItemReader& ir = ItemReader()) {
      uedgeset_reader.skipUMap(name, ir);
      return *this;
    }


    /// \brief Give a new edge map reading command to the reader.
    ///
    /// Give a new edge map reading command to the reader.
    template <typename Map>
    BpUGraphReader& readEdgeMap(std::string name, Map& map) { 
      uedgeset_reader.readEdgeMap(name, map);
      return *this;
    }

    template <typename Map>
    BpUGraphReader& readEdgeMap(std::string name, const Map& map) { 
      uedgeset_reader.readEdgeMap(name, map);
      return *this;
    }


    /// \brief Give a new edge map reading command to the reader.
    ///
    /// Give a new edge map reading command to the reader.
    template <typename ItemReader, typename Map>
    BpUGraphReader& readEdgeMap(std::string name, Map& map,
                              const ItemReader& ir = ItemReader()) {
      uedgeset_reader.readEdgeMap(name, map, ir);
      return *this;
    }

    template <typename ItemReader, typename Map>
    BpUGraphReader& readEdgeMap(std::string name, const Map& map,
                              const ItemReader& ir = ItemReader()) {
      uedgeset_reader.readEdgeMap(name, map, ir);
      return *this;
    }

    /// \brief Give a new edge map skipping command to the reader.
    ///
    /// Give a new edge map skipping command to the reader.
    template <typename ItemReader>
    BpUGraphReader& skipEdgeMap(std::string name,
                              const ItemReader& ir = ItemReader()) {
      uedgeset_reader.skipEdgeMap(name, ir);
      return *this;
    }

    /// \brief Give a new labeled node reading command to the reader.
    ///
    /// Give a new labeled node reading command to the reader.
    BpUGraphReader& readNode(std::string name, Node& node) {
      node_reader.readNode(name, node);
      return *this;
    }

    /// \brief Give a new labeled edge reading command to the reader.
    ///
    /// Give a new labeled edge reading command to the reader.
    BpUGraphReader& readEdge(std::string name, Edge& edge) {
      uedge_reader.readEdge(name, edge);
    }

    /// \brief Give a new labeled undirected edge reading command to the
    /// reader.
    ///
    /// Give a new labeled undirected edge reading command to the reader.
    BpUGraphReader& readUEdge(std::string name, UEdge& edge) {
      uedge_reader.readUEdge(name, edge);
    }

    /// \brief Give a new attribute reading command.
    ///
    ///  Give a new attribute reading command.
    template <typename Value>
    BpUGraphReader& readAttribute(std::string name, Value& value) {
      attribute_reader.readAttribute(name, value);
      return *this;
    }
    
    /// \brief Give a new attribute reading command.
    ///
    ///  Give a new attribute reading command.
    template <typename ItemReader, typename Value>
    BpUGraphReader& readAttribute(std::string name, Value& value, 
			       const ItemReader& ir = ItemReader()) {
      attribute_reader.readAttribute(name, value, ir);
      return *this;
    }

    /// \brief Conversion operator to LemonReader.
    ///
    /// Conversion operator to LemonReader. It make possible
    /// to access the encapsulated \e LemonReader, this way
    /// you can attach to this reader new instances of 
    /// \e LemonReader::SectionReader.
    operator LemonReader&() {
      return *reader;
    }

    /// \brief Executes the reading commands.
    ///
    /// Executes the reading commands.
    void run() {
      reader->run();
    }


    /// \brief Returns true if the reader can give back the items by its label.
    ///
    /// Returns true if the reader can give back the items by its label.
    bool isLabelReader() const {
      return nodeset_reader.isLabelReader() && 
        uedgeset_reader.isLabelReader();
    }

    /// \brief Gives back the node by its label.
    ///
    /// It reads an label from the stream and gives back which node belongs to
    /// it. It is possible only if there was read a "label" named node map.
    void readLabel(std::istream& is, Node& node) const {
      return nodeset_reader.readLabel(is, node);
    } 

    /// \brief Gives back the edge by its label
    ///
    /// It reads an label from the stream and gives back which edge belongs to
    /// it. It is possible only if there was read a "label" named edge map.
    void readLabel(std::istream& is, Edge& edge) const {
      return uedgeset_reader.readLabel(is, edge);
    } 

    /// \brief Gives back the undirected edge by its label.
    ///
    /// It reads an label from the stream and gives back which undirected edge 
    /// belongs to it. It is possible only if there was read a "label" named 
    /// edge map.
    void readLabel(std::istream& is, UEdge& uedge) const {
      return uedgeset_reader.readLabel(is, uedge);
    } 
    

  private:

    LemonReader* reader;
    bool own_reader;

    DefaultSkipper skipper;

    BpNodeSetReader<Graph, ReaderTraits> nodeset_reader;
    UEdgeSetReader<Graph, ReaderTraits> uedgeset_reader;

    NodeReader<Graph> node_reader;
    UEdgeReader<Graph> uedge_reader;
    
    AttributeReader<ReaderTraits> attribute_reader;
  };


  /// @}
}

#endif
