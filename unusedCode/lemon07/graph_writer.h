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
///\brief Lemon Graph Format writer.
///

#ifndef LEMON_GRAPH_WRITER_H
#define LEMON_GRAPH_WRITER_H

#include <iostream>

#include <lemon/error.h>
#include <lemon/lemon_writer.h>

namespace lemon {

  /// \addtogroup lemon_io
  /// @{

  /// \brief The graph writer class.
  ///
  /// The \c GraphWriter class provides the graph output.  Before you
  /// read this documentation it might be useful to read the general
  /// description of \ref graph-io-page "Graph Input-Output".
  ///
  /// To write a graph you should first give writing commands to the
  /// writer. You can declare write commands as \c NodeMap or \c
  /// EdgeMap writing and labeled Node and Edge writing.
  ///
  ///\code
  /// GraphWriter<ListGraph> writer(std::cout, graph);
  ///\endcode
  ///
  /// The \c writeNodeMap() function declares a \c NodeMap writing
  /// command in the \c GraphWriter. You should give as parameter the
  /// name of the map and the map object. The NodeMap writing command
  /// with name "label" should write a unique map because it is
  /// regarded as label map (such a map is essential if the graph has
  /// edges).
  ///
  ///\code
  /// IdMap<ListGraph, Node> nodeLabelMap;
  /// writer.writeNodeMap("label", nodeLabelMap);
  ///
  /// writer.writeNodeMap("coords", coords);
  /// writer.writeNodeMap("color", colorMap);
  ///\endcode
  ///
  /// With the \c writeEdgeMap() member function you can give an edge map
  /// writing command similar to the NodeMaps.
  ///
  ///\code
  /// DescriptorMap<ListGraph, Edge, ListGraph::EdgeMap<int> > 
  ///   edgeDescMap(graph);
  /// writer.writeEdgeMap("descriptor", edgeDescMap);
  ///
  /// writer.writeEdgeMap("weight", weightMap);
  /// writer.writeEdgeMap("label", labelMap);
  ///\endcode
  ///
  /// With \c writeNode() and \c writeEdge() functions you can 
  /// point out Nodes and Edges in the graph. For example, you can 
  /// write out the source and target of a maximum flow instance.
  ///
  ///\code
  /// writer.writeNode("source", sourceNode);
  /// writer.writeNode("target", targetNode);
  ///
  /// writer.writeEdge("observed", edge);
  ///\endcode
  ///
  /// After you give all write commands you must call the \c run() member
  /// function, which executes all the writing commands.
  ///
  ///\code
  /// writer.run();
  ///\endcode
  ///
  /// \see DefaultWriterTraits
  /// \see QuotedStringWriter
  /// \see IdMap
  /// \see DescriptorMap
  /// \see \ref GraphReader
  /// \see \ref graph-io-page
  /// \author Balazs Dezso
  template <typename _Graph, typename _WriterTraits = DefaultWriterTraits> 
  class GraphWriter {
  public:
    
    typedef _Graph Graph;
    typedef typename Graph::Node Node;
    typedef typename Graph::Edge Edge;

    typedef _WriterTraits WriterTraits;

    /// \brief Construct a new GraphWriter.
    ///
    /// This function constructs a new GraphWriter to write the given graph
    /// to the given stream.
    GraphWriter(std::ostream& _os, const Graph& _graph) 
      : writer(new LemonWriter(_os)), own_writer(true), 
	nodeset_writer(*writer, _graph, std::string()),
	edgeset_writer(*writer, _graph, nodeset_writer, std::string()),
	node_writer(*writer, nodeset_writer, std::string()),
	edge_writer(*writer, edgeset_writer, std::string()),
	attribute_writer(*writer, std::string()) {}

    /// \brief Construct a new GraphWriter.
    ///
    /// This function constructs a new GraphWriter to write the given graph
    /// to the given file.
    GraphWriter(const std::string& _filename, const Graph& _graph) 
      : writer(new LemonWriter(_filename)), own_writer(true), 
	nodeset_writer(*writer, _graph, std::string()),
	edgeset_writer(*writer, _graph, nodeset_writer, std::string()),
	node_writer(*writer, nodeset_writer, std::string()),
	edge_writer(*writer, edgeset_writer, std::string()),
	attribute_writer(*writer, std::string()) {}

    /// \brief Construct a new GraphWriter.
    ///
    /// This function constructs a new GraphWriter to write the given graph
    /// to the given LemonReader.
    GraphWriter(LemonWriter& _writer, const Graph& _graph)
      : writer(_writer), own_writer(false), 
	nodeset_writer(*writer, _graph, std::string()),
	edgeset_writer(*writer, _graph, nodeset_writer, std::string()),
	node_writer(*writer, nodeset_writer, std::string()),
	edge_writer(*writer, edgeset_writer, std::string()),
	attribute_writer(*writer, std::string()) {}

    /// \brief Destruct the graph writer.
    ///
    /// This function destructs the graph writer.
    ~GraphWriter() {
      if (own_writer) 
	delete writer;
    }

    /// \brief Issue a new node map writing command for the writer.
    ///
    /// This function issues a new <i> node map writing command</i> to the writer.
    template <typename Map>
    GraphWriter& writeNodeMap(std::string label, const Map& map) {
      nodeset_writer.writeNodeMap(label, map);
      return *this;
    }


    /// \brief Issue a new node map writing command for the writer.
    ///
    /// This function issues a new <i> node map writing command</i> to the writer.
    template <typename ItemWriter, typename Map>
    GraphWriter& writeNodeMap(std::string label, const Map& map, 
			      const ItemWriter& iw = ItemWriter()) {
      nodeset_writer.writeNodeMap(label, map, iw);
      return *this;
    }


    /// \brief Issue a new edge map writing command for the writer.
    ///
   /// This function issues a new <i> edge map writing command</i> to the writer.
    template <typename Map>
    GraphWriter& writeEdgeMap(std::string label, const Map& map) { 
      edgeset_writer.writeEdgeMap(label, map);
      return *this;
    }


    /// \brief Issue a new edge map writing command for the writer.
    ///
   /// This function issues a new <i> edge map writing command</i> to the writer.
    template <typename ItemWriter, typename Map>
    GraphWriter& writeEdgeMap(std::string label, const Map& map,
			      const ItemWriter& iw = ItemWriter()) {
      edgeset_writer.writeEdgeMap(label, map, iw);
      return *this;
    }

    /// \brief Issue a new labeled node writing command to the writer.
    ///
    /// This function issues a new <i> labeled node writing command</i> 
    /// to the writer.
    GraphWriter& writeNode(std::string label, const Node& node) {
      node_writer.writeNode(label, node);
      return *this;
    }

    /// \brief Issue a new labeled edge writing command to the writer.
    ///
    /// This function issues a new <i> labeled edge writing command</i> 
    /// to the writer.
    GraphWriter& writeEdge(std::string label, const Edge& edge) {
      edge_writer.writeEdge(label, edge);
    }

    /// \brief Issue a new attribute writing command.
    ///
    /// This function issues a new <i> attribute writing command</i> 
    /// to the writer.
    template <typename Value>
    GraphWriter& writeAttribute(std::string label, const Value& value) {
      attribute_writer.writeAttribute(label, value);
      return *this;
    }
    
    /// \brief Issue a new attribute writing command.
    ///
    /// This function issues a new <i> attribute writing command</i> 
    /// to the writer.
    template <typename ItemWriter, typename Value>
    GraphWriter& writeAttribute(std::string label, const Value& value, 
			       const ItemWriter& iw = ItemWriter()) {
      attribute_writer.writeAttribute(label, value, iw);
      return *this;
    }

    /// \brief Conversion operator to LemonWriter.
    ///
    /// Conversion operator to LemonWriter. It makes possible
    /// to access the encapsulated \e LemonWriter, this way
    /// you can attach to this writer new instances of 
    /// \e LemonWriter::SectionWriter. For more details see
    /// the \ref rwbackground "Background of Reading and Writing".
    operator LemonWriter&() {
      return *writer;
    }

    /// \brief Executes the writing commands.
    ///
    /// Executes the writing commands.
    void run() {
      writer->run();
    }

    /// \brief Returns true if the writer can give back the labels by the items.
    ///
    /// Returns true if the writer can give back the the labels by the items.
    bool isLabelWriter() const {
      return nodeset_writer.isLabelWriter() && 
        edgeset_writer.isLabelWriter();
    }

    /// \brief Write the label of the given node.
    ///
    /// It writes the label of the given node. If there was written a "label"
    /// named node map then it will write the map value belonging to the node.
    void writeLabel(std::ostream& os, const Node& item) const {
      nodeset_writer.writeLabel(os, item);
    } 

    /// \brief Write the label of the given edge.
    ///
    /// It writes the label of the given edge. If there was written a "label"
    /// named edge map then it will write the map value belonging to the edge.
    void writeLabel(std::ostream& os, const Edge& item) const {
      edgeset_writer.writeLabel(os, item);
    } 

    /// \brief Sorts the given node vector by label.
    ///
    /// Sorts the given node vector by label. If there was written an
    /// "label" named map then the vector will be sorted by the values
    /// of this map. Otherwise if the \c forceLabel parameter was true
    /// it will be sorted by its id in the graph.
    void sortByLabel(std::vector<Node>& nodes) const {
      nodeset_writer.sortByLabel(nodes);
    }

    /// \brief Sorts the given edge vector by label.
    ///
    /// Sorts the given edge vector by label. If there was written an
    /// "label" named map then the vector will be sorted by the values
    /// of this map. Otherwise if the \c forceLabel parameter was true
    /// it will be sorted by its id in the graph.
    void sortByLabel(std::vector<Edge>& edges) const {
      edgeset_writer.sortByLabel(edges);
    }

  private:

    LemonWriter* writer;
    bool own_writer;

    NodeSetWriter<Graph, WriterTraits> nodeset_writer;
    EdgeSetWriter<Graph, WriterTraits> edgeset_writer;

    NodeWriter<Graph> node_writer;
    EdgeWriter<Graph> edge_writer;
    
    AttributeWriter<WriterTraits> attribute_writer;
  };


  /// \brief The undirected graph writer class.
  ///
  /// The \c UGraphWriter class provides the ugraph output. To write 
  /// a graph you should first give writing commands to the writer. You can 
  /// declare write command as \c NodeMap, \c EdgeMap or \c UEdgeMap 
  /// writing and labeled Node, Edge or UEdge writing.
  ///
  ///\code
  /// UGraphWriter<ListUGraph> writer(std::cout, graph);
  ///\endcode
  ///
  /// The \c writeNodeMap() function declares a \c NodeMap writing 
  /// command in the \c UGraphWriter. You should give as parameter 
  /// the name of the map and the map object. The NodeMap writing 
  /// command with name "label" should write a unique map because it 
  /// is regarded as label map.
  ///
  ///\code
  /// IdMap<ListUGraph, Node> nodeLabelMap;
  /// writer.writeNodeMap("label", nodeLabelMap);
  ///
  /// writer.writeNodeMap("coords", coords);
  /// writer.writeNodeMap("color", colorMap);
  ///\endcode
  ///
  /// With the \c writeUEdgeMap() member function you can give an 
  /// undirected edge map writing command similar to the NodeMaps.
  ///
  ///\code
  /// DescriptorMap<ListGraph, Edge, ListGraph::EdgeMap<int> > 
  ///   edgeDescMap(graph);
  /// writer.writeUEdgeMap("descriptor", edgeDescMap);
  ///
  /// writer.writeUEdgeMap("weight", weightMap);
  /// writer.writeUEdgeMap("label", labelMap);
  ///\endcode
  /// 
  /// The EdgeMap handling is just a syntactical sugar. It writes
  /// two undirected edge map with '+' and '-' prefix in the name.
  ///
  ///\code
  /// writer.writeEdgeMap("capacity", capacityMap);
  ///\endcode
  ///
  ///
  /// With \c writeNode() and \c writeUEdge() functions you can 
  /// designate nodes and undirected edges in the graph. For example, you can 
  /// write out the source and target of the graph.
  ///
  ///\code
  /// writer.writeNode("source", sourceNode);
  /// writer.writeNode("target", targetNode);
  ///
  /// writer.writeUEdge("observed", uEdge);
  ///\endcode
  ///
  /// After you give all write commands you must call the \c run() member
  /// function, which executes all the writing commands.
  ///
  ///\code
  /// writer.run();
  ///\endcode
  ///
  /// \see DefaultWriterTraits
  /// \see QuotedStringWriter
  /// \see IdMap
  /// \see DescriptorMap
  /// \see \ref GraphWriter
  /// \see \ref graph-io-page
  /// \author Balazs Dezso
  template <typename _Graph, typename _WriterTraits = DefaultWriterTraits> 
  class UGraphWriter {
  public:
    
    typedef _Graph Graph;
    typedef typename Graph::Node Node;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::UEdge UEdge;

    typedef _WriterTraits WriterTraits;

    /// \brief Construct a new UGraphWriter.
    ///
    /// Construct a new UGraphWriter. It writes the given graph
    /// to the given stream.
    UGraphWriter(std::ostream& _os, const Graph& _graph) 
      : writer(new LemonWriter(_os)), own_writer(true), 
	nodeset_writer(*writer, _graph, std::string()),
	uedgeset_writer(*writer, _graph, nodeset_writer, std::string()),
	node_writer(*writer, nodeset_writer, std::string()),
	uedge_writer(*writer, uedgeset_writer, std::string()),
	attribute_writer(*writer, std::string()) {}

    /// \brief Construct a new UGraphWriter.
    ///
    /// Construct a new UGraphWriter. It writes the given graph
    /// to the given file.
    UGraphWriter(const std::string& _filename, const Graph& _graph) 
      : writer(new LemonWriter(_filename)), own_writer(true), 
	nodeset_writer(*writer, _graph, std::string()),
	uedgeset_writer(*writer, _graph, nodeset_writer, std::string()),
	node_writer(*writer, nodeset_writer, std::string()),
	uedge_writer(*writer, uedgeset_writer, std::string()),
	attribute_writer(*writer, std::string()) {}

    /// \brief Construct a new UGraphWriter.
    ///
    /// Construct a new UGraphWriter. It writes the given graph
    /// to given LemonWriter.
    UGraphWriter(LemonWriter& _writer, const Graph& _graph)
      : writer(_writer), own_writer(false), 
	nodeset_writer(*writer, _graph, std::string()),
	uedgeset_writer(*writer, _graph, nodeset_writer, std::string()),
	node_writer(*writer, nodeset_writer, std::string()),
	uedge_writer(*writer, uedgeset_writer, std::string()),
	attribute_writer(*writer, std::string()) {}

    /// \brief Destruct the graph writer.
    ///
    /// Destruct the graph writer.
    ~UGraphWriter() {
      if (own_writer) 
	delete writer;
    }

    /// \brief Issue a new node map writing command to the writer.
    ///
    /// This function issues a new <i> node map writing command</i> to
    /// the writer.
    template <typename Map>
    UGraphWriter& writeNodeMap(std::string label, const Map& map) {
      nodeset_writer.writeNodeMap(label, map);
      return *this;
    }

    /// \brief Issue a new node map writing command to the writer.
    ///
    /// This function issues a new <i> node map writing command</i> to
    /// the writer.
    template <typename ItemWriter, typename Map>
    UGraphWriter& writeNodeMap(std::string label, const Map& map, 
			      const ItemWriter& iw = ItemWriter()) {
      nodeset_writer.writeNodeMap(label, map, iw);
      return *this;
    }

    /// \brief Issue a new edge map writing command to the writer.
    ///
    /// This function issues a new <i> edge map writing command</i> to
    /// the writer.
    template <typename Map>
    UGraphWriter& writeEdgeMap(std::string label, const Map& map) { 
      uedgeset_writer.writeEdgeMap(label, map);
      return *this;
    }

    /// \brief Issue a new edge map writing command to the writer.
    ///
    /// This function issues a new <i> edge map writing command</i> to
    /// the writer.
    template <typename ItemWriter, typename Map>
    UGraphWriter& writeEdgeMap(std::string label, const Map& map,
				   const ItemWriter& iw = ItemWriter()) {
      uedgeset_writer.writeEdgeMap(label, map, iw);
      return *this;
    }

    /// \brief Issue a new undirected edge map writing command to the writer.
    ///
    /// This function issues a new <i> undirected edge map writing
    /// command</i> to the writer.
    template <typename Map>
    UGraphWriter& writeUEdgeMap(std::string label, const Map& map) { 
      uedgeset_writer.writeUEdgeMap(label, map);
      return *this;
    }

    /// \brief Issue a new undirected edge map writing command to the writer.
    ///
    /// This function issues a new <i> undirected edge map writing
    /// command</i> to the writer.
   template <typename ItemWriter, typename Map>
    UGraphWriter& writeUEdgeMap(std::string label, const Map& map,
					const ItemWriter& iw = ItemWriter()) {
      uedgeset_writer.writeUEdgeMap(label, map, iw);
      return *this;
    }

    /// \brief Issue a new labeled node writer to the writer.
    ///
    /// This function issues a new <i> labeled node writing
    /// command</i> to the writer.
    UGraphWriter& writeNode(std::string label, const Node& node) {
      node_writer.writeNode(label, node);
      return *this;
    }

    /// \brief Issue a new labeled edge writer to the writer.
    ///
    /// This function issues a new <i> labeled edge writing
    /// command</i> to the writer.
    UGraphWriter& writeEdge(std::string label, const Edge& edge) {
      uedge_writer.writeEdge(label, edge);
    }

    /// \brief Issue a new labeled undirected edge writing command to
    /// the writer.
    ///
    /// Issue a new <i>labeled undirected edge writing command</i> to
    /// the writer.
    UGraphWriter& writeUEdge(std::string label, const UEdge& edge) {
      uedge_writer.writeUEdge(label, edge);
    }

    /// \brief Issue a new attribute writing command.
    ///
    /// This function issues a new <i> attribute writing
    /// command</i> to the writer.
    template <typename Value>
    UGraphWriter& writeAttribute(std::string label, const Value& value) {
      attribute_writer.writeAttribute(label, value);
      return *this;
    }
    
    /// \brief Issue a new attribute writing command.
    ///
    /// This function issues a new <i> attribute writing
    /// command</i> to the writer.
    template <typename ItemWriter, typename Value>
    UGraphWriter& writeAttribute(std::string label, const Value& value, 
			       const ItemWriter& iw = ItemWriter()) {
      attribute_writer.writeAttribute(label, value, iw);
      return *this;
    }

    /// \brief Conversion operator to LemonWriter.
    ///
    /// Conversion operator to LemonWriter. It makes possible
    /// to access the encapsulated \e LemonWriter, this way
    /// you can attach to this writer new instances of 
    /// \e LemonWriter::SectionWriter.
    operator LemonWriter&() {
      return *writer;
    }

    /// \brief Executes the writing commands.
    ///
    /// Executes the writing commands.
    void run() {
      writer->run();
    }

    /// \brief Returns true if the writer can give back the labels by the items.
    ///
    /// Returns true if the writer can give back the the labels by the items.
    bool isLabelWriter() const {
      return nodeset_writer.isLabelWriter() && 
        uedgeset_writer.isLabelWriter();
    }

    /// \brief Write the label of the given node.
    ///
    /// It writes the label of the given node. If there was written a "label"
    /// named node map then it will write the map value belonging to the node.
    void writeLabel(std::ostream& os, const Node& item) const {
      nodeset_writer.writeLabel(os, item);
    } 

    /// \brief Write the label of the given edge.
    ///
    /// It writes the label of the given edge. If there was written a "label"
    /// named edge map then it will write the map value belonging to the edge.
    void writeLabel(std::ostream& os, const Edge& item) const {
      uedgeset_writer.writeLabel(os, item);
    } 

    /// \brief Write the label of the given undirected edge.
    ///
    /// It writes the label of the given undirected edge. If there was
    /// written a "label" named edge map then it will write the map
    /// value belonging to the edge.
    void writeLabel(std::ostream& os, const UEdge& item) const {
      uedgeset_writer.writeLabel(os, item);
    } 

    /// \brief Sorts the given node vector by label.
    ///
    /// Sorts the given node vector by label. If there was written an
    /// "label" named map then the vector will be sorted by the values
    /// of this map. Otherwise if the \c forceLabel parameter was true
    /// it will be sorted by its id in the graph.
    void sortByLabel(std::vector<Node>& nodes) const {
      nodeset_writer.sortByLabel(nodes);
    }

    /// \brief Sorts the given edge vector by label.
    ///
    /// Sorts the given edge vector by label. If there was written an
    /// "label" named map then the vector will be sorted by the values
    /// of this map. Otherwise if the \c forceLabel parameter was true
    /// it will be sorted by its id in the graph.
    void sortByLabel(std::vector<Edge>& edges) const {
      uedgeset_writer.sortByLabel(edges);
    }

    /// \brief Sorts the given undirected edge vector by label.
    ///
    /// Sorts the given undirected edge vector by label. If there was
    /// written an "label" named map then the vector will be sorted by
    /// the values of this map. Otherwise if the \c forceLabel
    /// parameter was true it will be sorted by its id in the graph.
    void sortByLabel(std::vector<UEdge>& uedges) const {
      uedgeset_writer.sortByLabel(uedges);
    }

  private:

    LemonWriter* writer;
    bool own_writer;

    NodeSetWriter<Graph, WriterTraits> nodeset_writer;
    UEdgeSetWriter<Graph, WriterTraits> uedgeset_writer;

    NodeWriter<Graph> node_writer;
    UEdgeWriter<Graph> uedge_writer;
    
    AttributeWriter<WriterTraits> attribute_writer;
  };

  /// \brief The bipartite graph writer class.
  ///
  /// The \c BpUGraphWriter class provides the ugraph output. To write 
  /// a graph you should first give writing commands to the writer. You can 
  /// declare write command as \c NodeMap, \c EdgeMap or \c UEdgeMap 
  /// writing and labeled Node, Edge or UEdge writing.
  ///
  ///\code
  /// BpUGraphWriter<ListUGraph> writer(std::cout, graph);
  ///\endcode
  ///
  /// The \c writeNodeMap() function declares a \c NodeMap writing 
  /// command in the \c BpUGraphWriter. You should give as parameter 
  /// the name of the map and the map object. The NodeMap writing 
  /// command with name "label" should write a unique map because it 
  /// is regarded as label map.
  ///
  ///\code
  /// IdMap<ListUGraph, Node> nodeLabelMap;
  /// writer.writeNodeMap("label", nodeLabelMap);
  ///
  /// writer.writeNodeMap("coords", coords);
  /// writer.writeNodeMap("color", colorMap);
  ///\endcode
  ///
  /// With the \c writeUEdgeMap() member function you can give an 
  /// undirected edge map writing command similar to the NodeMaps.
  ///
  ///\code
  /// DescriptorMap<ListGraph, Edge, ListGraph::EdgeMap<int> > 
  ///   edgeDescMap(graph);
  /// writer.writeUEdgeMap("descriptor", edgeDescMap);
  ///
  /// writer.writeUEdgeMap("weight", weightMap);
  /// writer.writeUEdgeMap("label", labelMap);
  ///\endcode
  /// 
  /// The EdgeMap handling is just a syntactical sugar. It writes
  /// two undirected edge map with '+' and '-' prefix in the name.
  ///
  ///\code
  /// writer.writeEdgeMap("capacity", capacityMap);
  ///\endcode
  ///
  ///
  /// With \c writeNode() and \c writeUEdge() functions you can 
  /// designate nodes and undirected edges in the graph. For example, you can 
  /// write out the source and target of the graph.
  ///
  ///\code
  /// writer.writeNode("source", sourceNode);
  /// writer.writeNode("target", targetNode);
  ///
  /// writer.writeUEdge("observed", uEdge);
  ///\endcode
  ///
  /// After you give all write commands you must call the \c run() member
  /// function, which executes all the writing commands.
  ///
  ///\code
  /// writer.run();
  ///\endcode
  ///
  /// \see DefaultWriterTraits
  /// \see QuotedStringWriter
  /// \see IdMap
  /// \see DescriptorMap
  /// \see \ref GraphWriter
  /// \see \ref graph-io-page
  /// \author Balazs Dezso
  template <typename _Graph, typename _WriterTraits = DefaultWriterTraits> 
  class BpUGraphWriter {
  public:
    
    typedef _Graph Graph;
    typedef typename Graph::Node Node;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::UEdge UEdge;

    typedef _WriterTraits WriterTraits;

    /// \brief Construct a new BpUGraphWriter.
    ///
    /// Construct a new BpUGraphWriter. It writes the given graph
    /// to the given stream.
    BpUGraphWriter(std::ostream& _os, const Graph& _graph) 
      : writer(new LemonWriter(_os)), own_writer(true), 
	nodeset_writer(*writer, _graph, std::string()),
	uedgeset_writer(*writer, _graph, nodeset_writer, std::string()),
	node_writer(*writer, nodeset_writer, std::string()),
	uedge_writer(*writer, uedgeset_writer, std::string()),
	attribute_writer(*writer, std::string()) {}

    /// \brief Construct a new BpUGraphWriter.
    ///
    /// Construct a new BpUGraphWriter. It writes the given graph
    /// to the given file.
    BpUGraphWriter(const std::string& _filename, const Graph& _graph) 
      : writer(new LemonWriter(_filename)), own_writer(true), 
	nodeset_writer(*writer, _graph, std::string()),
	uedgeset_writer(*writer, _graph, nodeset_writer, std::string()),
	node_writer(*writer, nodeset_writer, std::string()),
	uedge_writer(*writer, uedgeset_writer, std::string()),
	attribute_writer(*writer, std::string()) {}

    /// \brief Construct a new BpUGraphWriter.
    ///
    /// Construct a new BpUGraphWriter. It writes the given graph
    /// to given LemonWriter.
    BpUGraphWriter(LemonWriter& _writer, const Graph& _graph)
      : writer(_writer), own_writer(false), 
	nodeset_writer(*writer, _graph, std::string()),
	uedgeset_writer(*writer, _graph, nodeset_writer, std::string()),
	node_writer(*writer, nodeset_writer, std::string()),
	uedge_writer(*writer, uedgeset_writer, std::string()),
	attribute_writer(*writer, std::string()) {}

    /// \brief Destruct the graph writer.
    ///
    /// Destruct the graph writer.
    ~BpUGraphWriter() {
      if (own_writer) 
	delete writer;
    }

    /// \brief Issue a new node map writing command to the writer.
    ///
    /// This function issues a new <i> node map writing command</i> to
    /// the writer.
    template <typename Map>
    BpUGraphWriter& writeNodeMap(std::string label, const Map& map) {
      nodeset_writer.writeNodeMap(label, map);
      return *this;
    }

    /// \brief Issue a new node map writing command to the writer.
    ///
    /// This function issues a new <i> node map writing command</i> to
    /// the writer.
    template <typename ItemWriter, typename Map>
    BpUGraphWriter& writeNodeMap(std::string label, const Map& map, 
			      const ItemWriter& iw = ItemWriter()) {
      nodeset_writer.writeNodeMap(label, map, iw);
      return *this;
    }

    /// \brief Issue a new A-node map writing command to the writer.
    ///
    /// This function issues a new <i> A-node map writing command</i> to
    /// the writer.
    template <typename Map>
    BpUGraphWriter& writeANodeMap(std::string label, const Map& map) {
      nodeset_writer.writeANodeMap(label, map);
      return *this;
    }

    /// \brief Issue a new A-node map writing command to the writer.
    ///
    /// This function issues a new <i> A-node map writing command</i> to
    /// the writer.
    template <typename ItemWriter, typename Map>
    BpUGraphWriter& writeANodeMap(std::string label, const Map& map, 
			      const ItemWriter& iw = ItemWriter()) {
      nodeset_writer.writeANodeMap(label, map, iw);
      return *this;
    }
    /// \brief Issue a new B-node map writing command to the writer.
    ///
    /// This function issues a new <i> B-node map writing command</i> to
    /// the writer.
    template <typename Map>
    BpUGraphWriter& writeBNodeMap(std::string label, const Map& map) {
      nodeset_writer.writeBNodeMap(label, map);
      return *this;
    }

    /// \brief Issue a new B-node map writing command to the writer.
    ///
    /// This function issues a new <i> B-node map writing command</i> to
    /// the writer.
    template <typename ItemWriter, typename Map>
    BpUGraphWriter& writeBNodeMap(std::string label, const Map& map, 
			      const ItemWriter& iw = ItemWriter()) {
      nodeset_writer.writeBNodeMap(label, map, iw);
      return *this;
    }

    /// \brief Issue a new edge map writing command to the writer.
    ///
    /// This function issues a new <i> edge map writing command</i> to
    /// the writer.
    template <typename Map>
    BpUGraphWriter& writeEdgeMap(std::string label, const Map& map) { 
      uedgeset_writer.writeEdgeMap(label, map);
      return *this;
    }

    /// \brief Issue a new edge map writing command to the writer.
    ///
    /// This function issues a new <i> edge map writing command</i> to
    /// the writer.
    template <typename ItemWriter, typename Map>
    BpUGraphWriter& writeEdgeMap(std::string label, const Map& map,
				   const ItemWriter& iw = ItemWriter()) {
      uedgeset_writer.writeEdgeMap(label, map, iw);
      return *this;
    }

    /// \brief Issue a new undirected edge map writing command to the writer.
    ///
    /// This function issues a new <i> undirected edge map writing
    /// command</i> to the writer.
    template <typename Map>
    BpUGraphWriter& writeUEdgeMap(std::string label, const Map& map) { 
      uedgeset_writer.writeUEdgeMap(label, map);
      return *this;
    }

    /// \brief Issue a new undirected edge map writing command to the writer.
    ///
    /// This function issues a new <i> undirected edge map writing
    /// command</i> to the writer.
   template <typename ItemWriter, typename Map>
    BpUGraphWriter& writeUEdgeMap(std::string label, const Map& map,
					const ItemWriter& iw = ItemWriter()) {
      uedgeset_writer.writeUEdgeMap(label, map, iw);
      return *this;
    }

    /// \brief Issue a new labeled node writer to the writer.
    ///
    /// This function issues a new <i> labeled node writing
    /// command</i> to the writer.
    BpUGraphWriter& writeNode(std::string label, const Node& node) {
      node_writer.writeNode(label, node);
      return *this;
    }

    /// \brief Issue a new labeled edge writer to the writer.
    ///
    /// This function issues a new <i> labeled edge writing
    /// command</i> to the writer.
    BpUGraphWriter& writeEdge(std::string label, const Edge& edge) {
      uedge_writer.writeEdge(label, edge);
    }

    /// \brief Issue a new labeled undirected edge writing command to
    /// the writer.
    ///
    /// Issue a new <i>labeled undirected edge writing command</i> to
    /// the writer.
    BpUGraphWriter& writeUEdge(std::string label, const UEdge& edge) {
      uedge_writer.writeUEdge(label, edge);
    }

    /// \brief Issue a new attribute writing command.
    ///
    /// This function issues a new <i> attribute writing
    /// command</i> to the writer.
    template <typename Value>
    BpUGraphWriter& writeAttribute(std::string label, const Value& value) {
      attribute_writer.writeAttribute(label, value);
      return *this;
    }
    
    /// \brief Issue a new attribute writing command.
    ///
    /// This function issues a new <i> attribute writing
    /// command</i> to the writer.
    template <typename ItemWriter, typename Value>
    BpUGraphWriter& writeAttribute(std::string label, const Value& value, 
			       const ItemWriter& iw = ItemWriter()) {
      attribute_writer.writeAttribute(label, value, iw);
      return *this;
    }

    /// \brief Conversion operator to LemonWriter.
    ///
    /// Conversion operator to LemonWriter. It makes possible
    /// to access the encapsulated \e LemonWriter, this way
    /// you can attach to this writer new instances of 
    /// \e LemonWriter::SectionWriter.
    operator LemonWriter&() {
      return *writer;
    }

    /// \brief Executes the writing commands.
    ///
    /// Executes the writing commands.
    void run() {
      writer->run();
    }

    /// \brief Returns true if the writer can give back the labels by the items.
    ///
    /// Returns true if the writer can give back the the labels by the items.
    bool isLabelWriter() const {
      return nodeset_writer.isLabelWriter() && 
        uedgeset_writer.isLabelWriter();
    }

    /// \brief Write the label of the given node.
    ///
    /// It writes the label of the given node. If there was written a "label"
    /// named node map then it will write the map value belonging to the node.
    void writeLabel(std::ostream& os, const Node& item) const {
      nodeset_writer.writeLabel(os, item);
    } 

    /// \brief Write the label of the given edge.
    ///
    /// It writes the label of the given edge. If there was written a "label"
    /// named edge map then it will write the map value belonging to the edge.
    void writeLabel(std::ostream& os, const Edge& item) const {
      uedgeset_writer.writeLabel(os, item);
    } 

    /// \brief Write the label of the given undirected edge.
    ///
    /// It writes the label of the given undirected edge. If there was
    /// written a "label" named edge map then it will write the map
    /// value belonging to the edge.
    void writeLabel(std::ostream& os, const UEdge& item) const {
      uedgeset_writer.writeLabel(os, item);
    } 

    /// \brief Sorts the given node vector by label.
    ///
    /// Sorts the given node vector by label. If there was written an
    /// "label" named map then the vector will be sorted by the values
    /// of this map. Otherwise if the \c forceLabel parameter was true
    /// it will be sorted by its id in the graph.
    void sortByLabel(std::vector<Node>& nodes) const {
      nodeset_writer.sortByLabel(nodes);
    }

    /// \brief Sorts the given edge vector by label.
    ///
    /// Sorts the given edge vector by label. If there was written an
    /// "label" named map then the vector will be sorted by the values
    /// of this map. Otherwise if the \c forceLabel parameter was true
    /// it will be sorted by its id in the graph.
    void sortByLabel(std::vector<Edge>& edges) const {
      uedgeset_writer.sortByLabel(edges);
    }

    /// \brief Sorts the given undirected edge vector by label.
    ///
    /// Sorts the given undirected edge vector by label. If there was
    /// written an "label" named map then the vector will be sorted by
    /// the values of this map. Otherwise if the \c forceLabel
    /// parameter was true it will be sorted by its id in the graph.
    void sortByLabel(std::vector<UEdge>& uedges) const {
      uedgeset_writer.sortByLabel(uedges);
    }

  private:

    LemonWriter* writer;
    bool own_writer;

    BpNodeSetWriter<Graph, WriterTraits> nodeset_writer;
    UEdgeSetWriter<Graph, WriterTraits> uedgeset_writer;

    NodeWriter<Graph> node_writer;
    UEdgeWriter<Graph> uedge_writer;
    
    AttributeWriter<WriterTraits> attribute_writer;
  };

  /// @}

}

#endif
