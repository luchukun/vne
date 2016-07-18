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

///\ingroup paths
///\file
///\brief Classes for representing paths in graphs.
///

#ifndef LEMON_PATH_UTILS_H
#define LEMON_PATH_UTILS_H

#include <lemon/concepts/path.h>
#include <lemon/lemon_reader.h>
#include <lemon/lemon_writer.h>

namespace lemon {

  namespace _path_bits {

    template <typename Path, typename Enable = void>
    struct RevTagIndicator {
      static const bool value = false;
    };

    template <typename Graph>
    struct RevTagIndicator<
      Graph, 
      typename enable_if<typename Graph::RevTag, void>::type
    > {
      static const bool value = true;
    };

    template <typename Target, typename Source,
              typename BuildEnable = void, typename RevEnable = void>
    struct PathCopySelector {
      static void copy(Target& target, const Source& source) {
        target.clear();
        for (typename Source::EdgeIt it(source); it != INVALID; ++it) {
          target.addBack(it);
        }
      }
    };

    template <typename Target, typename Source, typename BuildEnable>
    struct PathCopySelector<
      Target, Source, BuildEnable, 
      typename enable_if<typename Source::RevPathTag, void>::type> {
      static void copy(Target& target, const Source& source) {
        target.clear();
        for (typename Source::RevEdgeIt it(source); it != INVALID; ++it) {
          target.addFront(it);
        }
      }
    };

    template <typename Target, typename Source, typename RevEnable>
    struct PathCopySelector<
      Target, Source, 
      typename enable_if<typename Target::BuildTag, void>::type, RevEnable> {
      static void copy(Target& target, const Source& source) {
        target.clear();
        target.build(source);
      }
    };

    template <typename Target, typename Source>
    struct PathCopySelector<
      Target, Source, 
      typename enable_if<typename Target::BuildTag, void>::type,
      typename enable_if<typename Source::RevPathTag, void>::type> {
      static void copy(Target& target, const Source& source) {
        target.clear();
        target.buildRev(source);
      }
    };

  }


  /// \brief Make of copy of a path.
  ///
  ///  Make of copy of a path.
  template <typename Target, typename Source>
  void copyPath(Target& target, const Source& source) {
    checkConcept<concepts::PathDumper<typename Source::Graph>, Source>();
    _path_bits::PathCopySelector<Target, Source>::copy(target, source);
  }

  /// \brief Checks the path's consistency.
  ///
  /// Checks that each edge's target is the next's source. 
  /// 
  template <typename Graph, typename Path>
  bool checkPath(const Graph& graph, const Path& path) {
    typename Path::EdgeIt it(path);
    if (it == INVALID) return true;
    typename Graph::Node node = graph.target(it);
    ++it;
    while (it != INVALID) {
      if (graph.source(it) != node) return false;
      node = graph.target(it);
      ++it;
    }
    return true;
  }

  /// \brief Gives back the source of the path
  ///
  /// Gives back the source of the path.
  template <typename Graph, typename Path>
  typename Graph::Node pathSource(const Graph& graph, const Path& path) {
    return graph.source(path.front());
  }

  /// \brief Gives back the target of the path
  ///
  /// Gives back the target of the path.
  template <typename Graph, typename Path>
  typename Graph::Node pathTarget(const Graph& graph, const Path& path) {
    return graph.target(path.back());
  }

  /// \brief Class which helps to iterate the nodes of a path
  ///
  /// In a sense, the path can be treated as a list of edges. The
  /// lemon path type stores just this list. As a consequence it
  /// cannot enumerate the nodes in the path and the zero length paths
  /// cannot store the node.
  ///
  /// This class implements the node iterator of a path structure. To
  /// provide this feature, the underlying graph should be given to
  /// the constructor of the iterator.
  template <typename Path>
  class PathNodeIt {
  private:
    const typename Path::Graph *_graph;
    typename Path::EdgeIt _it;
    typename Path::Graph::Node _nd;

  public:

    typedef typename Path::Graph Graph;
    typedef typename Graph::Node Node;
    
    /// Default constructor
    PathNodeIt() {}
    /// Invalid constructor
    PathNodeIt(Invalid) 
      : _graph(0), _it(INVALID), _nd(INVALID) {}
    /// Constructor
    PathNodeIt(const Graph& graph, const Path& path) 
      : _graph(&graph), _it(path) {
      _nd = (_it != INVALID ? _graph->source(_it) : INVALID);
    }
    /// Constructor
    PathNodeIt(const Graph& graph, const Path& path, const Node& src) 
      : _graph(&graph), _it(path), _nd(src) {}

    ///Conversion to Graph::Node
    operator Node() const {
      return _nd;
    }

    /// Next node
    PathNodeIt& operator++() {
      if (_it == INVALID) _nd = INVALID;
      else {
	_nd = _graph->target(_it);
	++_it;
      }
      return *this;
    }

    /// Comparison operator
    bool operator==(const PathNodeIt& n) const { 
      return _it == n._it && _nd == n._nd; 
    }
    /// Comparison operator
    bool operator!=(const PathNodeIt& n) const { 
      return _it != n._it || _nd != n._nd; 
    }
    /// Comparison operator
    bool operator<(const PathNodeIt& n) const { 
      return (_it < n._it && _nd != INVALID);
    }
    
  };

  /// \brief Item writer for paths
  ///
  /// This class can write paths into files. You can store paths in
  /// distinict mapset or in attributes section.
  ///
  ///\code
  /// GraphWriter<SmartGraph> gw(std::cout, g);
  /// NodeMapWriter<SmartGraph> nmw(gw, g, gw);
  ///
  /// SmartGraph::NodeMap<Path<SmartGraph> > pnm(g);
  /// for (SmartGraph::NodeIt n(g); n != INVALID; ++n) {
  ///   pnm[n] = bfs.path(n);
  /// }
  /// nmw.writeNodeMap("pnm", pnm, PathWriter<Path<SmartGraph> >(gw));
  ///
  /// gw.run();
  ///\endcode
  ///
  /// \warning Do not use this class to write node or edge map values
  /// into usual nodesets or edgesets. You will not be able to read
  /// back your paths. Rather use NodeMapWriter, EdgeSetWriter or
  /// UEdgeSetWriter to dump paths from maps to lemon file.
  template <typename Path>
  class PathWriter {
  private:

    typedef typename Path::Edge Edge;
    std::auto_ptr<_writer_bits::LabelWriterBase<Edge> > edgeLabelWriter;

  public:

    typedef Path Value;

    PathWriter(const PathWriter& pw) {
      edgeLabelWriter.reset(pw.edgeLabelWriter->clone());
    }

    /// \brief Constructor
    ///
    /// The paramter shold be an edge label writer which could
    /// be a GraphWriter or an EdgeSetWriter. 
    template <typename EdgeLabelWriter>
    explicit PathWriter(const EdgeLabelWriter& _edgeLabelWriter) {
      edgeLabelWriter.reset(new _writer_bits::
	LabelWriter<Edge, EdgeLabelWriter>(_edgeLabelWriter));
    }

    /// \brief Writer function
    ///
    /// Writes the path to the current stream. The representation
    /// is the edge labels beetween parentheses.
    void write(std::ostream& os, const Value& value) const {
      if (!edgeLabelWriter->isLabelWriter()) {
	throw DataFormatError("Cannot find edgeset or label map");
      }
      os << '(' << ' ';
      for (typename Path::EdgeIt e(value); e != INVALID; ++e) {
	edgeLabelWriter->write(os, e);
	os << ' ';
      }
      os << ')';
    }
    
  };

  namespace _path_bits {

    template <typename _Graph>
    class PathProxy {
    public:
      typedef False RevPathTag;

      typedef _Graph Graph;
      typedef typename Graph::Edge Edge;

      PathProxy(const std::vector<Edge>& edges)
	: _edges(edges) {}

      int length() const {
	return _edges.size();
      }

      bool empty() const {
	return _edges.size() == 0;
      }

      class EdgeIt {
      public:
	EdgeIt() {}
	EdgeIt(const PathProxy& path) 
	  : _path(&path), _index(0) {}
	
	operator const Edge() const {
	  return _path->_edges[_index];
	}

	EdgeIt& operator++() {
	  ++_index;
	  return *this;
	}

	bool operator==(Invalid) const { 
	  return int(_path->_edges.size()) == _index; 
	}
	bool operator!=(Invalid) const { 
	  return int(_path->_edges.size()) != _index; 
	}

      private:
	const PathProxy* _path;
	int _index;
      };
      
    private:
      const std::vector<Edge>& _edges;
      
    };

  }

  /// \brief Item reader for paths 
  ///
  /// This class can read paths from files. You can store paths in
  /// distinict mapset or in attributes section.
  ///
  ///\code
  /// GraphReader<SmartGraph> gr(std::cout, g);
  /// NodeMapReader<SmartGraph> nmr(gr, g, gr);
  ///
  /// SmartGraph::NodeMap<Path<SmartGraph> > pnm(g);
  /// nmr.readNodeMap("pnm", pnm, PathReader<Path<SmartGraph> >(gr));
  ///
  /// gr.run();
  ///\endcode
  ///
  /// \warning Do not use this class to read node or edge map values
  /// from nodesets or edgesets. The edges are not surely constructed
  /// when the edge list should be read. Rather use NodeMapReader,
  /// EdgeSetReader or UEdgeSetReader to read distinict map sets from file.
  template <typename Path>
  class PathReader {
  private:

    typedef typename Path::Edge Edge;
    std::auto_ptr<_reader_bits::LabelReaderBase<Edge> > edgeLabelReader;

  public:

    typedef Path Value;

    PathReader(const PathReader& pw) {
      edgeLabelReader.reset(pw.edgeLabelReader->clone());
    }

    /// \brief Constructor
    ///
    /// The paramter shold be an edge label reader which could
    /// be a GraphReader or an EdgeSetReader. 
    template <typename EdgeLabelReader>
    explicit PathReader(const EdgeLabelReader& _edgeLabelReader) {
      edgeLabelReader.reset(new _reader_bits::
	LabelReader<Edge, EdgeLabelReader>(_edgeLabelReader));
    }


    /// \brief Reader function
    ///
    /// Reads the path from the current stream. The representation
    /// is the edge labels beetween parentheses.
    void read(std::istream& is, Value& value) const {
      if (!edgeLabelReader->isLabelReader()) {
	throw DataFormatError("Cannot find edgeset or label map");
      }
      char c;
      if (!(is >> c) || c != '(') 
	throw DataFormatError("PathReader format error");
      std::vector<typename Path::Edge> v;
      while (is >> c && c != ')') {
	is.putback(c);
	Edge edge = edgeLabelReader->read(is);
	v.push_back(edge);
      }
      if (!is) throw DataFormatError("PathReader format error");
      copyPath(value, _path_bits::PathProxy<typename Path::Edge>(v));
    }
    
  };
  
}

#endif
