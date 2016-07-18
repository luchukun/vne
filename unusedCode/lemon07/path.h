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

#ifndef LEMON_PATH_H
#define LEMON_PATH_H

#include <vector>
#include <algorithm>

#include <lemon/path_utils.h>
#include <lemon/error.h>
#include <lemon/bits/invalid.h>

namespace lemon {

  /// \addtogroup paths
  /// @{


  /// \brief A structure for representing directed paths in a graph.
  ///
  /// A structure for representing directed path in a graph.
  /// \param Graph The graph type in which the path is.
  ///
  /// In a sense, the path can be treated as a list of edges. The
  /// lemon path type stores just this list. As a consequence it
  /// cannot enumerate the nodes in the path and the zero length paths
  /// cannot store the source.
  ///
  /// This implementation is a back and front insertable and erasable
  /// path type. It can be indexed in O(1) time. The front and back
  /// insertion and erasure is amortized O(1) time. The
  /// impelementation is based on two opposite organized vectors.
  template <typename _Graph>
  class Path {
  public:

    typedef _Graph Graph;
    typedef typename Graph::Edge Edge;

    /// \brief Default constructor
    ///
    /// Default constructor
    Path() {}

    /// \brief Template copy constructor
    ///
    /// This path can be initialized with any other path type. It just
    /// makes a copy of the given path.
    template <typename CPath>
    Path(const CPath& cpath) {
      copyPath(*this, cpath);
    }

    /// \brief Template copy assignment
    ///
    /// This path can be initialized with any other path type. It just
    /// makes a copy of the given path.
    template <typename CPath>
    Path& operator=(const CPath& cpath) {
      copyPath(*this, cpath);
      return *this;
    }

    /// \brief Lemon style iterator for path edges
    ///
    /// This class is used to iterate on the edges of the paths.
    class EdgeIt {
      friend class Path;
    public:
      /// \brief Default constructor
      EdgeIt() {}
      /// \brief Invalid constructor
      EdgeIt(Invalid) : path(0), idx(-1) {}
      /// \brief Initializate the constructor to the first edge of path
      EdgeIt(const Path &_path) 
        : path(&_path), idx(_path.empty() ? -1 : 0) {}

    private:

      EdgeIt(const Path &_path, int _idx) 
        : path(&_path), idx(_idx) {}

    public:

      /// \brief Conversion to Edge
      operator const Edge&() const {
        return path->nth(idx);
      }

      /// \brief Next edge
      EdgeIt& operator++() { 
        ++idx;
        if (idx >= path->length()) idx = -1; 
        return *this; 
      }

      /// \brief Comparison operator
      bool operator==(const EdgeIt& e) const { return idx==e.idx; }
      /// \brief Comparison operator
      bool operator!=(const EdgeIt& e) const { return idx!=e.idx; }
      /// \brief Comparison operator
      bool operator<(const EdgeIt& e) const { return idx<e.idx; }

    private:
      const Path *path;
      int idx;
    };

    /// \brief Length of the path.
    int length() const { return head.size() + tail.size(); }
    /// \brief Returns whether the path is empty.
    bool empty() const { return head.empty() && tail.empty(); }

    /// \brief Resets the path to an empty path.
    void clear() { head.clear(); tail.clear(); }

    /// \brief Gives back the nth edge.
    ///
    /// \pre n is in the [0..length() - 1] range
    const Edge& nth(int n) const {
      return n < int(head.size()) ? *(head.rbegin() + n) :
        *(tail.begin() + (n - head.size()));
    }

    /// \brief Initializes edge iterator to point to the nth edge
    ///
    /// \pre n is in the [0..length() - 1] range
    EdgeIt nthIt(int n) const {
      return EdgeIt(*this, n);
    }

    /// \brief Gives back the first edge of the path
    const Edge& front() const {
      return head.empty() ? tail.front() : head.back();
    }

    /// \brief Add a new edge before the current path
    void addFront(const Edge& edge) {
      head.push_back(edge);
    }

    /// \brief Erase the first edge of the path
    void eraseFront() {
      if (!head.empty()) {
        head.pop_back();
      } else {
        head.clear();
        int halfsize = tail.size() / 2;
        head.resize(halfsize);
        std::copy(tail.begin() + 1, tail.begin() + halfsize + 1,
                  head.rbegin());
        std::copy(tail.begin() + halfsize + 1, tail.end(), tail.begin());
        tail.resize(tail.size() - halfsize - 1);
      }
    }

    /// \brief Gives back the last edge of the path
    const Edge& back() const {
      return tail.empty() ? head.front() : tail.back();
    }

    /// \brief Add a new edge behind the current path
    void addBack(const Edge& edge) {
      tail.push_back(edge);
    }

    /// \brief Erase the last edge of the path
    void eraseBack() {
      if (!tail.empty()) {
        tail.pop_back();
      } else {
        int halfsize = head.size() / 2;
        tail.resize(halfsize);
        std::copy(head.begin() + 1, head.begin() + halfsize + 1,
                  tail.rbegin());
        std::copy(head.begin() + halfsize + 1, head.end(), head.begin());
        head.resize(head.size() - halfsize - 1);
      }
    }



    typedef True BuildTag;

    template <typename CPath>
    void build(const CPath& path) {
      int len = path.length();
      tail.reserve(len);
      for (typename CPath::EdgeIt it(path); it != INVALID; ++it) {
        tail.push_back(it);
      }
    }

    template <typename CPath>
    void buildRev(const CPath& path) {
      int len = path.length();
      head.reserve(len);
      for (typename CPath::RevEdgeIt it(path); it != INVALID; ++it) {
        head.push_back(it);
      }
    }

  protected:
    typedef std::vector<Edge> Container;
    Container head, tail;

  };

  /// \brief A structure for representing directed paths in a graph.
  ///
  /// A structure for representing directed path in a graph.
  /// \param Graph The graph type in which the path is.
  ///
  /// In a sense, the path can be treated as a list of edges. The
  /// lemon path type stores just this list. As a consequence it
  /// cannot enumerate the nodes in the path and the zero length paths
  /// cannot store the source.
  ///
  /// This implementation is a just back insertable and erasable path
  /// type. It can be indexed in O(1) time. The back insertion and
  /// erasure is amortized O(1) time. This implementation is faster
  /// then the \c Path type because it use just one vector for the
  /// edges.
  template <typename _Graph>
  class SimplePath {
  public:

    typedef _Graph Graph;
    typedef typename Graph::Edge Edge;

    /// \brief Default constructor
    ///
    /// Default constructor
    SimplePath() {}

    /// \brief Template copy constructor
    ///
    /// This path can be initialized with any other path type. It just
    /// makes a copy of the given path.
    template <typename CPath>
    SimplePath(const CPath& cpath) {
      copyPath(*this, cpath);
    }

    /// \brief Template copy assignment
    ///
    /// This path can be initialized with any other path type. It just
    /// makes a copy of the given path.
    template <typename CPath>
    SimplePath& operator=(const CPath& cpath) {
      copyPath(*this, cpath);
      return *this;
    }

    /// \brief Iterator class to iterate on the edges of the paths
    ///
    /// This class is used to iterate on the edges of the paths
    ///
    /// Of course it converts to Graph::Edge
    class EdgeIt {
      friend class SimplePath;
    public:
      /// Default constructor
      EdgeIt() {}
      /// Invalid constructor
      EdgeIt(Invalid) : path(0), idx(-1) {}
      /// \brief Initializate the constructor to the first edge of path
      EdgeIt(const SimplePath &_path) 
        : path(&_path), idx(_path.empty() ? -1 : 0) {}

    private:

      /// Constructor with starting point
      EdgeIt(const SimplePath &_path, int _idx) 
        : idx(_idx), path(&_path) {}

    public:

      ///Conversion to Graph::Edge
      operator const Edge&() const {
        return path->nth(idx);
      }

      /// Next edge
      EdgeIt& operator++() { 
        ++idx;
        if (idx >= path->length()) idx = -1; 
        return *this; 
      }

      /// Comparison operator
      bool operator==(const EdgeIt& e) const { return idx==e.idx; }
      /// Comparison operator
      bool operator!=(const EdgeIt& e) const { return idx!=e.idx; }
      /// Comparison operator
      bool operator<(const EdgeIt& e) const { return idx<e.idx; }

    private:
      const SimplePath *path;
      int idx;
    };

    /// \brief Length of the path.
    int length() const { return data.size(); }
    /// \brief Returns whether the path is empty.
    bool empty() const { return data.empty(); }

    /// \brief Resets the path to an empty path.
    void clear() { data.clear(); }

    /// \brief Gives back the nth edge.
    ///
    /// \pre n is in the [0..length() - 1] range
    const Edge& nth(int n) const {
      return data[n];
    }

    /// \brief  Initializes edge iterator to point to the nth edge.
    EdgeIt nthIt(int n) const {
      return EdgeIt(*this, n);
    }

    /// \brief Gives back the first edge of the path.
    const Edge& front() const {
      return data.front();
    }

    /// \brief Gives back the last edge of the path.
    const Edge& back() const {
      return data.back();
    }

    /// \brief Add a new edge behind the current path.
    void addBack(const Edge& edge) {
      data.push_back(edge);
    }

    /// \brief Erase the last edge of the path
    void eraseBack() {
      data.pop_back();
    }

    typedef True BuildTag;

    template <typename CPath>
    void build(const CPath& path) {
      int len = path.length();
      data.resize(len);
      int index = 0;
      for (typename CPath::EdgeIt it(path); it != INVALID; ++it) {
        data[index] = it;;
        ++index;
      }
    }

    template <typename CPath>
    void buildRev(const CPath& path) {
      int len = path.length();
      data.resize(len);
      int index = len;
      for (typename CPath::RevEdgeIt it(path); it != INVALID; ++it) {
        --index;
        data[index] = it;;
      }
    }

  protected:
    typedef std::vector<Edge> Container;
    Container data;

  };

  /// \brief A structure for representing directed paths in a graph.
  ///
  /// A structure for representing directed path in a graph.
  /// \param Graph The graph type in which the path is.
  ///
  /// In a sense, the path can be treated as a list of edges. The
  /// lemon path type stores just this list. As a consequence it
  /// cannot enumerate the nodes in the path and the zero length paths
  /// cannot store the source.
  ///
  /// This implementation is a back and front insertable and erasable
  /// path type. It can be indexed in O(k) time, where k is the rank
  /// of the edge in the path. The length can be computed in O(n)
  /// time. The front and back insertion and erasure is O(1) time
  /// and it can be splited and spliced in O(1) time.
  template <typename _Graph>
  class ListPath {
  public:

    typedef _Graph Graph;
    typedef typename Graph::Edge Edge;

  protected:

    // the std::list<> is incompatible 
    // hard to create invalid iterator
    struct Node {
      Edge edge;
      Node *next, *prev;
    };

    Node *first, *last;

    std::allocator<Node> alloc;

  public:
 
    /// \brief Default constructor
    ///
    /// Default constructor
    ListPath() : first(0), last(0) {}

    /// \brief Template copy constructor
    ///
    /// This path can be initialized with any other path type. It just
    /// makes a copy of the given path.
    template <typename CPath>
    ListPath(const CPath& cpath) : first(0), last(0) {
      copyPath(*this, cpath);
    }

    /// \brief Destructor of the path
    ///
    /// Destructor of the path
    ~ListPath() {
      clear();
    }

    /// \brief Template copy assignment
    ///
    /// This path can be initialized with any other path type. It just
    /// makes a copy of the given path.
    template <typename CPath>
    ListPath& operator=(const CPath& cpath) {
      copyPath(*this, cpath);
      return *this;
    }

    /// \brief Iterator class to iterate on the edges of the paths
    ///
    /// This class is used to iterate on the edges of the paths
    ///
    /// Of course it converts to Graph::Edge
    class EdgeIt {
      friend class ListPath;
    public:
      /// Default constructor
      EdgeIt() {}
      /// Invalid constructor
      EdgeIt(Invalid) : path(0), node(0) {}
      /// \brief Initializate the constructor to the first edge of path
      EdgeIt(const ListPath &_path) 
        : path(&_path), node(_path.first) {}

    protected:

      EdgeIt(const ListPath &_path, Node *_node) 
        : path(&_path), node(_node) {}


    public:

      ///Conversion to Graph::Edge
      operator const Edge&() const {
        return node->edge;
      }

      /// Next edge
      EdgeIt& operator++() { 
        node = node->next;
        return *this; 
      }

      /// Comparison operator
      bool operator==(const EdgeIt& e) const { return node==e.node; }
      /// Comparison operator
      bool operator!=(const EdgeIt& e) const { return node!=e.node; }
      /// Comparison operator
      bool operator<(const EdgeIt& e) const { return node<e.node; }

    private:
      const ListPath *path;
      Node *node;
    };

    /// \brief Gives back the nth edge.
    ///
    /// Gives back the nth edge in O(n) time.
    /// \pre n is in the [0..length() - 1] range
    const Edge& nth(int n) const {
      Node *node = first;
      for (int i = 0; i < n; ++i) {
        node = node->next;
      }
      return node->edge;
    }

    /// \brief Initializes edge iterator to point to the nth edge.
    EdgeIt nthIt(int n) const {
      Node *node = first;
      for (int i = 0; i < n; ++i) {
        node = node->next;
      }
      return EdgeIt(*this, node);
    }

    /// \brief Length of the path.
    int length() const {
      int len = 0;
      Node *node = first;
      while (node != 0) {
        node = node->next;
        ++len;
      }
      return len;
    }

    /// \brief Returns whether the path is empty.
    bool empty() const { return first == 0; }

    /// \brief Resets the path to an empty path.
    void clear() {
      while (first != 0) {
        last = first->next;
        alloc.destroy(first);
        alloc.deallocate(first, 1);
        first = last;
      }
    }

    /// \brief Gives back the first edge of the path
    const Edge& front() const {
      return first->edge;
    }

    /// \brief Add a new edge before the current path
    void addFront(const Edge& edge) {
      Node *node = alloc.allocate(1);
      alloc.construct(node, Node());
      node->prev = 0;
      node->next = first;
      node->edge = edge;
      if (first) {
        first->prev = node;
        first = node;
      } else {
        first = last = node;
      }
    }

    /// \brief Erase the first edge of the path
    void eraseFront() {
      Node *node = first;
      first = first->next;
      if (first) {
        first->prev = 0;
      } else {
        last = 0;
      }
      alloc.destroy(node);
      alloc.deallocate(node, 1);
    }

    /// \brief Gives back the last edge of the path.
    const Edge& back() const {
      return last->edge;
    }

    /// \brief Add a new edge behind the current path.
    void addBack(const Edge& edge) {
      Node *node = alloc.allocate(1);
      alloc.construct(node, Node());
      node->next = 0;
      node->prev = last;
      node->edge = edge;
      if (last) {
        last->next = node;
        last = node;
      } else {
        last = first = node;
      }
    }

    /// \brief Erase the last edge of the path
    void eraseBack() {
      Node *node = last;
      last = last->prev;
      if (last) {
        last->next = 0;
      } else {
        first = 0;
      }
      alloc.destroy(node);
      alloc.deallocate(node, 1);
    }

    /// \brief Splicing the given path to the current path.
    ///
    /// It splices the \c tpath to the back of the current path and \c
    /// tpath becomes empty. The time complexity of this function is
    /// O(1).
    void spliceBack(ListPath& tpath) {
      if (first) {
        if (tpath.first) {
          last->next = tpath.first;
          tpath.first->prev = last;
          last = tpath.last;
        }
      } else {
        first = tpath.first;
        last = tpath.last;
      }
      tpath.first = tpath.last = 0;
    }

    /// \brief Splicing the given path to the current path.
    ///
    /// It splices the \c tpath before the current path and \c tpath
    /// becomes empty. The time complexity of this function
    /// is O(1).
    void spliceFront(ListPath& tpath) {
      if (first) {
        if (tpath.first) {
          first->prev = tpath.last;
          tpath.last->next = first;
          first = tpath.first;
        }
      } else {
        first = tpath.first;
        last = tpath.last;
      }
      tpath.first = tpath.last = 0;
    }

    /// \brief Splicing the given path into the current path.
    ///
    /// It splices the \c tpath into the current path before the
    /// position of \c it iterator and \c tpath becomes empty. The
    /// time complexity of this function is O(1). If the \c it is \c
    /// INVALID then it will splice behind the current path.
    void splice(EdgeIt it, ListPath& tpath) {
      if (it.node) {
        if (tpath.first) {
          tpath.first->prev = it.node->prev;
          if (it.node->prev) {
            it.node->prev->next = tpath.first;
          } else {
            first = tpath.first;
          }
          it.node->prev = tpath.last;
          tpath.last->next = it.node;
        }
      } else {
        if (first) {
          if (tpath.first) {
            last->next = tpath.first;
            tpath.first->prev = last;
            last = tpath.last;
          }
        } else {
          first = tpath.first;
          last = tpath.last;
        }
      }
      tpath.first = tpath.last = 0;
    }

    /// \brief Spliting the current path.
    ///
    /// It splits the current path into two parts. The part before \c
    /// it iterator will remain in the current path and the part from
    /// the it will put into the \c tpath. If the \c tpath had edges
    /// before the operation they will be removed first.  The time
    /// complexity of this function is O(1) plus the clearing of \c
    /// tpath. If the \c it is \c INVALID then it just clears \c
    /// tpath.
    void split(EdgeIt it, ListPath& tpath) {
      tpath.clear();
      if (it.node) {
        tpath.first = it.node;
        tpath.last = last;
        if (it.node->prev) {
          last = it.node->prev;
          last->next = 0;
        } else {
          first = last = 0;
        }
        it.node->prev = 0;
      }
    }


    typedef True BuildTag;

    template <typename CPath>
    void build(const CPath& path) {
      for (typename CPath::EdgeIt it(path); it != INVALID; ++it) {
        addBack(it);
      }
    }

    template <typename CPath>
    void buildRev(const CPath& path) {
      for (typename CPath::RevEdgeIt it(path); it != INVALID; ++it) {
        addFront(it);
      }
    }

  };

  /// \brief A structure for representing directed paths in a graph.
  ///
  /// A structure for representing directed path in a graph.
  /// \param Graph The graph type in which the path is.
  ///
  /// In a sense, the path can be treated as a list of edges. The
  /// lemon path type stores just this list. As a consequence it
  /// cannot enumerate the nodes in the path and the zero length paths
  /// cannot store the source.
  ///
  /// This implementation is completly static, so it cannot be
  /// modified exclude the assign an other path. It is intented to be
  /// used when you want to store a large number of paths because it is
  /// the most memory efficient path type in the lemon.
  template <typename _Graph>
  class StaticPath {
  public:

    typedef _Graph Graph;
    typedef typename Graph::Edge Edge;

    /// \brief Default constructor
    ///
    /// Default constructor
    StaticPath() : len(0), edges(0) {}
    
    /// \brief Template copy constructor
    ///
    /// This path can be initialized with any other path type. It just
    /// makes a copy of the given path.
    template <typename CPath>
    StaticPath(const CPath& cpath) : edges(0) {
      copyPath(*this, cpath);
    }

    /// \brief Destructor of the path
    ///
    /// Destructor of the path
    ~StaticPath() {
      if (edges) delete[] edges;
    }

    /// \brief Template copy assignment
    ///
    /// This path can be initialized with any other path type. It just
    /// makes a copy of the given path.
    template <typename CPath>
    StaticPath& operator=(const CPath& cpath) {
      copyPath(*this, cpath);
      return *this;
    }

    /// \brief Iterator class to iterate on the edges of the paths
    ///
    /// This class is used to iterate on the edges of the paths
    ///
    /// Of course it converts to Graph::Edge
    class EdgeIt {
      friend class StaticPath;
    public:
      /// Default constructor
      EdgeIt() {}
      /// Invalid constructor
      EdgeIt(Invalid) : path(0), idx(-1) {}
      /// Initializate the constructor to the first edge of path
      EdgeIt(const StaticPath &_path) 
        : path(&_path), idx(_path.empty() ? -1 : 0) {}

    private:

      /// Constructor with starting point
      EdgeIt(const StaticPath &_path, int _idx) 
        : idx(_idx), path(&_path) {}

    public:

      ///Conversion to Graph::Edge
      operator const Edge&() const {
        return path->nth(idx);
      }

      /// Next edge
      EdgeIt& operator++() { 
        ++idx;
        if (idx >= path->length()) idx = -1; 
        return *this; 
      }

      /// Comparison operator
      bool operator==(const EdgeIt& e) const { return idx==e.idx; }
      /// Comparison operator
      bool operator!=(const EdgeIt& e) const { return idx!=e.idx; }
      /// Comparison operator
      bool operator<(const EdgeIt& e) const { return idx<e.idx; }

    private:
      const StaticPath *path;
      int idx;
    };

    /// \brief Gives back the nth edge.
    ///
    /// \pre n is in the [0..length() - 1] range
    const Edge& nth(int n) const {
      return edges[n];
    }

    /// \brief Initializes edge iterator to point to the nth edge.
    EdgeIt nthIt(int n) const {
      return EdgeIt(*this, n);
    }

    /// \brief Gives back the length of the path.
    int length() const { return len; }

    /// \brief Returns true when the path is empty.
    int empty() const { return len == 0; }

    /// \break Erase all edge in the graph.
    void clear() {
      len = 0;
      if (edges) delete[] edges;
      edges = 0;
    }

    /// \brief Gives back the first edge of the path.
    const Edge& front() const {
      return edges[0];
    }

    /// \brief Gives back the last edge of the path.
    const Edge& back() const {
      return edges[len - 1];
    }


    typedef True BuildTag;

    template <typename CPath>
    void build(const CPath& path) {
      len = path.length();
      edges = new Edge[len];
      int index = 0;
      for (typename CPath::EdgeIt it(path); it != INVALID; ++it) {
        edges[index] = it;
        ++index;
      }
    }

    template <typename CPath>
    void buildRev(const CPath& path) {
      len = path.length();
      edges = new Edge[len];
      int index = len;
      for (typename CPath::RevEdgeIt it(path); it != INVALID; ++it) {
        --index;
        edges[index] = it;
      }
    }

  private:
    int len;
    Edge* edges;
  };

  ///@}

} // namespace lemon

#endif // LEMON_PATH_H
