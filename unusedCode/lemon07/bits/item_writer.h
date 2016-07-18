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

/// \ingroup item_io
/// \file
/// \brief Item writer bits for lemon output.

#ifndef LEMON_BITS_ITEM_WRITER_H
#define LEMON_BITS_ITEM_WRITER_H

#include <iostream>
#include <sstream>
#include <string>

#include <vector>
#include <deque>
#include <list>
#include <set>
#include <map>

#include <lemon/error.h>

namespace lemon {
  
  template <typename Value>
  class DefaultWriter;

  /// \ingroup item_io
  /// \brief Writer class for quoted strings.
  ///
  /// Writer class for unformatted strings.
  /// \author Balazs Dezso
  class UnformattedWriter {
  public:
    typedef std::string Value;

    /// \brief Constructor for the writer.
    ///
    /// Constructor for the writer.
    UnformattedWriter() {}

    /// \brief Writes an unformatted string to the given stream.
    ///
    /// Writes an unformatted string to the given stream.
    void write(std::ostream& os, const std::string& value) const {
      os << value;
    }

    bool readable(const std::string& value) const {
      std::istringstream is(value);
      char c;
      while (is.get(c) && !whiteSpace(c)) {
        if (!processChar(c, is)) return false;
      }
      if (is) return false;
      return true;
    }

  private:

    bool processChar(char c, std::istream& is) const {
      switch (c) {
      case '(':
        is.putback(c);
        if (!readableParsed('(', ')', is)) return false;
        break;
      case '[':
        is.putback(c);
        if (!readableParsed('[', ']', is)) return false;
        break;
      case '{':
        is.putback(c);
        if (!readableParsed('{', '}', is)) return false;
        break;
      case '/':
        is.putback(c);
        if (!readableParsed('/', '/', is)) return false;
        break;
      case '\"':
        is.putback(c);
        if (!readableQuoted('\"', is)) return false;
        break;
      case '\'':
        is.putback(c);
        if (!readableQuoted('\'', is)) return false;
        break;
      default:
        break;
      }
      return true;
    }

    bool readableParsed(char open, char close, std::istream& is) const {
      char c;
      if (!is.get(c) || c != open) return false;
      while (is.get(c) && c != close) {
        if (!processChar(c, is)) return false;
      }
      if (!is) return false;
      return true;
    }

    bool readableQuoted(char quote, std::istream& is) const {
      char c;
      bool esc = false;
      if (!is.get(c) || c != quote) return false;
      while (is.get(c) && (c != quote || esc)) {
        if (c == '\\') esc = !esc;
        else esc = false;
      }
      if (!is) return false;
      return true;
    }

    static bool whiteSpace(char c) {
      return c == ' ' || c == '\t' || c == '\v' || 
        c == '\n' || c == '\r' || c == '\f'; 
    }

  };

  /// \ingroup item_io
  /// \brief Writer class for quoted strings.
  ///
  /// Writer class for quoted strings. It can process the escape
  /// sequences in the string.
  /// \author Balazs Dezso
  class QuotedStringWriter {
    friend class QuotedCharWriter;
  public:
    typedef std::string Value;

    /// \brief Constructor for the writer.
    ///
    /// Constructor for the writer. If the given parameter is true
    /// the writer creates escape sequences from special characters.
    QuotedStringWriter(bool _escaped = true, char _quote = '\"') 
      : escaped(_escaped), quote(_quote) {}

    /// \brief Writes a quoted string to the given stream.
    ///
    /// Writes a quoted string to the given stream.
    void write(std::ostream& os, const std::string& value) const {
      os << quote;
      if (escaped) {
	std::ostringstream ls;
	for (int i = 0; i < int(value.size()); ++i) {
	  writeEscape(ls, value[i]);
	}
	os << ls.str();
      } else {
	os << value;
      }
      os << quote;
    }

  private:
    
    static void writeEscape(std::ostream& os, char c) {
      switch (c) {
      case '\\':
	os << "\\\\";
	return;
      case '\"':
	os << "\\\"";
	return;
      case '\'':
	os << "\\\'";
	return;
      case '\?':
	os << "\\\?";
	return;
      case '\a':
	os << "\\a";
	return;
      case '\b':
	os << "\\b";
	return;
      case '\f':
	os << "\\f";
	return;
      case '\r':
	os << "\\r";
	return;
      case '\n':
	os << "\\n";
	return;
      case '\t':
	os << "\\t";
	return;
      case '\v':
	os << "\\v";
	return;
      default:
	if (c < 0x20) {
	  os << '\\' << std::oct << static_cast<int>(c);
	} else {
	  os << c;
	}
	return;
      }     
    }
  private:
    bool escaped;
    char quote;
  };

  /// \ingroup item_io
  /// \brief Writer class for quoted chars.
  ///
  /// Writer class for quoted char. It can process the escape
  /// sequences in the string.
  /// \author Balazs Dezso
  class QuotedCharWriter {
  public:
    typedef char Value;

    /// \brief Constructor for the writer.
    ///
    /// Constructor for the writer. If the given parameter is true
    /// the writer creates escape sequences from special characters.
    QuotedCharWriter(bool _escaped = true) : escaped(_escaped) {}

    /// \brief Writes a quoted char to the given stream.
    ///
    /// Writes a quoted char to the given stream.
    void write(std::ostream& os, const char& value) const {
      os << "\'";
      if (escaped) {
	std::ostringstream ls;
        QuotedStringWriter::writeEscape(ls, value);
	os << ls.str();
      } else {
	os << value;
      }
      os << "\'";
    }

  private:
    bool escaped;
  };

  /// \ingroup item_io
  /// \brief Writer class for quoted char array.
  ///
  /// Writer class for quoted char array. It can process the escape
  /// sequences in the char array.
  /// \author Balazs Dezso
  class QuotedCharArrayWriter {
  public:
    typedef const char* Value;

    /// \brief Constructor for the writer.
    ///
    /// Constructor for the writer. If the given parameter is true
    /// the writer creates escape sequences from special characters.
    QuotedCharArrayWriter(bool _escaped = true) : escaped(_escaped) {}

    /// \brief Writes a quoted char array to the given stream.
    ///
    /// Writes a quoted char array to the given stream.
    void write(std::ostream& os, const char* value) const {
      QuotedStringWriter(escaped).write(os, std::string(value));
    }

  private:    
    bool escaped;
  };


  /// \ingroup item_io
  ///
  /// \brief Writer for standard containers.
  ///
  /// Writer for each iterable standard containers. The representation
  /// of the container is the values enumerated between an open and a
  /// close parse. 
  ///
  /// \author Balazs Dezso
  template <
    typename _Container, 
    typename _ItemWriter = DefaultWriter<typename _Container::value_type> 
  >
  class IterableWriter {
  public:
    typedef _Container Value;
    typedef _ItemWriter ItemWriter;

  private:

    ItemWriter item_writer;

  public:

    IterableWriter(const ItemWriter& _item_writer = ItemWriter())
      : item_writer(_item_writer) {}

    /// \brief Writes the values of the container to the given stream.
    ///
    /// Writes the values of the container to the given stream.
    void write(std::ostream& os, const Value& value) const {
      typename Value::const_iterator it;
      os << '(';
      for (it = value.begin(); it != value.end(); ++it) {
	item_writer.write(os, *it);
	os << ' ';
      }
      os << ')';
    }

  };

  /// \ingroup item_io
  ///
  /// \brief Writer for standard pairs.
  ///
  /// Writer for standard pairs. The representation of a pair is
  ///\code ( first_value => second_value ) \endcode.
  /// \author Balazs Dezso
  template <typename _Pair, 
	    typename _FirstWriter = 
	    DefaultWriter<typename _Pair::first_type>,
	    typename _SecondWriter = 
	    DefaultWriter<typename _Pair::second_type> >
  class PairWriter {
  public:

    typedef _Pair Value;

    typedef _FirstWriter FirstWriter;
    typedef _SecondWriter SecondWriter;

  private:

    FirstWriter first_writer;
    SecondWriter second_writer;

  public:
    
    /// \brief Constructor.
    ///
    /// Constructor for the PairWriter.
    PairWriter(const FirstWriter& _first_writer = FirstWriter(), 
	       const SecondWriter& _second_writer = SecondWriter()) 
      : first_writer(_first_writer), second_writer(_second_writer) {}
    
    /// \brief Writes the pair from the given stream.
    ///
    /// Writes the pair from the given stream.
    void write(std::ostream& os, const Value& value) const {
      os << "( ";
      first_writer.write(os, value.first);
      os << " => ";
      second_writer.write(os, value.second);
      os << " )";
    }

  };

  /// \ingroup item_io
  /// 
  /// \brief The default item writer template class.
  ///
  /// The default item writer template class. If some section writer
  /// needs to write a value to the stream it will give the default way for it.
  ///
  /// \author Balazs Dezso
  template <typename _Value>
  class DefaultWriter {
  public:
    /// The value type.
    typedef _Value Value;
    /// \brief Writes the value to the given stream.
    ///
    /// Writes the value to the given stream.
    void write(std::ostream& os, const Value& value) const {
      os << value;
    }
  };

  template <>
  class DefaultWriter<std::string> {
  public:
    typedef std::string Value;
    
    void write(std::ostream& os, const Value& value) const {
      if (UnformattedWriter().readable(value)) {
        UnformattedWriter().write(os, value);
      } else {
        QuotedStringWriter().write(os, value);
      }
    }
      
  };

  template <>
  class DefaultWriter<char> 
    : public QuotedCharWriter {};

  template <>
  class DefaultWriter<bool> {
  public:
    typedef bool Value;
    
    void write(std::ostream& os, const Value& value) const {
      os << (value ? "1" : "0");
    }
      
  };

  template <int length>
  class DefaultWriter<char[length]> 
    : public QuotedCharArrayWriter {};

  template <int length>
  class DefaultWriter<const char[length]> 
    : public QuotedCharArrayWriter {};

  template <>
  class DefaultWriter<char*> 
    : public QuotedCharArrayWriter {};

  template <>
  class DefaultWriter<const char*> 
    : public QuotedCharArrayWriter {};

  template <typename Item>
  class DefaultWriter<std::vector<Item> > 
    : public IterableWriter<std::vector<Item> > {};

  template <typename Item>
  class DefaultWriter<std::deque<Item> > 
    : public IterableWriter<std::deque<Item> > {};

  template <typename Item>
  class DefaultWriter<std::list<Item> > 
    : public IterableWriter<std::list<Item> > {};
  
  template <typename Item>
  class DefaultWriter<std::set<Item> > 
    : public IterableWriter<std::set<Item> > {};

  template <typename Key, typename Value>
  class DefaultWriter<std::map<Key, Value> > 
    : public IterableWriter<std::map<Key, Value> > {};

  template <typename Item>
  class DefaultWriter<std::multiset<Item> > 
    : public IterableWriter<std::multiset<Item> > {};

  template <typename Key, typename Value>
  class DefaultWriter<std::multimap<Key, Value> > 
    : public IterableWriter<std::multimap<Key, Value> > {};

  template <typename First, typename Second>
  class DefaultWriter<std::pair<First, Second> > 
    : public PairWriter<std::pair<First, Second> > {};

  /// \ingroup item_io
  /// \brief Standard WriterTraits for the section writers.
  ///
  /// Standard WriterTraits for the section writers.
  /// It defines standard writing method for all type of value. 
  /// \author Balazs Dezso
  struct DefaultWriterTraits {

    template <typename _Value>
    struct Writer : DefaultWriter<_Value> {
      typedef DefaultWriter<_Value> Parent;
    };

  };

}

#endif
