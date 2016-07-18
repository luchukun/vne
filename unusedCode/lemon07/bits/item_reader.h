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
/// \brief Item reader bits for lemon input.

#ifndef LEMON_BITS_ITEM_READER_H
#define LEMON_BITS_ITEM_READER_H

#include <iostream>
#include <string>

#include <vector>
#include <deque>
#include <list>
#include <set>
#include <map>

#include <lemon/error.h>

namespace lemon {
  
  template <typename Value>
  class DefaultReader;

  /// \ingroup item_io
  ///
  /// \brief Reader class for unformatted strings.
  ///
  /// Reader class for unformatted strings. This class want to be
  /// a general reader type which can read the most 
  ///
  /// \author Balazs Dezso
  class UnformattedReader {
  public:
    /// \brief The value type of reader.
    ///
    /// The value type of reader.
    typedef std::string Value;
    
    /// \brief Constructor for the reader.
    ///
    /// Constructor for the reader.
    UnformattedReader() {} 
    
    /// \brief Reads an unformatted string from the given stream.
    ///
    /// Reads an unformatted string from the given stream.
    void read(std::istream& is, std::string& value) const {
      char c;
      value.clear();
      is >> std::ws;
      while (is.get(c) && !whiteSpace(c)) {
        processChar(c, is, value);
      }
    }

  private:

    void processChar(char c, std::istream& is, Value& value) const {
      switch (c) {
      case '(':
        is.putback(c);
        readParsed('(', ')', is, value);
        break;
      case '[':
        is.putback(c);
        readParsed('[', ']', is, value);
        break;
      case '{':
        is.putback(c);
        readParsed('{', '}', is, value);
        break;
      case '/':
        is.putback(c);
        readParsed('/', '/', is, value);
        break;
      case '\"':
        is.putback(c);
        readQuoted('\"', is, value);
        break;
      case '\'':
        is.putback(c);
        readQuoted('\'', is, value);
        break;
      default:
        value += c;
        break;
      }
    }

    void readParsed(char open, char close, 
                    std::istream& is, Value& value) const {
      char c;
      if (!is.get(c) || c != open)
	throw DataFormatError("Unformatted string format error");
      value += c;
      while (is.get(c) && c != close) {
        processChar(c, is, value);
      }
      if (!is) 
        throw DataFormatError("Unformatted string format error");
      value += c;      
    }

    void readQuoted(char quote, std::istream& is, Value& value) const {
      char c;
      bool esc = false;
      if (!is.get(c) || c != quote)
	throw DataFormatError("Unformatted string format error");
      value += c;
      while (is.get(c) && (c != quote || esc)) {
        if (c == '\\') esc = !esc;
        else esc = false;
        value += c;
      }
      if (!is) 
        throw DataFormatError("Unformatted string format error");
      value += c;
    }



    static bool whiteSpace(char c) {
      return c == ' ' || c == '\t' || c == '\v' || 
        c == '\n' || c == '\r' || c == '\f'; 
    }

    
  };

  /// \ingroup item_io
  ///
  /// \brief Reader class for quoted strings.
  ///
  /// Reader class for quoted strings. It can process the escape
  /// sequences in the string.
  ///
  /// \author Balazs Dezso
  class QuotedStringReader {
    friend class QuotedCharReader;
  public:
    /// \brief The value type of reader.
    ///
    /// The value type of reader.
    typedef std::string Value;
    
    /// \brief Constructor for the reader.
    ///
    /// Constructor for the reader. If the given parameter is true
    /// the reader processes the escape sequences.
    QuotedStringReader(bool _escaped = true) 
      : escaped(_escaped) {}
    
    /// \brief Reads a quoted string from the given stream.
    ///
    /// Reads a quoted string from the given stream.
    void read(std::istream& is, std::string& value) const {
      char c;
      value.clear();
      is >> std::ws;
      if (!is.get(c) || (c != '\"' && c != '\'')) 
	throw DataFormatError("Quoted format error");
      char quote = c;
      while (is.get(c) && c != quote) {
	if (escaped && c == '\\') {
	  value += readEscape(is);
	} else {
	  value += c;
	}
      }
      if (!is) throw DataFormatError("Quoted format error");
    }

  private:
    
    static char readEscape(std::istream& is) {
      char c;
      switch (is.get(c), c) {
      case '\\':
	return '\\';
      case '\"':
	return '\"';
      case '\'':
	return '\'';
      case '\?':
	return '\?';
      case 'a':
	return '\a';
      case 'b':
	return '\b';
      case 'f':
	return '\f';
      case 'n':
	return '\n';
      case 'r':
	return '\r';
      case 't':
	return '\t';
      case 'v':
	return '\v';
      case 'x':
	{
	  int code;
	  if (!is.get(c) || !isHex(c)) 
	    throw DataFormatError("Escape format error");
	  else if (code = valueHex(c), !is.get(c) || !isHex(c)) is.putback(c);
	  else code = code * 16 + valueHex(c);
	  return code;
	}
      default:
	{
	  int code;
	  if (!isOct(c)) 
	    throw DataFormatError("Escape format error");
	  else if (code = valueOct(c), !is.get(c) || !isOct(c)) 
	    is.putback(c);
	  else if (code = code * 8 + valueOct(c), !is.get(c) || !isOct(c)) 
	    is.putback(c);
	  else code = code * 8 + valueOct(c);
	  return code;
	}	      
      } 
    }

    static bool isOct(char c) {
      return '0' <= c && c <='7'; 
    }
    
    static int valueOct(char c) {
      return c - '0';
    }

   static bool isHex(char c) {
      return ('0' <= c && c <= '9') || 
	('a' <= c && c <= 'z') || 
	('A' <= c && c <= 'Z'); 
    }
    
    static int valueHex(char c) {
      if ('0' <= c && c <= '9') return c - '0';
      if ('a' <= c && c <= 'z') return c - 'a' + 10;
      return c - 'A' + 10;
    }

    bool escaped;
  };

  /// \ingroup item_io
  ///
  /// \brief Reader class for quoted strings.
  ///
  /// Reader class for quoted strings. It can process the escape
  /// sequences in the string.
  ///
  /// \author Balazs Dezso
  class QuotedCharReader {
  public:
    /// \brief The value type of reader.
    ///
    /// The value type of reader.
    typedef char Value;
    
    /// \brief Constructor for the reader.
    ///
    /// Constructor for the reader. If the given parameter is true
    /// the reader processes the escape sequences.
    QuotedCharReader(bool _escaped = true) 
      : escaped(_escaped) {}
    
    /// \brief Reads a quoted string from the given stream.
    ///
    /// Reads a quoted string from the given stream.
    void read(std::istream& is, char& value) const {
      char c;
      is >> std::ws;
      if (!is.get(c) || c != '\'') 
	throw DataFormatError("Quoted format error");
      if (!is.get(c)) 
        throw DataFormatError("Quoted format error");
      if (escaped && c == '\\') {
        value = QuotedStringReader::readEscape(is);
      } else {
        value = c;
      }
      if (!is.get(c) || c != '\'') 
	throw DataFormatError("Quoted format error");
    }

  private:
    bool escaped;
  };

  /// \ingroup item_io
  /// \brief Reader for standard containers.
  ///
  /// Reader for back insertable standard containers. The representation
  /// of the container is the values enumerated between an open and a
  /// close parse. 
  ///
  /// \author Balazs Dezso
  template <
    typename _Container, 
    typename _ItemReader = DefaultReader<typename _Container::value_type> 
  >
  class PushBackReader {
  public:
    typedef _Container Value;
    typedef _ItemReader ItemReader;

  private:

    ItemReader item_reader;

  public:

    /// \brief Constructor for InsertReader
    ///
    /// Constructor for InsertReader
    PushBackReader(const ItemReader& _item_reader = ItemReader())
      : item_reader(_item_reader) {}

    /// \brief Reads the values into the container from the given stream.
    ///
    /// Reads the values into the container from the given stream.
    void read(std::istream& is, Value& value) const {
      char c;
      if (!(is >> c) || c != '(') 
	throw DataFormatError("PushBackReader format error");
      while (is >> c && c != ')') {
	is.putback(c);
	typename ItemReader::Value item;
	item_reader.read(is, item);
	value.push_back(item);
      }
      if (!is) throw DataFormatError("PushBackReader format error");
    }

  };

  /// \ingroup item_io
  ///
  /// \brief Reader for standard containers.
  ///
  /// Reader for insertable standard containers. The representation
  /// of the container is the values enumerated between an open and a
  /// close parse. 
  ///
  /// \author Balazs Dezso
  template <
    typename _Container, 
    typename _ItemReader = DefaultReader<typename _Container::value_type> 
  >
  class InsertReader {
  public:
    typedef _Container Value;
    typedef _ItemReader ItemReader;

  private:

    ItemReader item_reader;

  public:

    /// \brief Constructor for InsertReader
    ///
    /// Constructor for InsertReader
    InsertReader(const ItemReader& _item_reader = ItemReader())
      : item_reader(_item_reader) {}

    /// \brief Reads the values into the container from the given stream.
    ///
    /// Reads the values into the container from the given stream.
    void read(std::istream& is, Value& value) const {
      char c;
      if (!(is >> c) || c != '(') 
	throw DataFormatError("InsertReader format error");
      while (is >> c && c != ')') {
	is.putback(c);
	typename ItemReader::Value item;
	item_reader.read(is, item);
	value.insert(item);
      }
      if (!is) throw DataFormatError("PushBackReader format error");
    }

  };

  /// \ingroup item_io
  /// \brief Reader for parsed string.
  ///
  /// Reader for parsed strings. You can define the open and close
  /// parse characters. It reads from the input a character sequence
  /// which is right parsed.
  ///
  /// \author Balazs Dezso
  class ParsedStringReader {
  public:
    typedef std::string Value;

    /// \brief Constructor.
    ///
    /// Constructor for ParsedStringReader. You can give as parameter
    /// the open and close parse characters.
    ParsedStringReader(char _open = '(', char _close = ')')
      : open(_open), close(_close) {}
    
    
    /// \brief Reads the parsed string from the given stream.
    ///
    /// Reads the parsed string from the given stream.
    void read(std::istream& is, Value& value) const {
      char c;
      value.clear();
      if (!(is >> c) || c != open) {
	throw DataFormatError("ParsedStringReader format error");
      }
      value += c;
      int counter = 1;
      while (counter > 0 && is >> c) {
	if (c == close) {
	  --counter;
	} else if (c == open) {
	  ++counter;
	}
	value += c;
      }
      if (!is) {
	throw DataFormatError("ParsedStrinReader format error");
      }
    }

  private:
    char open, close;

  };

  /// \ingroup item_io
  /// \brief Reader for read the whole line.
  ///
  /// Reader for read the whole line.
  ///
  /// \author Balazs Dezso
  class LineReader {
  public:
    typedef std::string Value;

    /// \brief Constructor.
    ///
    /// Constructor for the LineReader. If the given parameter is
    /// true then the spaces before the first not space character are
    /// skipped.
    LineReader(bool _skipSpaces = true) : skipSpaces(_skipSpaces) {}
    
    /// \brief Reads the line from the given stream.
    ///
    /// Reads the line from the given stream.
    void read(std::istream& is, Value& value) const {
      if (skipSpaces) is >> std::ws;
      if (!getline(is, value)) {
	throw DataFormatError("LineReader format error");
      }
    }
  private:
    bool skipSpaces;
  };

  /// \ingroup item_io
  /// \brief Reader for std::pair.
  ///
  /// Reader for std::pair.
  ///
  /// \author Balazs Dezso
  template <typename _Pair, 
	    typename _FirstReader = 
	    DefaultReader<typename _Pair::first_type>,
	    typename _SecondReader = 
	    DefaultReader<typename _Pair::second_type> >
  class PairReader {
  public:
    typedef _Pair Value;

    typedef _FirstReader FirstReader;
    typedef _SecondReader SecondReader;

  private:

    FirstReader first_reader;
    SecondReader second_reader;

  public:
    
    /// \brief Constructor.
    ///
    /// Constructor for the PairReader.
    PairReader(const FirstReader& _first_reader = FirstReader(), 
	       const SecondReader& _second_reader = SecondReader()) 
      : first_reader(_first_reader), second_reader(_second_reader) {}
    
    /// \brief Reads the pair from the given stream.
    ///
    /// Reads the pair from the given stream.
    void read(std::istream& is, Value& value) const {
      char c;
      if (!(is >> c) || c != '(') {
	throw DataFormatError("PairReader format error");
      }
      first_reader.read(is, value.first);
      if (!(is >> c) || c != '=') {
	throw DataFormatError("PairReader format error");
      }
      if (!(is >> c) || c != '>') {
	throw DataFormatError("PairReader format error");
      }
      second_reader.read(is, value.second);
      if (!(is >> c) || c != ')') {
	throw DataFormatError("PairReader format error");
      }
    }
  };

  /// \ingroup item_io
  /// 
  /// \brief The default item reader template class.
  ///
  /// The default item reader template class. If some section reader
  /// needs to read a value from a stream it will give the default way for it.
  ///
  /// \author Balazs Dezso
  template <typename _Value>
  class DefaultReader {
  public:
    /// The value type.
    typedef _Value Value;
    /// \brief Reads a value from the given stream.
    ///
    /// Reads a value from the given stream.
    void read(std::istream& is, Value& value) const {
      if (!(is >> value)) 
	throw DataFormatError("DefaultReader format error");
    }
  };

  template <>
  class DefaultReader<std::string> {
  public:
    typedef std::string Value;
    
    void read(std::istream& is, Value& value) const {
      char c;
      if (!(is >> std::ws >> c)) 
        throw DataFormatError("DefaultReader<string> format error");
      is.putback(c);
      switch (c) {
      case '\"':
	QuotedStringReader().read(is, value);
	break;
      default:
        UnformattedReader().read(is, value);
	break;
      }
    }
    
  };

  template <>
  class DefaultReader<char> {
  public:
    typedef char Value;
    
    void read(std::istream& is, Value& value) const {
      char c;
      if (!(is >> std::ws >> c))
        throw DataFormatError("DefaultReader<char> format error");
      is.putback(c);
      switch (c) {
      case '\'':
	QuotedCharReader().read(is, value);
	break;
      default:
        { 
          int temp;          
          if (!(is >> temp)) 
            throw DataFormatError("DefaultReader<char> format error");
          value = static_cast<char>(temp);
          break;
        }
      }
    }    
  };

  template <>
  class DefaultReader<bool> {
  public:
    typedef bool Value;
    
    void read(std::istream& is, Value& value) const {
      std::string rep;
      if (!(is >> rep))
        throw DataFormatError("DefaultReader<bool> format error");
      if (rep == "true" || rep == "t" || rep == "1") {
        value = true;
      } else if (rep == "false" || rep == "f" || rep == "0") {
        value = false;
      } else throw DataFormatError("DefaultReader<bool> format error");
    }    
  };

  template <typename Item>
  class DefaultReader<std::vector<Item> > 
    : public PushBackReader<std::vector<Item> > {};

  template <typename Item>
  class DefaultReader<std::deque<Item> > 
    : public PushBackReader<std::deque<Item> > {};

  template <typename Item>
  class DefaultReader<std::list<Item> > 
    : public PushBackReader<std::list<Item> > {};

  template <typename Item>
  class DefaultReader<std::set<Item> > 
    : public InsertReader<std::set<Item> > {};

  template <typename Key, typename Value>
  class DefaultReader<std::map<Key, Value> > 
    : public InsertReader<std::map<Key, Value>,
			  DefaultReader<std::pair<Key, Value> > > {};

  template <typename Item>
  class DefaultReader<std::multiset<Item> > 
    : public InsertReader<std::multiset<Item> > {};

  template <typename Key, typename Value>
  class DefaultReader<std::multimap<Key, Value> > 
    : public InsertReader<std::multimap<Key, Value>,
			  DefaultReader<std::pair<Key, Value> > > {};

  template <typename First, typename Second>
  class DefaultReader<std::pair<First, Second> > 
    : public PairReader<std::pair<First, Second> > {};

  /// \ingroup item_io
  /// 
  /// \brief The default item reader for skipping a value in the stream.
  ///
  /// The default item reader for skipping a value in the stream.
  ///
  /// \author Balazs Dezso
  class DefaultSkipper : public UnformattedReader {};

  /// \ingroup item_io  
  /// \brief Standard ReaderTraits for the GraphReader class.
  ///
  /// Standard ReaderTraits for the GraphReader class.
  /// It defines standard reading method for all type of value. 
  /// \author Balazs Dezso
  struct DefaultReaderTraits {

    template <typename _Value>
    struct Reader : DefaultReader<_Value> {};

    typedef DefaultSkipper Skipper;

  };

}

#endif
