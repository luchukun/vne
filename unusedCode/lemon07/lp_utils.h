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

#ifndef LEMON_LP_UTILS_H
#define LEMON_LP_UTILS_H

#include <lemon/lp_base.h>

#include <lemon/lemon_reader.h>
#include <lemon/lemon_writer.h>

#include <map>
#include <set>


namespace lemon {

  /// \ingroup lp_utils
  ///
  /// \brief The map for the result of the lp variables.
  ///
  /// The map for the result of the lp variables.
  class LpResultMap {
  public:
    typedef LpSolverBase::Col Key;
    typedef LpSolverBase::Value Value;

    /// \brief Constructor
    ///
    /// Constructor
    LpResultMap(const LpSolverBase& _lp) : lp(_lp) {}
    LpSolverBase::Value operator[](const LpSolverBase::Col col) const {
      return lp.primal(col);
    }
  private:
    const LpSolverBase& lp;
  };

  /// \ingroup lp_utils
  ///
  /// \brief Returns an \ref LpResultMap class
  ///
  /// This function just returns an \ref LpResultMap class.
  /// \relates LpResultMap
  LpResultMap lpResultMap(const LpSolverBase& lp);

  /// \ingroup lp_utils
  ///
  /// \brief The map for the name of the lp variables.
  ///
  /// The map for the name of the lp variables.
  class LpColNameMap {
  public:
    typedef LpSolverBase::Col Key;
    typedef std::string Value;

    /// \brief Constructor
    ///
    /// Constructor
    LpColNameMap(const LpSolverBase& _lp) : lp(_lp) {}
    std::string operator[](const LpSolverBase::Col col) const {
      return lp.colName(col);
    }
  private:
    const LpSolverBase& lp;
  };

  /// \ingroup lp_utils
  ///
  /// \brief Writeable map for the name of the lp variables.
  ///
  /// Writeable map for the name of the lp variables.
  class LpColNameWriteMap {
  public:
    typedef LpSolverBase::Col Key;
    typedef std::string Value;

    /// \brief Constructor
    ///
    /// Constructor
    LpColNameWriteMap(LpSolverBase& _lp) : lp(_lp) {}
    std::string operator[](const LpSolverBase::Col col) const {
      return lp.colName(col);
    }
    void set(const LpSolverBase::Col col, const std::string& name) {
      lp.colName(col, name);
    }
  private:
    LpSolverBase& lp;
  };

  /// \ingroup lp_utils
  ///
  /// \brief Returns an \ref LpColNameMap class
  ///
  /// This function just returns an \ref LpColNameMap class.
  /// \relates LpColNameMap
  LpColNameMap lpColNameMap(const LpSolverBase& lp);

  LpColNameWriteMap lpColNameMap(LpSolverBase& lp);

  /// \ingroup lp_utils
  ///
  /// \brief Lp variable item reader for Lemon IO
  ///
  /// This class is an Lp variable item reader for Lemon IO. It makes
  /// possible to store lp variables in lemon file. The usage of this
  /// class is very simple:
  ///
  ///\code
  /// Graph graph;
  /// Lp lp;
  /// Graph::EdgeMap<Lp::Col> var(graph); 
  ///
  /// GraphReader<Graph> reader(cin, graph);
  /// reader.readEdgeMap("lpvar", var, LpColReader(lp));
  /// reader.run();
  ///\endcode
  ///
  /// If there is already a variable with the same name in the lp then
  /// it will be used for the return value. If there is no such variable
  /// then it will be constructed. The file could contain \c '-' as value
  /// which means that no lp variable associated to the current item and
  /// in this case INVALID will be returned.
  class LpColReader {
  public:

    /// \brief The value type of reader.
    ///
    /// The value type of reader.
    typedef LpSolverBase::Col Value;

    /// \brief Constructor for the reader.
    ///
    /// Constructor for the reader.
    LpColReader(LpSolverBase& _lp) : lp(_lp) {}

    /// \brief Reads an lp variable from the given stream.
    ///
    /// Reads an lp variable string from the given stream.
    void read(std::istream& is, LpSolverBase::Col& col) const {
      char x;
      is >> std::ws;
      if (!is.get(x))
        throw DataFormatError("Wrong Lp Column format!");
      if (varFirstChar(x)) {
        std::string var;
        var += x;
        
        while (is.get(x) && varChar(x)) {
          var += x;
        }
        if (!is) {
          is.clear();
        } else {
          is.putback(x);
        }
        col = lp.colByName(var);
        if (col == INVALID) {
          col = lp.addCol();
          lp.colName(col, var);
        }
      } else if (x == '-') {
        col = INVALID;
        is >> std::ws;
      } else {
        std::cerr << x << std::endl;
        throw DataFormatError("Wrong Lp Column format");
      }
    }

  private:

    static bool varFirstChar(char c) {
      return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
    }

    static bool varChar(char c) {
      return (c >= '0' && c <= '9') || 
        (c >= 'a' && c <= 'z') || 
        (c >= 'A' && c <= 'Z') ||
        c == '[' || c == ']';
    }
    
    LpSolverBase& lp;

  };

  /// \ingroup lp_utils
  ///
  /// \brief Lp variable item writer for Lemon IO
  ///
  /// This class is an Lp variable item writer for Lemon IO. It makes
  /// possible to store lp variables in lemon file. The usage of this
  /// class is very simple:
  ///
  ///\code
  /// Graph graph;
  /// Lp lp;
  /// Graph::EdgeMap<Lp::Col> var(graph); 
  ///
  /// GraphWriter<Graph> writer(cin, graph);
  /// writer.writeEdgeMap("lpvar", var, LpColWriter(lp));
  /// writer.run();
  ///\endcode
  ///
  /// If there is no name associated to the current item then the name
  /// will automatically constructed. If the value is INVALID then it
  /// will write an \c '-' value to the file.
  class LpColWriter {
  public:

    /// \brief The value type of writer.
    ///
    /// The value type of writer.
    typedef LpSolverBase::Col Value;

    /// \brief Constructor for the writer.
    ///
    /// Constructor for the writer.
    LpColWriter(const LpSolverBase& _lp) : lp(_lp), num(0) {}

    /// \brief Writes an lp variable to the given stream.
    ///
    /// Writes an lp variable to the given stream.
    void write(std::ostream& os, const LpSolverBase::Col& col) const {
      if (col != INVALID) {
        std::string name = lp.colName(col);
        if (name.empty()) {
          std::ostringstream ls;
          ls << "x" << num;
          ++num;
          while (lp.colByName(ls.str()) != INVALID) {
            ls.str(std::string());
            ls << "x" << num;
            ++num;
          }
          const_cast<LpSolverBase&>(lp).colName(col, ls.str());
          os << ls.str();
        } else {
          os << name;
        }
      } else {
        os << "-";
      }
    }

  private:

    const LpSolverBase& lp;
    mutable int num;

  };

  /// \ingroup lp_utils
  ///
  /// \brief Lp section reader for lemon IO
  ///
  /// This section reader provides an easy way to read an Lp problem
  /// from a file. The lemon lp section format contains three parts.
  ///
  ///\code
  /// @lp
  /// constraints
  ///   7 == x1 - 1 * x2
  ///   2 * x1 + x3 / 2 <= 7
  ///   x2 + -3 * x3 >= 8
  ///   3 <= x2 - 2 * x1 <= 8
  /// bounds
  ///   x1 >= 3
  ///   2 <= x2 <= 5
  ///   0 <= x3 <= 8
  /// objective
  ///   min x1 + 2 * x2 - x3
  ///\endcode
  ///
  /// The first part gives the constraints to the lp. The constraints
  /// could be equality, lower bound, upper bound or range for an
  /// expression or equality, less or more for two expressions.
  ///
  /// The second part gives the bounds for the variables. The bounds
  /// could be the same as for the single expression so equality,
  /// lower bound, upper bound or range.
  ///
  /// The third part is the objective function, it should start with
  /// \c min or \c max and then a valid linear expression.
  ///
  /// The reader constructs automatically the variable for a name if
  /// there is not already a variable with the same name.
  ///
  /// The basic way to read an LP problem is made by the next code:
  ///\code
  /// Lp lp;
  ///
  /// LemonReader reader(filename_or_istream);
  /// LpReader lpreader(reader, lp);
  /// 
  /// reader.run();
  ///\endcode
  ///
  /// Of course instead of \c LemonReader you can give a \c
  /// GraphReader to the LpReader constructor. Also useful that you
  /// can mix lp problems and graphs in the same file.
  class LpReader : public LemonReader::SectionReader {
    typedef LemonReader::SectionReader Parent;
  public:


    /// \brief Constructor.
    ///
    /// Constructor for LpReader. It creates the LpReader and attachs
    /// it into the given LemonReader. The lp reader will add
    /// variables, constraints and objective function to the
    /// given lp solver.
    LpReader(LemonReader& _reader, LpSolverBase& _lp, 
             const std::string& _name = std::string())
      : Parent(_reader), lp(_lp), name(_name) {} 


    /// \brief Destructor.
    ///
    /// Destructor for LpReader.
    virtual ~LpReader() {}

  private:
    LpReader(const LpReader&);
    void operator=(const LpReader&);
  
  protected:

    /// \brief Gives back true when the SectionReader can process 
    /// the section with the given header line.
    ///
    /// It gives back true when the header line starts with \c \@lp,
    /// and the header line's name and the nodeset's name are the same.
    virtual bool header(const std::string& line) {
      std::istringstream ls(line);
      std::string command;
      std::string id;
      ls >> command >> id;
      return command == "@lp" && name == id;
    }

  private:

    enum Part {
      CONSTRAINTS, BOUNDS, OBJECTIVE
    };

    enum Relation {
      LE, EQ, GE
    };

    std::istream& readConstraint(std::istream& is) {
      LpSolverBase::Constr c;
      char x;
      LpSolverBase::Expr e1, e2;
      Relation op1;
      is >> std::ws;
      readExpression(is, e1);
      is >> std::ws;
      readRelation(is, op1);
      is >> std::ws;
      readExpression(is, e2);
      is >> std::ws;
      if (!is.get(x)) {
        if (op1 == LE) {
          if (e1.size() == 0) {
            c = e2 >= e1;
          } else {
            c = e1 <= e2;
          }
          c = e1 <= e2;
        } else if (op1 == GE) {
          if (e1.size() == 0) {
            c = e2 <= e1;
          } else {
            c = e1 >= e2;
          }
        } else {
          if (e1.size() == 0) {
            c = e2 == e1;
          } else {
            c = e1 == e2;
          }
        }
      } else {
        is.putback(x);
        LpSolverBase::Expr e3;
        Relation op2;
        readRelation(is, op2);
        is >> std::ws;
        readExpression(is, e3);
        if (!e1.empty() || !e3.empty()) {
          throw DataFormatError("Wrong range format");
        }
        if (op2 != op1 || op1 == EQ) {
          throw DataFormatError("Wrong range format");
        }
        if (op1 == LE) {
          c = e1.constComp() <= e2 <= e3.constComp();
        } else {
          c = e1.constComp() >= e2 >= e3.constComp();
        }
        is >> std::ws;
        if (is.get(x)) {
          throw DataFormatError("Wrong variable bounds format");
        }
      }
      lp.addRow(c);
      return is;
    }

    std::istream& readObjective(std::istream& is) {
      is >> std::ws;
      std::string sense;
      is >> sense;
      if (sense != "min" && sense != "max") {
        throw DataFormatError("Wrong objective function format");
      }
      LpSolverBase::Expr expr;
      is >> std::ws;
      readExpression(is, expr);
      lp.obj(expr);
      if (sense == "min") {
        lp.min();
      } else {
        lp.max();
      }
      return is;
    }

    std::istream& readBounds(std::istream& is) {
      LpSolverBase::Col c;
      char x;
      is >> std::ws;
      if (!is.get(x)) {
        throw DataFormatError("Wrong variable bounds format");
      }
      if (varFirstChar(x)) {
        is.putback(x);
        readCol(is, c);
        is >> std::ws;
        Relation op1;
        readRelation(is, op1);
        double v;
        readNum(is, v);
        if (op1 == EQ) {
          lp.colLowerBound(c, v);
          lp.colUpperBound(c, v);
        } else if (op1 == LE) {
          lp.colUpperBound(c, v);
        } else {
          lp.colLowerBound(c, v);
        }
        is >> std::ws;
        if (is.get(x)) {
          throw DataFormatError("Wrong variable bounds format");
        }
      } else {
        is.putback(x);
        double v;
        readNum(is, v);
        is >> std::ws;
        Relation op1;
        readRelation(is, op1);
        is >> std::ws;
        readCol(is, c);
        is >> std::ws;
        if (is.get(x)) {
          is.putback(x);
          is >> std::ws;
          Relation op2;
          readRelation(is, op2);
          double w;
          is >> std::ws;
          readNum(is, w);
          if (op2 != op1 || op1 == EQ) {
            throw DataFormatError("Wrong range format");
          }
          if (op1 == LE) {
            lp.colLowerBound(c, v);
            lp.colUpperBound(c, w);
          } else {
            lp.colLowerBound(c, w);
            lp.colUpperBound(c, v);
          }
          is >> std::ws;
          if (is.get(x)) {
            throw DataFormatError("Wrong variable bounds format");
          }
        } else {
          if (op1 == EQ) {
            lp.colLowerBound(c, v);
            lp.colUpperBound(c, v);
          } else if (op1 == LE) {
            lp.colLowerBound(c, v);
          } else {
            lp.colUpperBound(c, v);
          }
        }
      }
      return is;
    }

    std::istream& readExpression(std::istream& is, LpSolverBase::Expr& e) {
      LpSolverBase::Col c;
      double d;
      char x;
      readElement(is, c, d);
      if (c != INVALID) {
        e += d * c;
      } else {
        e += d;
      }
      is >> std::ws;
      while (is.get(x) && (x == '+' || x == '-')) {
        is >> std::ws;
        readElement(is, c, d);
        if (c != INVALID) {
          e += (x == '+' ? d : -d) * c;
        } else {
          e += (x == '+' ? d : -d);
        }
        is >> std::ws;
      }
      if (!is) {
        is.clear();
      } else {
        is.putback(x);
      }
      return is;
    }

    std::istream& readElement(std::istream& is, 
                              LpSolverBase::Col& c, double& d) { 
      d = 1.0;
      c = INVALID;
      char x, y;
      if (!is.get(x)) throw DataFormatError("Cannot find lp element");
      if (x == '+' || x == '-') {
        is >> std::ws;
        d *= x == '-' ? -1 : 1;
        while (is.get(x) && (x == '+' || x == '-')) {
          d *= x == '-' ? -1 : 1;
          is >> std::ws;
        }
        if (!is) throw DataFormatError("Cannot find lp element");
      }
      if (numFirstChar(x)) {
        is.putback(x);
        double e;
        readNum(is, e);
        d *= e;
      } else if (varFirstChar(x)) {
        is.putback(x);
        LpSolverBase::Col f;
        readCol(is, f);
        c = f;
      } else {
        throw DataFormatError("Invalid expression format");          
      }
      is >> std::ws;
      while (is.get(y) && (y == '*' || y == '/')) {
        is >> std::ws;
        if (!is.get(x)) throw DataFormatError("Cannot find lp element");
        if (x == '+' || x == '-') {
          is >> std::ws;
          d *= x == '-' ? -1 : 1;
          while (is.get(x) && (x == '+' || x == '-')) {
            d *= x == '-' ? -1 : 1;
            is >> std::ws;
          }
          if (!is) throw DataFormatError("Cannot find lp element");
        }
        if (numFirstChar(x)) {
          is.putback(x);
          double e;
          readNum(is, e);
          if (y == '*') {
            d *= e;
          } else {
            d /= e;
          }
        } else if (varFirstChar(x)) {
          is.putback(x);
          LpSolverBase::Col f;
          readCol(is, f);
          if (y == '*') {
            if (c == INVALID) {
              c = f;
            } else {
              throw DataFormatError("Quadratic element in expression");
            }
          } else {
            throw DataFormatError("Division by variable");
          }
        } else {
          throw DataFormatError("Invalid expression format");          
        }
        is >> std::ws;
      }
      if (!is) {
        is.clear();
      } else {
        is.putback(y);
      }
      return is;
    }

    std::istream& readCol(std::istream& is, LpSolverBase::Col& c) {
      char x;
      std::string var;
      while (is.get(x) && varChar(x)) {
        var += x;
      }
      if (!is) {
        is.clear();
      } else {
        is.putback(x);
      }
      c = lp.colByName(var);
      if (c == INVALID) {
        c = lp.addCol();
        lp.colName(c, var);
      }
      return is;
    }

    std::istream& readNum(std::istream& is, double& d) {
      is >> d;
      if (!is) throw DataFormatError("Wrong number format");
      return is;
    }

    std::istream& readRelation(std::istream& is, Relation& op) {
      char x, y;
      if (!is.get(x) || !(x == '<' || x == '=' || x == '>')) {
        throw DataFormatError("Wrong relation operator");
      }
      if (!is.get(y) || y != '=') {
        throw DataFormatError("Wrong relation operator");
      }
      switch (x) {
      case '<': op = LE; 
        break;
      case '=': op = EQ; 
        break;
      case '>': op = GE; 
        break;
      }
      return is;
    }

    static bool relationFirstChar(char c) {
      return c == '<' || c == '=' || c == '>';
    }

    static bool varFirstChar(char c) {
      return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
    }

    static bool numFirstChar(char c) {
      return (c >= '0' && c <= '9') || c == '.';
    }
    
    static bool varChar(char c) {
      return (c >= '0' && c <= '9') || 
        (c >= 'a' && c <= 'z') || 
        (c >= 'A' && c <= 'Z') ||
        c == '[' || c == ']';
    }
    
  protected:

    /// \brief Reader function of the section.
    ///
    /// It reads the content of the section.
    virtual void read(std::istream& is) {
      std::string line;
      Part part = CONSTRAINTS;
      while (getline(is, line)) {
        std::istringstream ls(line);
        std::string type;
        ls >> type;
        if (type == "constraints") {
          part = CONSTRAINTS;
          ls >> std::ws;
          char x;
          if (ls.get(x))
            throw DataFormatError("Wrong Lp format");
        } else if (type == "bounds") {
          part = BOUNDS;
          ls >> std::ws;
          char x;
          if (ls.get(x))
            throw DataFormatError("Wrong Lp format");
        } else if (type == "objective") {
          part = OBJECTIVE;
          ls >> std::ws;
          char x;
          if (ls.get(x))
            throw DataFormatError("Wrong Lp format");
        } else {
          ls.str(line);
          switch (part) {
          case CONSTRAINTS:
            readConstraint(ls);
            break;
          case BOUNDS:
            readBounds(ls);
            break;
          case OBJECTIVE:
            readObjective(ls);
            break;
          }
        }
      }
    }
      
    virtual void missing() {
      ErrorMessage msg;
      msg << "Lp section not found in file: @lp " << name;
      throw IoParameterError(msg.message());
    }

  private:

      
    LpSolverBase& lp;
    std::string name;
  };


  /// \ingroup lp_utils
  ///
  /// \brief Lp section writer for lemon IO
  ///
  /// This section reader provides an easy way to write an Lp problem
  /// to a file. The lemon lp section format contains three parts.
  ///
  ///\code
  /// @lp
  /// constraints
  ///   7 == x1 - 1 * x2
  ///   2 * x1 + x3 / 2 <= 7
  ///   x2 + -3 * x3 >= 8
  ///   3 <= x2 - 2 * x1 <= 8
  /// bounds
  ///   x1 >= 3
  ///   2 <= x2 <= 5
  ///   0 <= x3 <= 8
  /// objective
  ///   min x1 + 2 * x2 - x3
  ///\endcode
  ///
  /// The first part gives the constraints to the lp. The constraints
  /// could be equality, lower bound, upper bound or range for an
  /// expression or equality, less or more for two expressions.
  ///
  /// The second part gives the bounds for the variables. The bounds
  /// could be the same as for the single expression so equality,
  /// lower bound, upper bound or range.
  ///
  /// The third part is the objective function, it should start with
  /// \c min or \c max and then a valid linear expression.
  ///
  /// If an LP variable does not have name in the writer, then it will
  /// automatically created in the writer. This make a slight change
  /// in the \c const variable.
  ///
  /// The basic way to write an LP problem is made by the next code:
  ///\code
  /// Lp lp;
  ///
  /// LemonWriter writer(filename_or_ostream);
  /// LpWriter lpwriter(writer, lp);
  /// 
  /// writer.run();
  ///\endcode
  ///
  /// Of course instead of \c LemonWriter you can give a \c
  /// GraphWriter to the LpWriter constructor. Also useful that you
  /// can mix lp problems and graphs in the same file.
  class LpWriter : public LemonWriter::SectionWriter {
    typedef LemonWriter::SectionWriter Parent;
  public:


    /// \brief Constructor.
    ///
    /// Constructor for LpWriter. It creates the LpWriter and attach
    /// it into the given LemonWriter.
    LpWriter(LemonWriter& _writer, const LpSolverBase& _lp, 
             const std::string& _name = std::string())
      : Parent(_writer), lp(_lp), name(_name) {} 


    /// \brief Destructor.
    ///
    /// Destructor for LpWriter.
    virtual ~LpWriter() {}

  private:
    LpWriter(const LpWriter&);
    void operator=(const LpWriter&);
  
  protected:

    /// \brief Gives back true when the SectionWriter can process 
    /// the section with the given header line.
    ///
    /// It gives back the header line of the \c \@lp section.
    virtual std::string header() {
      std::ostringstream ls;
      ls << "@lp " << name;
      return ls.str();
    }

  private:

    void createCols() {
      int num = 0;

      for (LpSolverBase::ColIt it(lp); it != INVALID; ++it) {
        std::string name = lp.colName(it);
        if (name.empty()) {
          std::ostringstream ls;
          ls << "x" << num;
          ++num;
          while (lp.colByName(ls.str()) != INVALID) {
            ls.str(std::string());
            ls << "x" << num;
            ++num;
          }
          const_cast<LpSolverBase&>(lp).colName(it, ls.str());
        }
      }
    }
    
    void writeExpression(std::ostream& os, const LpSolverBase::Expr& e) {
      bool first = true;
      for (LpSolverBase::Expr::const_iterator jt = e.begin(); 
           jt != e.end(); ++jt) {
        if (jt->second < 0.0) {
          os << "- ";
          if (jt->second != -1.0) {
            os << -jt->second << " * ";
          }
        } else if (jt->second > 0.0) {
          if (!first) {
            os << "+ ";
          }
          if (jt->second != 1.0) {
            os << jt->second << " * ";
          }          
        }
        first = false;
        os << lp.colName(jt->first) << " ";
      }
      if (e.constComp() < 0.0) {
        os << "- " << -e.constComp() << " ";
      } else if (e.constComp() > 0.0) {
        if (!first) {
          os << "+ ";
        }
        os << e.constComp() << " ";
      }
      if (e.begin() == e.end() && e.constComp() == 0.0) {
        os << "0 ";
      }
    }

  protected:

    /// \brief Writer function of the section.
    ///
    /// It writes the content of the section.
    virtual void write(std::ostream& os) {
      createCols();

      os << "constraints" << std::endl;
      
      for (LpSolverBase::RowIt it(lp); it != INVALID; ++it) {
        double lower, upper;
        lp.getRowBounds(it, lower, upper);
        if (lower == -LpSolverBase::INF && upper == LpSolverBase::INF) 
          continue;
        os << "   ";
        
        if (lower != upper) {
          if (lower != -LpSolverBase::INF) {
            os << lower << " <= ";
          }
          
          writeExpression(os, lp.row(it));
          
          if (upper != LpSolverBase::INF) {
            os << "<= " << upper << " ";
          }
        } else {

          writeExpression(os, lp.row(it));
          os << "== " << upper << " ";
          
        }

        os << std::endl;
      }

      os << "bounds" << std::endl;
      
      for (LpSolverBase::ColIt it(lp); it != INVALID; ++it) {
        double lower = lp.colLowerBound(it);
        double upper = lp.colUpperBound(it);
        if (lower == -LpSolverBase::INF && upper == LpSolverBase::INF) 
          continue;
        os << "   ";

        if (lower != upper) {
          if (lower != -LpSolverBase::INF) {
            os << lower << " <= ";
          }
          
          os << lp.colName(it) << " ";
          
          if (upper != LpSolverBase::INF) {
            os << "<= " << upper << " ";
          }
        } else {
          os << lp.colName(it) << " == " << upper << " ";
        }
        os << std::endl;
      }

      os << "objective" << std::endl;
      os << "   ";
      if (lp.isMin()) {
        os << "min ";
      } else {
        os << "max ";
      }
      writeExpression(os, lp.obj());
      os << std::endl;
    }
      

  private:

    const LpSolverBase& lp;
    std::string name;
  };

}

#endif
