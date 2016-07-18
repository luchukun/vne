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

#ifndef LEMON_SIMANN_H
#define LEMON_SIMANN_H

/// \ingroup experimental
/// \file
/// \brief Simulated annealing framework.
///
/// \todo A test and some demo should be added
/// \todo Doc should be improved
/// \author Akos Ladanyi

#include <cstdlib>
#include <lemon/math.h>
#include <limits>
#include <lemon/time_measure.h>
#include <lemon/random.h>

namespace lemon {

  class SimAnnBase;

  /// \brief A base class for controllers.
  class ControllerBase {
  public:
    friend class SimAnnBase;
    /// \brief Pointer to the simulated annealing base class.
    SimAnnBase *simann;
    /// \brief Initializes the controller.
    virtual void init() {}
    /// \brief This is called by the simulated annealing class when a
    /// neighbouring state gets accepted.
    virtual void acceptEvent() {}
    /// \brief This is called by the simulated annealing class when the
    /// accepted neighbouring state's cost is less than the best found one's.
    virtual void improveEvent() {}
    /// \brief This is called by the simulated annealing class when a
    /// neighbouring state gets rejected.
    virtual void rejectEvent() {}
    /// \brief Decides whether to continue the annealing process or not.
    virtual bool next() = 0;
    /// \brief Decides whether to accept the current solution or not.
    virtual bool accept() = 0;
    /// \brief Destructor.
    virtual ~ControllerBase() {}
  };

  /// \brief Skeleton of an entity class.
  class EntityBase {
  public:
    /// \brief Makes a minor change to the entity.
    /// \return the new cost
    virtual double mutate() = 0;
    /// \brief Restores the entity to its previous state i.e. reverts the
    /// effects of the last mutate().
    virtual void revert() = 0;
    /// \brief Makes a copy of the entity.
    virtual EntityBase* clone() = 0;
    /// \brief Makes a major change to the entity.
    virtual void randomize() = 0;
    /// \brief Destructor.
    virtual ~EntityBase() {}
  };

  /// \brief Simulated annealing abstract base class.
  ///
  /// Can be used to derive a custom simulated annealing class if \ref SimAnn
  /// doesn't fit your needs.
  class SimAnnBase {
  private:
    /// \brief Pointer to the controller.
    ControllerBase *controller;
    /// \brief Cost of the current solution.
    double curr_cost;
    /// \brief Cost of the best solution.
    double best_cost;
    /// \brief Cost of the previous solution.
    double prev_cost;
    /// \brief Cost of the solution preceding the previous one.
    double prev_prev_cost;
    /// \brief Number of iterations.
    long iter;
    /// \brief Number of iterations which did not improve the solution since
    /// the last improvement.
    long last_impr;
  protected:
    /// \brief Step to a neighbouring state.
    virtual double mutate() = 0;
    /// \brief Reverts the last mutate().
    virtual void revert() = 0;
    /// \brief Saves the current solution as the best one.
    virtual void saveAsBest() = 0;
    /// \brief Does initializations before each run.
    virtual void init() {
      controller->init();
      curr_cost = prev_cost = prev_prev_cost = best_cost =
        std::numeric_limits<double>::infinity();
      iter = last_impr = 0;
    }
  public:
    /// \brief Sets the controller class to use.
    void setController(ControllerBase &_controller) {
      controller = &_controller;
      controller->simann = this;
    }
    /// \brief Returns the cost of the current solution.
    double getCurrCost() const { return curr_cost; }
    /// \brief Returns the cost of the previous solution.
    double getPrevCost() const { return prev_cost; }
    /// \brief Returns the cost of the best solution.
    double getBestCost() const { return best_cost; }
    /// \brief Returns the number of iterations done.
    long getIter() const { return iter; }
    /// \brief Returns the ordinal number of the last iteration when the
    /// solution was improved.
    long getLastImpr() const { return last_impr; }
    /// \brief Performs one iteration.
    bool step() {
      iter++;
      prev_prev_cost = prev_cost;
      prev_cost = curr_cost;
      curr_cost = mutate();
      if (controller->accept()) {
        controller->acceptEvent();
        last_impr = iter;
        if (curr_cost < best_cost) {
          best_cost = curr_cost;
          saveAsBest();
          controller->improveEvent();
        }
      }
      else {
        revert();
        curr_cost = prev_cost;
        prev_cost = prev_prev_cost;
        controller->rejectEvent();
      }
      return controller->next();
    }
    /// \brief Performs a given number of iterations.
    /// \param n the number of iterations
    bool step(int n) {
      for(; n > 0 && step(); --n) ;
      return !n;
    }
    /// \brief Starts the annealing process.
    void run() {
      init();
      do { } while (step());
    }
    /// \brief Destructor.
    virtual ~SimAnnBase() {}
  };

  /// \ingroup metah
  ///
  /// \brief Simulated annealing class.
  class SimAnn : public SimAnnBase {
  private:
    /// \brief Pointer to the current entity.
    EntityBase *curr_ent;
    /// \brief Pointer to the best entity.
    EntityBase *best_ent;
    /// \brief Does initializations before each run.
    void init() {
      SimAnnBase::init();
      if (best_ent) delete best_ent;
      best_ent = NULL;
      curr_ent->randomize();
    }
  public:
    /// \brief Constructor.
    SimAnn() : curr_ent(NULL), best_ent(NULL) {}
    /// \brief Destructor.
    virtual ~SimAnn() {
      if (best_ent) delete best_ent;
    }
    /// \brief Step to a neighbouring state.
    double mutate() {
      return curr_ent->mutate();
    }
    /// \brief Reverts the last mutate().
    void revert() {
      curr_ent->revert();
    }
    /// \brief Saves the current solution as the best one.
    void saveAsBest() { 
      if (best_ent) delete best_ent;
      best_ent = curr_ent->clone();
    }
    /// \brief Sets the current entity.
    void setEntity(EntityBase &_ent) {
      curr_ent = &_ent;
    }
    /// \brief Returns a copy of the best found entity.
    EntityBase* getBestEntity() { return best_ent->clone(); }
  };

  /// \brief A simple controller for the simulated annealing class.
  ///
  /// This controller starts from a given initial temperature and evenly
  /// decreases it.
  class SimpleController : public ControllerBase {
  private:
    /// \brief Maximum number of iterations.
    long max_iter;
    /// \brief Maximum number of iterations which do not improve the
    /// solution.
    long max_no_impr;
    /// \brief Temperature.
    double temp;
    /// \brief Annealing factor.
    double ann_fact;
    /// \brief Constructor.
    /// \param _max_iter maximum number of iterations
    /// \param _max_no_impr maximum number of consecutive iterations which do
    ///        not yield a better solution
    /// \param _temp initial temperature
    /// \param _ann_fact annealing factor
  public:
    SimpleController(long _max_iter = 500000, long _max_no_impr = 20000,
    double _temp = 1000.0, double _ann_fact = 0.9999) : max_iter(_max_iter),
      max_no_impr(_max_no_impr), temp(_temp), ann_fact(_ann_fact)
    {
    }
    /// \brief This is called when a neighbouring state gets accepted.
    void acceptEvent() {}
    /// \brief This is called when the accepted neighbouring state's cost is
    /// less than the best found one's.
    void improveEvent() {}
    /// \brief This is called when a neighbouring state gets rejected.
    void rejectEvent() {}
    /// \brief Decides whether to continue the annealing process or not. Also
    /// decreases the temperature.
    bool next() {
      temp *= ann_fact;
      bool quit = (simann->getIter() > max_iter) ||
        (simann->getIter() - simann->getLastImpr() > max_no_impr);
      return !quit;
    }
    /// \brief Decides whether to accept the current solution or not.
    bool accept() {
      double cost_diff = simann->getCurrCost() - simann->getPrevCost();
      return (rnd() <= exp(-(cost_diff / temp)));
    }
    /// \brief Destructor.
    virtual ~SimpleController() {}
  };

  /// \brief A controller with preset running time for the simulated annealing
  /// class.
  ///
  /// With this controller you can set the running time of the annealing
  /// process in advance. It works the following way: the controller measures
  /// a kind of divergence. The divergence is the difference of the average
  /// cost of the recently found solutions the cost of the best found one. In
  /// case this divergence is greater than a given threshold, then we decrease
  /// the annealing factor, that is we cool the system faster. In case the
  /// divergence is lower than the threshold, then we increase the temperature.
  /// The threshold is a function of the elapsed time which reaches zero at the
  /// desired end time.
  class AdvancedController : public ControllerBase {
  private:
    /// \brief Timer class to measure the elapsed time.
    Timer timer;
    /// \brief Calculates the threshold value.
    /// \param time the elapsed time in seconds
    virtual double threshold(double time) {
      return (-1.0) * start_threshold / end_time * time + start_threshold;
    }
    /// \brief Parameter used to calculate the running average.
    double alpha;
    /// \brief Parameter used to decrease the annealing factor.
    double beta;
    /// \brief Parameter used to increase the temperature.
    double gamma;
    /// \brief The time at the end of the algorithm.
    double end_time;
    /// \brief The time at the start of the algorithm.
    double start_time;
    /// \brief Starting threshold.
    double start_threshold;
    /// \brief Average cost of recent solutions.
    double avg_cost;
    /// \brief Temperature.
    double temp;
    /// \brief Annealing factor.
    double ann_fact;
    /// \brief Initial annealing factor.
    double init_ann_fact;
    /// \brief True when the annealing process has been started.
    bool start;
  public:
    /// \brief Constructor.
    /// \param _end_time running time in seconds
    /// \param _alpha parameter used to calculate the running average
    /// \param _beta parameter used to decrease the annealing factor
    /// \param _gamma parameter used to increase the temperature
    /// \param _ann_fact initial annealing factor
    AdvancedController(double _end_time, double _alpha = 0.2,
    double _beta = 0.9, double _gamma = 1.6, double _ann_fact = 0.9999) :
    alpha(_alpha), beta(_beta), gamma(_gamma), end_time(_end_time),
    ann_fact(_ann_fact), init_ann_fact(_ann_fact), start(false)
    {
    }
    /// \brief Does initializations before each run.
    void init() {
      avg_cost = simann->getCurrCost();
    }
    /// \brief This is called when a neighbouring state gets accepted.
    void acceptEvent() {
      avg_cost = alpha * simann->getCurrCost() + (1.0 - alpha) * avg_cost;
      if (!start) {
        static int cnt = 0;
        cnt++;
        if (cnt >= 100) {
          // calculate starting threshold and starting temperature
          start_threshold = 5.0 * fabs(simann->getBestCost() - avg_cost);
          temp = 10000.0;
          start = true;
          timer.restart();
        }
      }
    }
    /// \brief Decides whether to continue the annealing process or not.
    bool next() {
      if (!start) {
        return true;
      }
      else {
        double elapsed_time = timer.realTime();
        if (fabs(avg_cost - simann->getBestCost()) > threshold(elapsed_time)) {
          // decrease the annealing factor
          ann_fact *= beta;
        }
        else {
          // increase the temperature
          temp *= gamma;
          // reset the annealing factor
          ann_fact = init_ann_fact;
        }
        temp *= ann_fact;
        return elapsed_time < end_time;
      }
    }
    /// \brief Decides whether to accept the current solution or not.
    bool accept() {
      if (!start) {
        return true;
      }
      else {
        double cost_diff = simann->getCurrCost() - simann->getPrevCost();
        return (rnd() <= exp(-(cost_diff / temp)));
      }
    }
    /// \brief Destructor.
    virtual ~AdvancedController() {}
  };

}

#endif
