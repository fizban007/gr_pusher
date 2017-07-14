#ifndef _METRICS_H_
#define _METRICS_H_

enum class MetricType : char {
  Cartesian,
    spherical,
    log_spherical,
    Schwarzschild,
    // log_Schwarzschild,
    Kerr_Schild,
    // log_Kerr_Schild,
    Boyer_Lindquist,
    log_Boyer_Lindquist
    // rot_log_Schwarzschild,
    // rot_Minkowski
      };

#include "metrics/metric_cartesian.h"
#include "metrics/metric_spherical.h"
#include "metrics/metric_log_spherical.h"
#include "metrics/metric_schwarzschild.h"
#include "metrics/metric_kerr_schild.h"
#include "metrics/metric_boyer_lindquist.h"
#include "metrics/metric_log_boyer_lindquist.h"

namespace Aperture {

/// This function parses a metric description string, initializes the global
/// metric according to the description string and returns the initialized
/// metric type. This way the simulation environment simply calls parse_metric
/// and it will initialize the correct metric to the desired configuration, and
/// know which metric type it is.
MetricType parse_metric(const std::string& str);

/// This function selects the metric to use for a given templated function
/// object. The way to use this function is as follows. Assume you have a
/// templated functor called `func` which accepts a template parameter `Metric`,
/// and you want to use a dynamically given `MetricType` to determine which
/// version of `func` to call. Assuming `func` is defined in the following way
///
/// \code{.cpp}
/// struct Func {
///   template <typename Metric>
///   void operator() (const Metric& metric, int i) {
///     // do stuff here
///   }
/// };
/// \endcode
///
/// Now you can pass an instance of Func into this function, for example:
///
/// \code{.cpp}
/// SelectMetric(WeightType::CLOUD_IN_CELL, Func(), i);
/// \endcode
///
/// where `i` is the argument that will be passed to `func`. This way one
/// doesn't need to write out a if-else statement each time when
/// selecting the weight type dynamically.
template <typename Functor, typename... Args>
void select_metric(MetricType metric_type, Functor func, Args && ...args) {
  if (metric_type == MetricType::Cartesian) {
    func(metric::g_metric_cartesian, std::forward<Args>(args)...);
  // } else if (metric_type == MetricType::spherical) {
  //   func(metric::g_metric_spherical, std::forward<Args>(args)...);
  } else if (metric_type == MetricType::log_spherical) {
    func(metric::g_metric_log_spherical, std::forward<Args>(args)...);
 // } else if (metric_type == MetricType::Schwarzschild) {
  //   func(metric::g_metric_schwarzschild, std::forward<Args>(args)...);
  } else if (metric_type == MetricType::Kerr_Schild) {
    func(metric::g_metric_kerr_schild, std::forward<Args>(args)...);
  } else if (metric_type == MetricType::Boyer_Lindquist) {
    func(metric::g_metric_boyer_lindquist, std::forward<Args>(args)...);
  } else if (metric_type == MetricType::log_Boyer_Lindquist) {
    func(metric::g_metric_log_boyer_lindquist, std::forward<Args>(args)...);
  }
}


}

#endif  // _METRICS_H_
