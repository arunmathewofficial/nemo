
#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H

#ifdef NEMO_OMP
#include <memory>
#include <mutex>
#include <omp.h>
#include <spdlog/sinks/base_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#endif

#include "nemo-datastructs.h"
#include "nemo/nemo.h"
using namespace NEMO;


#ifdef NEMO_OMP
// Passes on messages from OpenMP thread 0, and warnings/errors from any thread.
class thread0_sink : public spdlog::sinks::base_sink<std::mutex> {
public:
  thread0_sink() : out(std::make_shared<spdlog::sinks::stdout_color_sink_mt>())
  {
  }

protected:
  void sink_it_(const spdlog::details::log_msg &msg) override
  {
    if (omp_get_thread_num() == 0 || msg.level >= spdlog::level::warn)
      out->log(msg);
  }
  void flush_() override { out->flush(); }

private:
  std::shared_ptr<spdlog::sinks::stdout_color_sink_mt> out;
};
#endif


// Calculate optical depth to the inner edge of every shell.
void spherical_optical_depth(
    class nemo &,                  ///< pointer to NEMO class
    const nemo_star_properties &,  ///< struct containing stellar properties
    const std::vector<std::vector<double> > &,  ///< vector of state vectors
    std::vector<std::vector<struct nemo_cell_rt_data> > &  ///< vec of cell data
);


#endif  // COMMON_FUNCTIONS_H
