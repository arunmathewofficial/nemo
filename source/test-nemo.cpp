/// test-nemo.cpp
///
/// For testing nemo module
///



using namespace std;
#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/spdlog.h>

#include "git_version.h"
#include "nemo/nemo.h"
#include "datastructs.h"

int main(int argc, char **argv)
{
#ifdef NDEBUG
  spdlog::set_level(spdlog::level::info);
  spdlog::flush_on(spdlog::level::err);
#else
  spdlog::set_level(spdlog::level::trace);
  spdlog::flush_on(spdlog::level::trace);
#endif

  spdlog::info("NEMO tests: git {} running", git_version);


  std::string tr[5] = {"X_H", "X_He", "H", "He", "He1+"};
  struct which_physics wp;
  wp.dynamics           =1;
  wp.raytracing         =1;
  wp.cooling            =1;
  wp.chemistry          =1;
  wp.coll_ionisation    =1;
  wp.rad_recombination  =1;
  wp.phot_ionisation    =1;
  wp.charge_exchange    =1;
  wp.sat_thermal_cond   =1;
  wp.tc_strength        =1;
  wp.update_erg         =1;
  wp.MP_timestep_limit  =1;
  wp.MinTemperature     =10.0;
  wp.MaxTemperature     =1.0e9;
  wp.H_MassFrac         =0.0;
  wp.Helium_MassFrac    =0.0; 
  wp.Metal_MassFrac     =0.0;
  wp.compton_cool       =0;
  wp.wind_acceleration  =0;
  wp.shocked_tracer     =0;
  
  struct rad_sources;
  rs.Nsources=0;

  class nemo MP(1,1,15,5, tr, &wp, &rs, 5.0/3.0); 
  
  return 0;
}


