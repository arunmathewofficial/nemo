if(PION_NEMO)
  add_library(
    nemo

      nemo.cpp
      mellema_cooling.cpp
      atomic_physics_database.cpp
      coll_ionise_recomb.cpp
      chianti_cooling.cpp
      stellar_atmosphere_models.cpp
      multi_ion_photo_xsections.cpp
      photoionisation.cpp
      charge_exchange.cpp
  )
endif(PION_NEMO)


