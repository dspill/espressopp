/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include "python.hpp"
#include "VolumeScaler.hpp"
#include "System.hpp"
#include "esutil/Error.hpp"
#include "bc/BC.hpp"

namespace espressopp {
  
  using namespace analysis;
  using namespace esutil;
  using namespace std;
  
  namespace integrator {

    LOG4ESPP_LOGGER(VolumeScaler::theLogger, "VolumeScaler");

    VolumeScaler::VolumeScaler(shared_ptr<System> system): Extension(system){
      perc  = 0.;
      L0 = 1.0;
      
      type = Extension::Barostat;

      LOG4ESPP_INFO(theLogger, "VolumeScaler constructed");
    }

    VolumeScaler::~VolumeScaler(){
      LOG4ESPP_INFO(theLogger, "~VolumeScaler");
      disconnect();
    }
    
    void VolumeScaler::disconnect(){
      _runInit.disconnect();
      _aftIntV.disconnect();
    }

    void VolumeScaler::connect(){
      // connection to initialisation
      _runInit = integrator->runInit.connect( boost::bind(&VolumeScaler::initialize, this));

      // connection to the signal at the end of the run
      _aftIntV = integrator->aftIntV.connect( boost::bind(&VolumeScaler::barostat, this));
    }

    // set and get time constant for VolumeScaler barostat
    void VolumeScaler::setPerc(real _perc) {
      if(perc < 0) throw std::runtime_error("perc must be positive");
      perc = _perc;
    }
    real VolumeScaler::getPerc() {
      return perc;
    }
    // set and get external pressure
    void VolumeScaler::setL0(real _L0) {
      L0 = _L0; // TODO make multidim
    }
    real VolumeScaler::getL0() {
      return L0;
    }

    void VolumeScaler::barostat(){
      LOG4ESPP_DEBUG(theLogger, "scaling volume");
      System& system = getSystemRef();
      Real3D L = system.bc->getBoxL();

      Real3D signs;
      signs[0] = ((L0 - L[0]) > 0) - ((L0 - L[0]) < 0);
      signs[1] = ((L0 - L[1]) > 0) - ((L0 - L[1]) < 0);
      signs[2] = ((L0 - L[2]) > 0) - ((L0 - L[2]) < 0);

      signs[0] *= std::min(fabs(L0 - L[0]), 1.);
      signs[1] *= std::min(fabs(L0 - L[1]), 1.);
      signs[2] *= std::min(fabs(L0 - L[2]), 1.);

      system.scaleVolume(Real3D(1.) + signs * perc, true);
    }

     // calculate the prefactors
    void VolumeScaler::initialize(){
      LOG4ESPP_INFO(theLogger, "init, perc = " << perc << 
                               ", external pressure = " << L0);
      real dt = integrator->getTimeStep();
      //pref = dt / perc;
    }


    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void VolumeScaler::registerPython() {

      using namespace espressopp::python;

      class_<VolumeScaler, shared_ptr<VolumeScaler>, bases<Extension> >

        ("integrator_VolumeScaler", init< shared_ptr<System> >())

        .add_property("perc",
              &VolumeScaler::getPerc,
              &VolumeScaler::setPerc)
        .add_property("L0",
              &VolumeScaler::getL0, 
              &VolumeScaler::setL0)
      
        .def("connect", &VolumeScaler::connect)
        .def("disconnect", &VolumeScaler::disconnect)
      ;
    }

  }
}

