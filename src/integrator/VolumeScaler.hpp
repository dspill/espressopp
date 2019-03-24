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

#ifndef _INTEGRATOR_VOLUMESCALER_HPP
#define	_INTEGRATOR_VOLUMESCALER_HPP

#include "types.hpp"
#include "logging.hpp"

#include "analysis/Pressure.hpp"
#include "Extension.hpp"

#include "boost/signals2.hpp"
#include "Int3D.hpp"

namespace espressopp {
  
  using namespace analysis;

  namespace integrator {

    class VolumeScaler: public Extension {

      public:
        VolumeScaler(shared_ptr< System > system);
        
        void setFixed(Int3D);
        Int3D getFixed();
        void setPerc(real);
        real getPerc();
        void setL0(real);
        real getL0();

        ~VolumeScaler();

        void connect();
        void disconnect();
        
        /* Register in Python. */
        static void registerPython();

      private:
        boost::signals2::connection _runInit, _aftIntV;
        
        real perc;  // time constant
        real L0;    // external pressure

        void initialize();

        /* rescale the system size and coord. of particles */
        void barostat();

        /* Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif

