/*
  Copyright (C) 2014
      Pierre de Buyl
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
#include "ZConfinement.hpp"
#include "Real3D.hpp"

namespace espressopp {
  namespace interaction {


    typedef class SingleParticleInteractionTemplate <ZConfinement>
    SingleParticleZConfinement;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    ZConfinement::registerPython() {
      using namespace espressopp::python;

      class_< ZConfinement, bases< SingleParticlePotential > >
        ("interaction_ZConfinement", init< real >())
        .add_property("K", &ZConfinement::getK, &ZConfinement::setK)
        ;

      class_< SingleParticleZConfinement, bases< Interaction > >
        ("interaction_SingleParticleZConfinement", init< shared_ptr<System>, shared_ptr<ZConfinement> >())
        .def("setPotential", &SingleParticleZConfinement::setPotential)
        .def("getPotential", &SingleParticleZConfinement::getPotential)
      ;
    }
  }
}
