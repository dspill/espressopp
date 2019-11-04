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

// ESPP_CLASS
#ifndef _INTERACTION_ZCONFINEMENT_HPP
#define _INTERACTION_ZCONFINEMENT_HPP

#include "SingleParticlePotential.hpp"
#include "SingleParticleInteractionTemplate.hpp"
#include "Real3D.hpp"
#include "bc/BC.hpp"
#include <cmath>
#include <cstdio>

namespace espressopp {
    namespace interaction {

        /** This class provides methods to compute forces and energies for an
          harmonic well potential.
          */
        class ZConfinement : public SingleParticlePotentialTemplate<ZConfinement> {
            private:
                real K;

            public:
                static void registerPython();

                ZConfinement() : K(0.0) {
                }

                ZConfinement(real _K) : K(_K) {
                }

                //~ZConfinement() {};
                int bondType() { return Single; }
                real getMaxCutoff() { return 0.; }

                // Setter and getter
                void setK(real _K) {
                    K = _K;
                }
                real getK() const { return K; }


                real _computeEnergyRaw(const Particle& p, const bc::BC& bc) const {
                    real L = bc.getBoxL()[2];
                    real z = p.position()[2];

                    return K/2*(1 - cos(2*M_PI*z/L));
                } 

                bool _computeForceRaw(Real3D& force,
                        const Particle& p,
                        const bc::BC& bc) const {
                    real L = bc.getBoxL()[2];
                    real z = p.position()[2];
                    real f = -K*M_PI/L*sin(2*M_PI*z/L);

                    force = Real3D(0, 0, f);

                    return true;
                }
        };
    }
}

#endif
