#  Copyright (C) 2014
#      Pierre de Buyl
#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

r"""
*************************************************
**espressopp.interaction.ZConfinement**
*************************************************

.. math::

    U = \frac{k}{2}(z - L/2)^2


.. function:: espressopp.interaction.ZConfinement(k)

        :param K:
        :type K:

.. function:: espressopp.interaction.SingleParticleZConfinement(system, potential)

        :param system:
        :param potential:
        :type system:
        :type potential:

.. function:: espressopp.interaction.SingleParticleZConfinement.setPotential(potential)

        :param potential:
        :type potential:
"""
from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.SingleParticlePotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_ZConfinement, interaction_SingleParticleZConfinement


class ZConfinementLocal(SingleParticlePotentialLocal, interaction_ZConfinement):

    def __init__(self, K=0.0):
        """Initialize the local ZConfinement object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_ZConfinement, K)

class SingleParticleZConfinementLocal(InteractionLocal, interaction_SingleParticleZConfinement):

    def __init__(self, system, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_SingleParticleZConfinement, system, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

    def getPotential(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self)

if pmi.isController:
    class ZConfinement(SingleParticlePotential):
        'The ZConfinement potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.ZConfinementLocal',
            pmiproperty = ['K']
            )

    class SingleParticleZConfinement(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.SingleParticleZConfinementLocal',
            pmicall = ['setPotential', 'getPotential']
            )
