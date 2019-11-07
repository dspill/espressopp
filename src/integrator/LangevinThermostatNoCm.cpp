/*
   Copyright (C) 2012,2013,2014,2015,2016
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
#include "LangevinThermostatNoCm.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"
#include "storage/DomainDecomposition.hpp"

namespace espressopp {

    namespace integrator {

        using namespace espressopp::iterator;

        LangevinThermostatNoCm::LangevinThermostatNoCm(shared_ptr<System> system,
                int _number_of_chains)
            :Extension(system) {

                type = Extension::Thermostat;

                gamma  = 0.0;
                temperature = 0.0;
                number_of_chains = _number_of_chains;

                adress = false;

                // Determine total number of particles in the system
                int Npart = system->storage->getNRealParticles();
                totNPart = 0;
                boost::mpi::all_reduce(*getSystem()->comm, Npart, totNPart, std::plus<int>());
                if(totNPart == 0)
                {
                    throw std::runtime_error("You need to fill the system with "
                            "particles before you set up LangevinThermostatNoCm.");
                }
                if(totNPart % number_of_chains != 0)
                {
                    throw std::runtime_error("Number of chains is not"
                            " commensurable with number of particles.");
                }
                degree_of_polymerization = totNPart / number_of_chains;

                if (!system->rng) {
                    throw std::runtime_error("system has no RNG");
                }
                rng = system->rng;

                excess_forces = (Real3D*) calloc(number_of_chains, sizeof(Real3D));
                total_excess_forces = (Real3D*) calloc(number_of_chains, sizeof(Real3D));

                LOG4ESPP_INFO(theLogger, "Langevin constructed");
            }

        void LangevinThermostatNoCm::setGamma(real _gamma)
        {
            gamma = _gamma;
        }

        real LangevinThermostatNoCm::getGamma()
        {
            return gamma;
        }

        void LangevinThermostatNoCm::setAdress(bool _adress)
        {
            adress = _adress;
        }

        bool LangevinThermostatNoCm::getAdress()
        {
            return adress;
        }

        void LangevinThermostatNoCm::setTemperature(real _temperature)
        {
            temperature = _temperature;
        }

        real LangevinThermostatNoCm::getTemperature()
        {
            return temperature;
        }

        LangevinThermostatNoCm::~LangevinThermostatNoCm()
        {
            disconnect();
        }

        void LangevinThermostatNoCm::disconnect() 
        {

            _initialize.disconnect();
            _heatUp.disconnect();
            _coolDown.disconnect();
            _thermalize.disconnect();
            _thermalizeAdr.disconnect();

        }

        void LangevinThermostatNoCm::connect() 
        {

            // connect to initialization inside run()
            _initialize = integrator->runInit.connect(
                    boost::bind(&LangevinThermostatNoCm::initialize, this));

            _heatUp = integrator->recalc1.connect(
                    boost::bind(&LangevinThermostatNoCm::heatUp, this));

            _coolDown = integrator->recalc2.connect(
                    boost::bind(&LangevinThermostatNoCm::coolDown, this));

            if (adress) {
                _thermalizeAdr = integrator->aftCalcF.connect(
                        boost::bind(&LangevinThermostatNoCm::thermalizeAdr, this));
            }
            else {
                _thermalize = integrator->aftCalcF.connect(
                        boost::bind(&LangevinThermostatNoCm::thermalize, this));
            }
        }

        void LangevinThermostatNoCm::thermalize()
        {
            LOG4ESPP_DEBUG(theLogger, "thermalize");

            // reinitialize arrays
            std::fill(excess_forces, excess_forces + number_of_chains, 0);
            std::fill(total_excess_forces, total_excess_forces + number_of_chains, 0);

            System& system = getSystemRef();
            int i_part, i_chain;

            CellList cells = system.storage->getRealCells();

            for(CellListIterator cit(cells); !cit.isDone(); ++cit) 
            {
                i_part = cit->id();
                i_chain = i_part / degree_of_polymerization;

                excess_forces[i_chain] += frictionThermo(*cit);
            }

            boost::mpi::all_reduce(*getSystem()->comm, excess_forces,
                    number_of_chains, total_excess_forces, std::plus<Real3D>());

            //int my_rank = 0;
            //MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
            //printf("nchains %d\n", number_of_chains);
            //printf("dop = %d\n", degree_of_polymerization);

            //for(int i = 0; i < number_of_chains; ++i)
            //{
                //std::cout << total_excess_forces[i] << ' ';
            //}
            //std::cout << std::endl;

            for(CellListIterator cit(cells); !cit.isDone(); ++cit) 
            {
                i_part = cit->id();
                i_chain = (i_part - 1) / degree_of_polymerization;
                //printf("%3d %2d %4.2f\n", i_part, i_chain, excess_forces[i_chain].abs());

                cit->force() -= total_excess_forces[i_chain];
            }
        }

        // for AdResS
        void LangevinThermostatNoCm::thermalizeAdr()
        {
            LOG4ESPP_DEBUG(theLogger, "thermalize");

            System& system = getSystemRef();

            // thermalize AT particles
            ParticleList& adrATparticles = system.storage->getAdrATParticles();
            for (std::vector<Particle>::iterator it = adrATparticles.begin();
                    it != adrATparticles.end(); it++) {

                frictionThermo(*it);

            }
        }

        Real3D LangevinThermostatNoCm::frictionThermo(Particle& p)
        {
            Real3D excess_force;
            real massf = sqrt(p.mass());

            // get a random value for each vector component
            Real3D ranval((*rng)() - 0.5, (*rng)() - 0.5, (*rng)() - 0.5);

            excess_force = pref1 * p.velocity() * p.mass() +
                pref2 * ranval * massf;

            p.force() += excess_force;

            LOG4ESPP_TRACE(theLogger, "new force of p = " << p.force());

            return excess_force;
        }

        void LangevinThermostatNoCm::initialize()
        { // calculate the prefactors

            real timestep = integrator->getTimeStep();

            LOG4ESPP_INFO(theLogger, "init, timestep = " << timestep <<
                    ", gamma = " << gamma <<
                    ", temperature = " << temperature);

            pref1 = -gamma;
            pref2 = sqrt(24.0 * temperature * gamma / timestep);

        }

        /** very nasty: if we recalculate force when leaving/reentering the integrator,
          a(t) and a((t-dt)+dt) are NOT equal in the vv algorithm. The random
          numbers are drawn twice, resulting in a different variance of the random force.
          This is corrected by additional heat when restarting the integrator here.
          Currently only works for the Langevin thermostat, although probably also others
          are affected.
          */

        void LangevinThermostatNoCm::heatUp()
        {
            LOG4ESPP_INFO(theLogger, "heatUp");

            pref2buffer = pref2;
            pref2       *= sqrt(3.0);
        }

        /** Opposite to heatUp */

        void LangevinThermostatNoCm::coolDown()
        {
            LOG4ESPP_INFO(theLogger, "coolDown");

            pref2 = pref2buffer;
        }

        /****************************************************
         ** REGISTRATION WITH PYTHON
         ****************************************************/

        void LangevinThermostatNoCm::registerPython() {


            using namespace espressopp::python;


            class_<LangevinThermostatNoCm, shared_ptr<LangevinThermostatNoCm>, bases<Extension> >
                ("integrator_LangevinThermostatNoCm", init< shared_ptr<System>, int >())
                .def("connect", &LangevinThermostatNoCm::connect)
                .def("disconnect", &LangevinThermostatNoCm::disconnect)
                .add_property("adress", &LangevinThermostatNoCm::getAdress, &LangevinThermostatNoCm::setAdress)
                .add_property("gamma", &LangevinThermostatNoCm::getGamma, &LangevinThermostatNoCm::setGamma)
                .add_property("temperature", &LangevinThermostatNoCm::getTemperature, &LangevinThermostatNoCm::setTemperature)
                ;


        }

    }
}

