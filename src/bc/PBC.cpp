// this macro removes all LOG statements with level < WARN at source level
#define LOG4ESPP_LEVEL_WARN

#include "PBC.hpp"
#include <python.hpp>
#include <iostream>

using namespace espresso::bc;

/* ---------------------------------------------------------------------- */

LOG4ESPP_LOGGER(PBC::theLogger, "bc.PBC");

/* ---------------------------------------------------------------------- */

static inline real dround(real x) { return floor(x + 0.5); }

PMI_REGISTER_CLASS(PBC);

PBC::PBC() {}

PBC::PBC(real _length) {
  setLocal(_length);
}

PBC::~PBC() {}

void PBC::set(real _length)
#ifndef HAVE_MPI
{  setLocal(_length); }
#else
{
  //real v = _length;
  real v[3] = { 1.0, 2.0, 3.0 };
  pmiObject.invoke<&PBC::setWorker>();
  boost::mpi::communicator world;
  boost::mpi::broadcast(world, v, 3, pmi::getControllerMPIRank());
  setLocal(0.0);
}

void PBC::setWorker() {
  //real v;
  real v[3];
  boost::mpi::communicator world;
  boost::mpi::broadcast(world, v, 3, pmi::getControllerMPIRank());
  setLocal(0.0);
}

PMI_REGISTER_METHOD(PBC, setWorker);
#endif

void PBC::setLocal(real _length) {
  length = _length;
  lengthInverse = 1.0 / length;
}

/** Routine delivers the distance vector between two positions.
    \sa bc::BC::getDist */
Real3D PBC::getDist(const Real3D& pos1, const Real3D& pos2) const {

  real xij;
  real yij;
  real zij;

  xij = pos1.getX() - pos2.getX();
  yij = pos1.getY() - pos2.getY();
  zij = pos1.getZ() - pos2.getZ();

  xij -= dround(xij*lengthInverse)*length;
  yij -= dround(yij*lengthInverse)*length;
  zij -= dround(zij*lengthInverse)*length;

  return Real3D(xij, yij, zij);
}

#ifdef HAVE_PYTHON
//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void
PBC::registerPython() {
  using namespace boost::python;

  class_<PBC>("bc_PBC", init<>())
    .def("set", &PBC::set)
    .def("getDist", &PBC::getDist);
}
#endif
