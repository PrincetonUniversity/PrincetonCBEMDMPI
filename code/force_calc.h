/*!
 \file force_calc.h
 \brief Header for force_calc function
*/

#ifndef FORCECALC_H_
#define FORCECALC_H_

#include "system.h"
#include "atom.h"
#include "interaction.h"
using namespace sim_system;

//! Calculates the forces between the particles in the system
int force_calc(System *sys);

#endif
