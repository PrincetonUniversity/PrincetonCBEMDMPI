/*!
 \file force_calc.h
 \brief Header for force_calc function
 \author {Frank Ricci, Jun Park}
*/

#ifndef FORCECALC_H_
#define FORCECALC_H_

#include "system.h"
#include "atom.h"
#include "interaction.h"

//! Calculates the forces between the particles in the system
int force_calc(System *sys);

//! Move atoms between processors (domains)
int send_atoms(System *sys);

#endif
