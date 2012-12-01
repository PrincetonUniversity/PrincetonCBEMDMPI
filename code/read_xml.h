/**
 Read HOOMD XML file format into System object
 \author Nathan A. Mahynski
 **/

#ifndef READ_XML_H_
#define READ_XML_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <map>
#include "system.h"
#include "misc.h"
#include "atom.h"
#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>

using namespace std;
using namespace system;
using namespace misc;
using namespace atom;

/*!
 Parse an XML file to obtain atom information. Returns 0 if successful, -1 if failure. This only reads on the main node when MPI is used, MPI_Send's 
 according to domain decomposition.
 \param [in] filename Name of file to open and read.
 \param [in] nprocs Number of processors total
 \param [in,out] \*sys System object for the main node to stores its information at.
 */
int read_xml (const string filename, const int nprocs, System *sys) {
	FILE *input = mfopen(filename, "r");
	if (input == NULL) {
		sprintf(err_msg, "Could not initialize from %s", filename);
		flag_error (err_msg, __FILE__, __LINE__);
		return -1;
	}
	
	// read in data
	sprintf(err_msg, "Reading initial configuration from %s", filename);
	flag_notify (err_msg, __FILE__, __LINE__);
	
	// header
	check = 0;
	int natoms = -1;
	while (fgets(buff, buffsize, fp1) != NULL) {
		if (strstr(buff, "natoms") != NULL) {
			check = 1;
			vector <string> fields;
			split(fields, buff, is_any_of("=,\", "), token_compress_on);
			for (int i = 0; i < fields.size(); ++i) {
				if (fields[i] == "natoms") {
					if (i+1 < fields.size()) {
						natoms = atoi(fields[i+1].c_str());
					}
					break;
				}
			}
			break;
		}
	}
	if (check != 1) {
		sprintf(err_msg, "Could not locate natoms in %s", filename);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(fp1);
		return -1;
	}
	if (natoms < 0) {
		sprintf(err_msg, "Number of atoms < 0 in %s", filename);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(fp1);
		return -1;
	}
	sprintf(err_msg, "Found %d atoms in %s", natoms, filename);
	flag_notify (err_msg, __FILE__, __LINE__);
	vector<Atom> new_atoms(natoms);
	
	rewind(fp1);
	// read box
	check = 0;
	vector <double> box(3);
	while (fgets(buff, buffsize, fp1) != NULL) {
		if (strstr(buff, "box") != NULL) {
			check = 1;
			vector <string> fields;
			split(fields, buff, is_any_of("=,\", "), token_compress_on);
			for (int i = 0; i < fields.size(); ++i) {
				if (fields[i] == "lx") {
					if (i+1 < fields.size()) {
						box[0] = atof(fields[i+1].c_str());
					}
				}
				if (fields[i] == "ly") {
					if (i+1 < fields.size()) {
						box[1] = atof(fields[i+1].c_str());
					}
				}
				if (fields[i] == "lz") {
					if (i+1 < fields.size()) {
						box[2] = atof(fields[i+1].c_str());
					}
				}
				
			}
			break;
		}
	}
	if (check != 1) {
		sprintf(err_msg, "Could not locate box specification in %s", filename);
		flag_error (err_msg, __FILE__, __LINE__);
		return -1;
	}
	for (int i = 0; i < 3; ++i) {
		if (box[i] <= 0.0) {
			sprintf(err_msg, "Could not read box dimensions (%g,%g,%g) properly from %s", box[0], box[1], box[2], filename);
			flag_error (err_msg, __FILE__, __LINE__);
			return -1;
		}
	}
	
	rewind(fp1);
	// also read bonds now
	check = 0;
	int nbonds = 0;
	while (fgets(buff, buffsize, fp1) != NULL) {
		if (strstr(buff, "bond") != NULL) {
			check = 1;
			vector <string> fields;
			split(fields, buff, is_any_of("=,\", "), token_compress_on);
			for (int i = 0; i < fields.size(); ++i) {
				if (fields[i] == "num") {
					if (i+1 < fields.size()) {
						nbonds = atoi(fields[i+1].c_str());
					}
					break;
				}
			}
			break;
		}
	}
	if (check != 1) {
		sprintf(err_msg, "Could not locate bond information in %s", filename);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(fp1);
		return -1;
	}
	if (nbonds < 0) {
		sprintf(err_msg, "Number of bonds < 0 in %s", filename);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(fp1);
		return -1;
	}
	sprintf(err_msg, "Found %d bonds in %s", nbonds, filename);
	flag_notify (err_msg, __FILE__, __LINE__);
	
	// allocate space to temporarily store atom coordinates, velocities, etc.
	vector <double> atom_coords(3), atom_vel(3);
	double atom_mass, atom_diameter;
	vector <int> lbond (nbonds, -1);
	vector <int> rbond (nbonds, -1);
	vector <int> bond_type (nbonds, -1);
	
	rewind(fp1);
	// read in positions, remember that by default this program normalizes the box corner to 0,0,0 and hoomd is shifted
	check = 0;
	while (fgets(buff, buffsize, fp1) != NULL) {
		if (strstr(buff, "position") != NULL) {
			check = 1;
			for (int i = 0; i < natoms; ++i) {
				fgets(buff, buffsize, fp1);
				if (sscanf(buff, "%lf %lf %lf", &atom_coords[0], &atom_coords[1], &atom_coords[2]) != 3) {
					sprintf(err_msg, "Could not read atom index %d coordinates from %s", i, filename);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(fp1);
					return -1;
				}
				for (int j = 0; j < 3; ++j) {
					// hoomd format places box at origin, we normalize corner to 0,0,0
					atom_coords[j] += box[j]/2.0;
				}
				new_atoms[i].set_pos(&atom_coords);
			}
			break;
		}
	}	
	if (check != 1) {
		sprintf(err_msg, "Could not find atom coordinates properly from %s", filename);
		flag_error (err_msg, __FILE__, __LINE__);
		return -1;
	}
	sprintf(err_msg, "Read atom coordinates from %s", filename);
	flag_notify (err_msg, __FILE__, __LINE__);
	
	rewind(fp1);
	// can optionally read velocities
	check = 0;
	while (fgets(buff, buffsize, fp1) != NULL) {
		if (strstr(buff, "velocity") != NULL) {
			check = 1;
			for (int i = 0; i < natoms; ++i) {
				fgets(buff, buffsize, fp1);
				if (sscanf(buff, "%lf %lf %lf", &new_atom[i].vel[0], &new_atom[i].vel[1], &new_atom[i].vel[2]) != 3) {
					sprintf(err_msg, "Could not read atom index %d velocity from %s", i, filename);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(fp1);
					return -1;
				}
			}
			break;
		}
	}
	if (check == 0) {
		sprintf(err_msg, "No velocities found in %s", filename);
		flag_notify (err_msg, __FILE__, __LINE__);
	}
	else {
		sprintf(err_msg, "Read velocities from %s", filename);
		flag_notify (err_msg, __FILE__, __LINE__);
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	sys->set_box(box);
	
	// organize atoms to where they need to be sent
	
	// send to slaves: box, atoms
	
	
	
}

/*!
 \param [in] \*filename Name of file to initialize from
 \param [in[ rank Rank of processor this is.  Rank 0 reads, other wait to recieve information
 \param [in,out] \*sys System object to store information this processor is responsible for
 */
int initialize (const char *filename, int *rank, int *nprocs, System *sys) {
	int argc, rc;
	char *argv[];
	
	// set up MPI
	rc = MPI_Init(&argc, &argv);
	if (rc != MPI_SUCCESS) {
		printf ("Error starting MPI. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
		return MPI_FAIL;
	}
	
	MPI_Comm_size(MPI_COMM_WORLD,nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,rank);
	
	// Create MPI_ATOM datatype
	create_MPI_ATOM();
	
	// check that interactions are such that ONLY neighboring procs need to interact
	
	// handle domain decomp --> look at MPI_graph?
	
	if (rank == 0) {
		int check = read_xml(	...	);
	} else {
		int check = worker_recv_sys(	...	);
	}
	
	return 0;
}

int finalize () {
	// free atom type after we are done running
	delete_MPI_atom();
	MPI_Finalize();
	return flag;
}
#endif