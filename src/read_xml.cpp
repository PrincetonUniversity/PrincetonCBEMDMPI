/**
 \file read_xml.cpp
 \brief I/O for XML file format 
 \author Nathan A. Mahynski
 **/

#include "read_xml.h"

/*!
 Parse an XML file to obtain atom information. Initializes a System object but only stores
 information that belongs to this processor rank according to domain decomposition.
 Returns SAFE_EXIT if successful, else returns an error flag.
 \param [in] filename Name of file to open and read.
 \param [in,out] \*sys Pointer to System object to store its information at.
 */
int read_xml (const string filename, System *sys) {
	int nprocs, rank;
	MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	
	const char *filename_cstr=filename.c_str();
	char err_msg[MYERR_FLAG_SIZE];
	FILE *input = mfopen(filename_cstr, "r");
	if (input == NULL) {
		sprintf(err_msg, "Could not initialize from %s on rank %d", filename_cstr, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		return FILE_ERROR;
	}
	
	// Read header
	int check = 0;
	int natoms = -1;
	const int buffsize=1000;
	char buff[buffsize], *dummy_char, *dummy_char2;
	while ((dummy_char = fgets(buff, buffsize, input)) != NULL) {
		if (strstr(buff, "natoms") != NULL) {
			vector <string> fields;
			split(fields, buff, is_any_of("=,\", "), token_compress_on);
			for (unsigned int i = 0; i < fields.size(); ++i) {
				if (fields[i] == "natoms") {
					if (i+1 < fields.size()) {
						natoms = atoi(fields[i+1].c_str());
						check = 1;
					}
					break;
				}
			}
			break;
		}
	}
	if (check != 1) {
		sprintf(err_msg, "Could not locate natoms in %s on rank %d", filename_cstr, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(input);
		return FILE_ERROR;
	}
	if (natoms <= 0) {
		sprintf(err_msg, "Number of atoms < 0 in %s on rank %d", filename_cstr, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(input);
		return ILLEGAL_VALUE;
	}
	vector<Atom> new_atoms(natoms);
	
	try {
		sys->global_atom_types.resize(natoms);
	}
	catch (bad_alloc& ba) {
		sprintf(err_msg, "Out of memory");
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(input);
		return BAD_MEM;
	}

	// Read box size
	rewind(input);
	check = 0;
	vector<double> box(NDIM, -1);
	while ((dummy_char = fgets(buff, buffsize, input)) != NULL) {
		if (strstr(buff, "box") != NULL) {
			check = 1;
			vector <string> fields;
			split(fields, buff, is_any_of("=,\", "), token_compress_on);
			for (unsigned int i = 0; i < fields.size(); ++i) {
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
		sprintf(err_msg, "Could not locate box specification in %s on rank %d", filename_cstr, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(input);
		return FILE_ERROR;
	}
	for (int i = 0; i < NDIM; ++i) {
		if (box[i] <= 0.0) {
			sprintf(err_msg, "Could not read box dimensions (%g,%g,%g) properly from %s on rank %d", box[0], box[1], box[2], filename_cstr, rank);
			flag_error (err_msg, __FILE__, __LINE__);
			fclose (input);
			return ILLEGAL_VALUE;
		}
	}
	sys->set_box(box);

	// Read bond information
	rewind(input);
	check = 0;
	int nbonds = 0;
	while ((dummy_char = fgets(buff, buffsize, input)) != NULL) {
		if (strstr(buff, "bond") != NULL) {
			vector <string> fields;
			split(fields, buff, is_any_of("=,\", "), token_compress_on);
			for (unsigned int i = 0; i < fields.size(); ++i) {
				if (fields[i] == "num") {
					if (i+1 < fields.size()) {
						nbonds = atoi(fields[i+1].c_str());
						check = 1;
					}
					break;
				}
			}
			break;
		}
	}
	if (check != 1) {
		sprintf(err_msg, "Could not locate bond information in %s on rank %d", filename_cstr, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(input);
		return FILE_ERROR;
	}
	if (nbonds < 0) {
		sprintf(err_msg, "Number of bonds < 0 in %s on rank %d", filename_cstr, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(input);
		return ILLEGAL_VALUE;
	}
	
	// Allocate space to temporarily store atom coordinates, velocities, etc.
	vector <double> atom_coords(NDIM);
	vector <int> lbond (nbonds, -1);
	vector <int> rbond (nbonds, -1);
	
	// Read in positions
	rewind(input);
	check = 0;
	vector <int> atom_belongs;
	while ((dummy_char = fgets(buff, buffsize, input)) != NULL) {
		if (strstr(buff, "position") != NULL) {
			check = 1;
			for (int i = 0; i < natoms; ++i) {
			    if ((dummy_char2 = fgets(buff, buffsize, input)) == NULL) {
				sprintf(err_msg, "Unexpected read error on rank %d", rank);
				flag_error (err_msg, __FILE__, __LINE__);
				fclose(input);
				return FILE_ERROR;
			    }				
				if (sscanf(buff, "%lf %lf %lf", &atom_coords[0], &atom_coords[1], &atom_coords[2]) != NDIM) {
					sprintf(err_msg, "Could not read atom index %d coordinates from %s on rank %d", i, filename_cstr, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(input);
					return FILE_ERROR;
				}
				for (int j = 0; j < NDIM; ++j) {
					if (atom_coords[j] < 0.0 || atom_coords[j] >= box[j]) {
						sprintf(err_msg, "Atom index %d coordinates (%g,%g,%g) out of range", i, atom_coords[0], atom_coords[1], atom_coords[2]);
						flag_error (err_msg, __FILE__, __LINE__);
						return ILLEGAL_VALUE;
					}
					new_atoms[i].pos[j] = atom_coords[j];
				}

				// See if this atom belongs on this processor
				if (floor(new_atoms[i].pos[PARALLELDIM] / sys->box()[PARALLELDIM] * nprocs) == rank) {
					atom_belongs.push_back(i);
				}
				
				new_atoms[i].sys_index = i;
			}
			break;
		}
	}	
	if (check != 1) {
		sprintf(err_msg, "Could not find atom coordinates properly from %s on rank %d", filename_cstr, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(input);
		return FILE_ERROR;
	}
	
	// Read velocities
	rewind(input);
	check = 0;
	while ((dummy_char = fgets(buff, buffsize, input)) != NULL) {
		if (strstr(buff, "velocity") != NULL) {
			check = 1;
			for (int i = 0; i < natoms; ++i) {
			    if ((dummy_char2 = fgets(buff, buffsize, input)) == NULL) {
				sprintf(err_msg, "Unexpected read error on rank %d", rank);
				flag_error (err_msg, __FILE__, __LINE__);
				fclose(input);
				return FILE_ERROR;
			    }				
				if (sscanf(buff, "%lf %lf %lf", &new_atoms[i].vel[0], &new_atoms[i].vel[1], &new_atoms[i].vel[2]) != NDIM) {
					sprintf(err_msg, "Could not read atom index %d velocity from %s on rank %d", i, filename_cstr, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(input);
					return FILE_ERROR;
				}
			}
			break;
		}
	}
	if (check != 1) {
		sprintf(err_msg, "No velocities found in %s on rank %d", filename_cstr, rank);
		flag_notify (err_msg, __FILE__, __LINE__);
		fclose(input);
		return FILE_ERROR;
	}
	
	// Read mass
	rewind(input);
	check = 0;
	while ((dummy_char = fgets(buff, buffsize, input)) != NULL) {
		if (strstr(buff, "mass") != NULL) {
			check = 1;
			for (int i = 0; i < natoms; ++i) {
			    if ((dummy_char2 = fgets(buff, buffsize, input)) == NULL) {
				sprintf(err_msg, "Unexpected read error on rank %d", rank);
				flag_error (err_msg, __FILE__, __LINE__);
				fclose(input);
				return FILE_ERROR;
			    }				
				if (sscanf(buff, "%lf", &new_atoms[i].mass) != 1) {
					sprintf(err_msg, "Could not read atom index %d mass from %s on rank %d", i, filename_cstr, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(input);
					return ILLEGAL_VALUE;
				}
				if (new_atoms[i].mass < 0.0) {
					sprintf(err_msg, "Atom index %d mass = %g < 0.0 in %s on rank %d", i, new_atoms[i].mass, filename_cstr, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(input);
					return ILLEGAL_VALUE;
				}
			}
			break;
		}
	}
	if (check != 1) {
		sprintf(err_msg, "No masses found in %s on rank %d", filename_cstr, rank);
		flag_notify (err_msg, __FILE__, __LINE__);
		fclose(input);
		return FILE_ERROR;
	}
	
	// Read sizes (diameters)
	rewind(input);
	check = 0;
	while ((dummy_char = fgets(buff, buffsize, input)) != NULL) {
		if (strstr(buff, "diameter") != NULL) {
			check = 1;
			for (int i = 0; i < natoms; ++i) {
			    if ((dummy_char2 = fgets(buff, buffsize, input)) == NULL) {
				sprintf(err_msg, "Unexpected read error on rank %d", rank);
				flag_error (err_msg, __FILE__, __LINE__);
				fclose(input);
				return FILE_ERROR;
			    }				
				if (sscanf(buff, "%lf", &new_atoms[i].diam) != 1) {
					sprintf(err_msg, "Could not read atom index %d diameter from %s on rank %d", i, filename_cstr, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(input);
					return FILE_ERROR;
				}
				if (new_atoms[i].diam < 0.0) {
					sprintf(err_msg, "Atom index %d diameter = %g < 0.0 in %s on rank %d", i, new_atoms[i].diam, filename_cstr, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(input);
					return ILLEGAL_VALUE;
				}
			}
			break;
		}
	}
	if (check != 1) {
		sprintf(err_msg, "No atom diameters found in %s on rank %d", filename_cstr, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(input);
		return FILE_ERROR;
	}

	// Read types
	rewind(input);
	char aname[ATOM_NAME_LENGTH];
	check = 0;
	while ((dummy_char = fgets(buff, buffsize, input)) != NULL) {
		if (strstr(buff, "type") != NULL) {
			check = 1;
			for (int i = 0; i < natoms; ++i) {
			    if ((dummy_char2 = fgets(buff, buffsize, input)) == NULL) {
				sprintf(err_msg, "Unexpected read error on rank %d", rank);
				flag_error (err_msg, __FILE__, __LINE__);
				fclose(input);
				return FILE_ERROR;
			    }				
				if (sscanf(buff, "%s", aname) != 1) {
					sprintf(err_msg, "Could not read atom index %d type name from %s on rank %d", i, filename_cstr, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(input);
					return FILE_ERROR;
				}
				if (sys->add_atom_type(aname) == -1) {
					sprintf(err_msg, "Could not convert atom index %d type name from %s to an internal index on rank %d", i, filename_cstr, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(input);
					return ILLEGAL_VALUE;
				}
				else {
					new_atoms[i].type = sys->atom_type(aname);
					sys->global_atom_types[i] = aname;
				}
			}
			break;
		}
	}
	if (check != 1) {
		sprintf(err_msg, "No atom types found in %s on rank %d", filename_cstr, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(input);
		return FILE_ERROR;
	}

	// Read bonds
	rewind(input);
	char bname[BOND_NAME_LENGTH];
	check = 0;
	while ((dummy_char = fgets(buff, buffsize, input)) != NULL) {
		if (strstr(buff, "bond") != NULL) {
			check = 1;
			for (int i = 0; i < nbonds; ++i) {
			    if ((dummy_char2 = fgets(buff, buffsize, input)) == NULL) {
				sprintf(err_msg, "Unexpected read error on rank %d", rank);
				flag_error (err_msg, __FILE__, __LINE__);
				fclose(input);
				return FILE_ERROR;
			    }				
				if (sscanf(buff, "%s %d %d", bname, &lbond[i], &rbond[i]) != NDIM) {
					sprintf(err_msg, "Could not read bond number %d from %s on rank %d", i+1, filename_cstr, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(input);
					return FILE_ERROR;
				}
				if (lbond[i] < 0 || rbond[i] < 0) {
					sprintf(err_msg, "Bond number %d is between illegal values for atoms (%d,%d) in %s on rank %d", i+1, lbond[i], rbond[i], filename_cstr, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(input);
					return ILLEGAL_VALUE;
				}
				if (lbond[i] == rbond[i]) {
					sprintf(err_msg, "Atom number %d appears bonded to itself in %s on rank %d", i+1, filename_cstr, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(input);
					return ILLEGAL_VALUE;
				}
				if (sys->add_bond_type(bname) == -1) {
					sprintf(err_msg, "Could not convert bond index %d type name from %s to an internal index on rank %d", i, filename_cstr, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(input);
					return ILLEGAL_VALUE;
				} else {
					sys->add_bond(lbond[i], rbond[i], sys->bond_type(bname));
				}	
			}
			break;
		}
	}
	if (check != 1) {
		sprintf(err_msg, "No bond information found in %s on rank %d", filename_cstr, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(input);
		return FILE_ERROR;
	}

	// Now add the atoms that belong to this domain to the System object
	vector<Atom> atom_array;
	atom_array.resize(atom_belongs.size());
	for (unsigned int i = 0; i < atom_belongs.size(); ++i) {
		atom_array[i] = new_atoms[atom_belongs[i]];
	}
	vector <int> dummy_vec = sys->add_atoms(&atom_array);
	sys->set_num_atoms(atom_belongs.size());
	
	sprintf(err_msg, "Successfully read coordinates from %s on rank %d", filename_cstr, rank);
	flag_notify (err_msg, __FILE__, __LINE__);
	return SAFE_EXIT;
}

/*!
 Print atom information to an xml file. Returns SAFE_EXIT if successful, else an error flag if failure. This only operates on the master node when multiple processors are used, 
 the rest pause and are sequentially informed to send information as needed.
 \param [in] filename Name of file to open and read.
 \param [in,out] \*sys Pointer to System object for the main node to stores its information at.
 */
int print_xml (const string filename, const System *sys) {
    char err_msg[MYERR_FLAG_SIZE];
	MPI_Status Stat;
	Atom *atom_vec;
	int status = 0, rank, nprocs;
	
	MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Barrier (MPI_COMM_WORLD);
	
	if (rank == 0) {
		sprintf(err_msg, "Printing configuration to %s", filename.c_str());
		flag_notify (err_msg, __FILE__, __LINE__);
		
		FILE *fp1 = mfopen (filename.c_str(), "w");
		int i_need_to_print = 0;
		if (fp1 == NULL) {
			for (int i = 1; i < nprocs; ++i) {
				MPI_Send(&i_need_to_print, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			}
			return -1;
		}

		// Get number of incoming atoms from workers
		const int nprocs_less_one = nprocs-1;
		int worker_atoms[nprocs_less_one];
		for (int i=0; i<nprocs-1; i++) {
		    worker_atoms[i] = 0;
		}
		i_need_to_print = 1;
		MPI_Request send_reqs[nprocs_less_one], recv_reqs[nprocs_less_one], recv_reqs2[nprocs_less_one];
		MPI_Status worker_stats[nprocs_less_one], worker_stats2[nprocs_less_one];

		for (int i = 1; i < nprocs; ++i) {
			MPI_Isend (&i_need_to_print, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &send_reqs[i-1]);
			MPI_Irecv (&worker_atoms[i-1], 1, MPI_INT, i, i, MPI_COMM_WORLD, &recv_reqs[i-1]);
		}
		if (nprocs > 1) {
			MPI_Waitall (nprocs-1, recv_reqs, worker_stats);
		}

		int tot_atoms = sys->natoms();

		// Get memory "offsets" so we know where to store incoming atoms
		vector <int> offsets (nprocs-1);
		for (int i = 0; i < nprocs-1; ++i) {
			tot_atoms += worker_atoms[i];
			if (i > 0) {
				offsets[i] = worker_atoms[i-1]+offsets[i-1];
			} else {
				offsets[i] = 0;
			}
			MPI_Isend (&i_need_to_print, 1, MPI_INT, i+1, 0, MPI_COMM_WORLD, &send_reqs[i]);
		}

		// Allocate room for incoming atoms
		Atom *system_atoms = new Atom[tot_atoms-sys->natoms()];
		
		for (int i = 0; i < nprocs-1; ++i) {
			MPI_Irecv(&system_atoms[offsets[i]], worker_atoms[i], MPI_ATOM, i+1, i+1, MPI_COMM_WORLD, &recv_reqs2[i]);
		}

		if (nprocs > 1) {
			MPI_Waitall (nprocs-1, recv_reqs2, worker_stats2);
		}

		// Now main node has all atoms, use pointers to print atoms in order
		vector <Atom *> atom_ptr(tot_atoms);
		for (int i = 0; i < sys->natoms(); ++i) {
		    atom_ptr[((System *)sys)->get_atom(i)->sys_index] = ((System *)sys)->get_atom(i);
		}

		for (int i = 0; i < tot_atoms-sys->natoms(); ++i) {
			atom_ptr[system_atoms[i].sys_index] = &system_atoms[i];
		}

		// Print the file
		vector <double> box = sys->box(), npos(NDIM);
		fprintf(fp1, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<hoomd_xml version=\"1.4\">\n");
		fprintf(fp1, "<configuration time_step=\"0\" dimensions=\"3\" natoms=\"%d\">\n<box lx=\"%12.12g\" ly=\"%12.12g\" lz=\"%12.12g\"/>\n", tot_atoms, box[0], box[1], box[2]);
		fprintf(fp1, "<position num=\"%d\">\n", tot_atoms);
		for (int i = 0; i < tot_atoms; ++i) {
			npos = pbc (&atom_ptr[i]->pos[0], box);
			for (int m = 0; m < NDIM; ++m) {
				if (npos[m] != npos[m]) {
					sprintf(err_msg, "Atom %d has nan for its coordinates, cannot print", i);
					flag_error (err_msg, __FILE__, __LINE__);
					return ILLEGAL_VALUE;
				}
				fprintf(fp1, "%12.12g\t", npos[m]);
			}
			fprintf(fp1, "\n");
		}	
		fprintf(fp1, "</position>\n");

		fprintf(fp1, "<velocity num=\"%d\">\n", tot_atoms);
		for (int i = 0; i < tot_atoms; ++i) {
			for (int m = 0; m < NDIM; ++m) {
				fprintf(fp1, "%12.12g\t", atom_ptr[i]->vel[m]);
			}
			fprintf(fp1, "\n");
		}
		fprintf(fp1, "</velocity>\n");
		
		fprintf(fp1, "<acceleration num=\"%d\">\n", tot_atoms);
		for (int i = 0; i < tot_atoms; ++i) {
			for (int m = 0; m < NDIM; ++m) {
				fprintf(fp1, "%12.12g\t", atom_ptr[i]->force[m]/atom_ptr[i]->mass);
			}
			fprintf(fp1, "\n");
		}	
		fprintf(fp1, "</acceleration>\n");
		
		fprintf(fp1, "<mass num=\"%d\">\n", tot_atoms);
		for (int i = 0; i < tot_atoms; ++i) {
			fprintf(fp1, "%12.12g\n", atom_ptr[i]->mass);
		}	
		fprintf(fp1, "</mass>\n");
		
		fprintf(fp1, "<diameter num=\"%d\">\n", tot_atoms);
		for (int i = 0; i < tot_atoms; ++i) {
			fprintf(fp1, "%12.12g\n", atom_ptr[i]->diam);
		}	
		fprintf(fp1, "</diameter>\n");
		
		fprintf(fp1, "<type num=\"%d\">\n", tot_atoms);
		for (int i = 0; i < tot_atoms; ++i) {
			fprintf(fp1, "%s\n", ((System *)sys)->atom_name(atom_ptr[i]->type).c_str());
		}	
		
		fprintf(fp1, "<bond num=\"%d\">\n", ((System *)sys)->nbonds());
		pair <int, int> ibond;
		for (int i = 0; i < ((System *)sys)->nbonds(); ++i) {
			ibond = ((System *)sys)->get_bond(i);
			fprintf(fp1, "%s\t%d\t%d\n", ((System *)sys)->bond_name(((System *)sys)->get_bond_type(i)).c_str(), ibond.first, ibond.second);
		}
		fprintf(fp1, "</bond>\n");
		fprintf(fp1, "</configuration>\n</hoomd_xml>");

		if (nprocs > 1) {
			delete [] system_atoms;
		}
	} else {
		// Wait to recieve a signal to report atoms on this processor
		int print_flag, dummy, natoms = sys->natoms();
		MPI_Recv (&print_flag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &Stat);
		if (print_flag == 1) {
			MPI_Send (&natoms, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
			MPI_Recv (&dummy, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &Stat);
			atom_vec = new Atom[natoms]; 
			for (int i = 0; i < natoms; ++i) {
				atom_vec[i] = ((System *)sys)->copy_atom(i);
			}
			MPI_Send (atom_vec, natoms, MPI_ATOM, 0, rank, MPI_COMM_WORLD);
			delete [] atom_vec;
		} else {
			// If there is a problem, report that to all subsequent processors and fail
			sprintf(err_msg, "Rank %d received bad signal to report its atom to master node, print failure", rank);
			flag_error (err_msg, __FILE__, __LINE__);
			status = FILE_ERROR;
		}
	}
	
	MPI_Barrier (MPI_COMM_WORLD);
	return status;
}

/*!
 Write an xyz file that stores the coordinates of all the atoms in the simulation
 If timestep=0, then the file is created, or overwritten if it exists already.
 If timestep>0, then the current information is appended to the existing file.
 This file can then be read by VMD or another visualization program to produce animations.
 \param [in] filename Name of file to open and write to.
 \param [in] \*sys System object where the atoms are stored
 \param [in] timestep Current timestep of the simulation
 \param [in] wrap_pos Whether the positions should be written in unwrapped (wrap_pos=false) or wrapped (wrap_pos=true) coordinates
 */
int write_xyz (const string filename, const System *sys, const int timestep, const bool wrap_pos) {
	char err_msg[MYERR_FLAG_SIZE];
	MPI_Status Stat;
	Atom *atom_vec;
	int status = 0, rank, nprocs;
	
	MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Barrier (MPI_COMM_WORLD);
	
	if (rank == 0) {
		FILE *fp1;
		if (timestep == 0) {
			fp1 = mfopen (filename.c_str(), "w");
		}
		else {
			fp1 = mfopen (filename.c_str(), "a");
		}
		int i_need_to_print = 0;
		if (fp1 == NULL) {
			for (int i = 1; i < nprocs; ++i) {
				MPI_Send(&i_need_to_print, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			}
			fclose(fp1);
			return -1;
		}

		// Get number of incoming atoms from workers
		const int nprocs_less_one = nprocs-1;
		int worker_atoms[nprocs_less_one];
		for (int i=0; i<nprocs-1; i++) {
		    worker_atoms[i] = 0;
		}
		i_need_to_print = 1;
		MPI_Request send_reqs[nprocs_less_one], recv_reqs[nprocs_less_one], recv_reqs2[nprocs_less_one];
		MPI_Status worker_stats[nprocs_less_one], worker_stats2[nprocs_less_one];

		for (int i = 1; i < nprocs; ++i) {
			MPI_Isend (&i_need_to_print, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &send_reqs[i-1]);
			MPI_Irecv (&worker_atoms[i-1], 1, MPI_INT, i, i, MPI_COMM_WORLD, &recv_reqs[i-1]);
		}
		if (nprocs > 1) {
			MPI_Waitall (nprocs-1, recv_reqs, worker_stats);
		}

		int tot_atoms = sys->natoms();

		// Get memory "offsets" so we know where to store incoming atoms
		vector <int> offsets (nprocs-1);
		for (int i = 0; i < nprocs-1; ++i) {
			tot_atoms += worker_atoms[i];
			if (i > 0) {
				offsets[i] = worker_atoms[i-1]+offsets[i-1];
			} else {
				offsets[i] = 0;
			}
			MPI_Isend (&i_need_to_print, 1, MPI_INT, i+1, 0, MPI_COMM_WORLD, &send_reqs[i]);
		}

		// Allocate room for incoming atoms
		Atom *system_atoms = new Atom[tot_atoms-sys->natoms()];
		
		for (int i = 0; i < nprocs-1; ++i) {
			MPI_Irecv(&system_atoms[offsets[i]], worker_atoms[i], MPI_ATOM, i+1, i+1, MPI_COMM_WORLD, &recv_reqs2[i]);
		}

		if (nprocs > 1) {
			MPI_Waitall (nprocs-1, recv_reqs2, worker_stats2);
		}

		// Now main node has all atoms, use pointers to print atoms in order
		vector <Atom *> atom_ptr(tot_atoms);
		for (int i = 0; i < sys->natoms(); ++i) {
		    atom_ptr[((System *)sys)->get_atom(i)->sys_index] = ((System *)sys)->get_atom(i);
		}

		for (int i = 0; i < tot_atoms-sys->natoms(); ++i) {
			atom_ptr[system_atoms[i].sys_index] = &system_atoms[i];
		}

		// Print record
		vector <double> box = sys->box(), npos(NDIM);
		fprintf(fp1, "%d\n", tot_atoms);
		fprintf(fp1, "Atoms. Timestep: %d\n", timestep);
		for (int i = 0; i < tot_atoms; ++i) {
			fprintf(fp1, "%s ", ((System *)sys)->atom_name(atom_ptr[i]->type).c_str());
			if (wrap_pos) {
				npos = pbc (&atom_ptr[i]->pos[0], box);
			} else {
				npos.assign (atom_ptr[i]->pos, atom_ptr[i]->pos + sizeof(atom_ptr[i]->pos) / sizeof(double)); 
			}
			for (int m = 0; m < NDIM; ++m) {
				if (npos[m] != npos[m]) {
					sprintf(err_msg, "Atom %d has nan for its coordinates, cannot print", i);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(fp1);
					return ILLEGAL_VALUE;
				}
				fprintf(fp1, "%g", npos[m]);
				if (m < 2) {
					fprintf(fp1, " ");
				}

			}
			fprintf(fp1, "\n");
		}	
		if (nprocs > 1) {
			delete [] system_atoms;
		}
		fclose(fp1);

	} else {
		// Wait to report atoms to main node
		int print_flag, dummy, natoms = sys->natoms();
		MPI_Recv (&print_flag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &Stat);
		if (print_flag == 1) {
			MPI_Send (&natoms, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
			MPI_Recv (&dummy, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &Stat);
			atom_vec = new Atom[natoms]; 
			for (int i = 0; i < natoms; ++i) {
				atom_vec[i] = ((System *)sys)->copy_atom(i);
			}
			MPI_Send (atom_vec, natoms, MPI_ATOM, 0, rank, MPI_COMM_WORLD);
			delete [] atom_vec;
		} else {
			// If failure, report failure to subsequent processes
			sprintf(err_msg, "Rank %d received bad signal to report its atom to master node, print failure", rank);
			flag_error (err_msg, __FILE__, __LINE__);
			status = FILE_ERROR;
		}
	}

	MPI_Barrier (MPI_COMM_WORLD);
	return status;
}
