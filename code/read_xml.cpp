/**
 Read HOOMD XML file format into System object
 \author Nathan A. Mahynski
 **/

#include "read_xml.h"

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
 Print atom information to an xml file. Returns 0 if successful, -1 if failure. This only operates on the node it was called from when MPI is used, 
 the rest pause and are sequentially informed to send information as needed.
 \param [in] filename Name of file to open and read.
 \param [in] nprocs Number of processors total
 \param [in] rank Rank of process calling the routine
 \param [in] \*MCOMM MPI Communicator
 \param [in,out] \*sys System object for the main node to stores its information at.
 */
int print_xml (const string filename, const int nprocs, const int rank, MPI_COMM *MCOMM, System *sys) {
	MPI_Status Stat;
	Atom *atom_vec;
	MPI_Barrier (*MCOMM);
	int status = 0;
	
	if (rank == 0) {
		sprintf(err_msg, "Printing configuration to %s", filename.c_str());
		flag_notify (err_msg, __FILE__, __LINE__);
		
		FILE *fp1 = mfopen (filename.c_str(), "w");
		int i_need_to_print = 0, i_finished;
		if (fp1 == NULL) {
			for (int i = 1; i < nprocs; ++i) {
				MPI_Send(&i_need_to_print, 1, MPI_INT, i, 0, *MCOMM);
			}
			return -1;
		}
		
		// get atoms from workers
		int worker_atoms[nprocs-1] = {0};
		i_need_to_print = 1;
		MPI_Request send_reqs[nprocs-1], recv_reqs[nprocs-1], recv_reqs2[nprocs-1];
		MPI_Status worker_stats[nprocs-1], worker_stats2[nprocs-1];
		
		for (int i = 1; i < nprocs; ++i) {
			MPI_Isend (&i_need_to_print, 1, MPI_INT, i, 0, *MCOMM, &send_reqs[i-1]);
			MPI_Irecv (&worker_atoms[i-1], 1, MPI_INT, i, i, *MCOMM, &recv_reqs[i-1]);
		}
		if (nprocs > 1) {
			MPI_Waitall (nprocs-1, recv_reqs, worker_stats);
		}
		
		int tot_atoms = sys->num_atoms();
		Atom** system_atoms;
		if (nprocs > 1) {
			system_atoms = new (Atom*) [nprocs-1];
			for (int i = 0; i < nprocs-1; ++i) {
				system_atoms[i] = new (Atom) [worker_atoms[i]];
			}
		}
		
		for (int i = 0; i < nprocs-1; ++i) {
			tot_atoms += worker_atoms[i];
			MPI_Isend (&i_need_to_print, 1, MPI_INT, i+1, 0, *MCOMM, &send_reqs[i]);
		}
		
		for (int i = 1; i < nprocs; ++i) {
			MPI_Irecv(&system_atoms[i-1], worker_atoms[i-1], MPI_ATOM, i, i, *MCOMM, &recv_reqs2[i-1]);
		}
		if (nprocs > 1) {
			MPI_Waitall (nprocs-1, recv_reqs2, worker_stats2);
		}
		
		// Now main node has all atoms, use pointers to print atoms in order
		Atom *atom_ptr[tot_atoms];
		for (int i = 0; i < sys->num_atoms(); ++i) {
			atom_ptr[sys->get_atom(i)->sys_index] = sys->get_atom(i);
		}
		for (int i = 1; i < nprocs; ++i) {
			for (int j = 0; j < worker_atoms[i-1]; ++j) {
				atom_ptr[system_atoms[i-1][j].sys_index] = &system_atoms[i-1][j];
			}
		}
		
		// print header
		vector <double> box = sys->box(), npos(3);
		fprintf(fp1, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<hoomd_xml version=\"1.4\">\n");
		fprintf(fp1, "<configuration time_step=\"0\" dimensions=\"3\" natoms=\"%d\">\n<box lx=\"%12.12g\" ly=\"%12.12g\" lz=\"%12.12g\"/>\n", tot_atoms, box[0], box[1], box[2]);
		fprintf(fp1, "<position num=\"%d\">\n", tot_atoms);
		for (int i = 0; i < tot_atoms; ++i) {
			npos = pbc (&atom_ptr[i]->pos[0], box);
			// now shift for xml convention to have origin at center of box
			for (int m = 0; m < 3; ++m) {
				fprintf(fp1, "%12.12g\t", npos[m]-box[m]/2.0);
			}
			fprintf(fp1, "\n");
		}	
		fprintf(fp1, "</position>\n");
		
		fprintf(fp1, "<velocity num=\"%d\">\n", tot_atoms);
		for (int i = 0; i < tot_atoms; ++i) {
			for (int m = 0; m < 3; ++m) {
				fprintf(fp1, "%12.12g\t", atom_ptr[i]->vel[m]);
			}
			fprintf(fp1, "\n");
		}
		fprintf(fp1, "</velocity>\n");
		
		fprintf(fp1, "<acceleration num=\"%d\">\n", tot_atoms);
		for (int i = 0; i < tot_atoms; ++i) {
			for (int m = 0; m < 3; ++m) {
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
		
		fprintf(fp1, "<type num=\"%d\">\n", natoms);
		for (int i = 0; i < tot_atoms; ++i) {
			fprintf(fp1, "%s\n", sys->atom_name[atom_ptr[i]->type].c_str());
		}	
		
		fprintf(fp1, "<bond num=\"%d\">\n", nbonds);
		pair <int, int> ibond;
		for (int i = 0; i < sys->nbonds(); ++i) {
			ibond = sys->get_bond(i);
			fprintf(fp1, "%s\t%d\t%d\n", sys->bond_name(sys->get_bond_type(i)).c_str(), ibond.first, ibond.second);
		}
		fprintf(fp1, "</bond>\n");
		fprintf(fp1, "</configuration>\n</hoomd_xml>");

		if (nprocs > 1) {
			for (int i = 0; i < nprocs-1; ++i) {
				delete [] system_atoms[i];
			}
			delete [] system_atoms;
		}
	} else {
		int print_flag, dummy, natoms = sys->natoms();
		MPI_Recv (&print_flag, 1, MPI_INT, 0, 0, *MCOMM, &Stat);
		if (print_flag == 1) {
			MPI_Send (&natoms, 1, MPI_INT, rank, 0, *MCOMM);
			MPI_Recv (&dummy, 1, MPI_INT, 0, 0, *MCOMM, &Stat);
			atom_vec = new Atom[natoms]; 
			MPI_Send (atom_vec, natoms, MPI_ATOM, rank, 0, *MCOMM);
		} else {
			status = -1;
		}
	}
				  
    MPI_Barrier (*MCOMM);
	if (rank > 0 && status == 0) {
		delete [] atom_vec;
	}
	MPI_Barrier (*MCOMM);
	return status;
}


