/**
 Read HOOMD XML file format into System object
 \author Nathan A. Mahynski
 **/

#include "read_xml.h"

/*!
 Parse an XML file to obtain atom information. Returns 0 if successful, -1 if failure. Operates in 
 a "cascade" between ranks so that each processor (if MPI is used) opens and initializes from the
 file in order. Returns 0 if successful, else integer flag for failure.  Does domain decomposition.
 \param [in] filename Name of file to open and read.
 \param [in] nprocs Number of processors total.
 \param [in] rank Rank of this process.
 \param [in,out] \*sys System object to store its information at.
 */
int initialize_from_xml (const string filename, const int nprocs, const int rank, System *sys) {
	int iread = 1, isignal, check;
	MPI_Status Stat;
	vector <double> box = sys->box();
	
	// do domain decomposition before initialization
	check = init_domain_decomp (box, nprocs, sys->widths, &sys->final_breakup);
	
	MPI_Barrier (MPI_WORLD_COMM);
	
	// Cascade of reading statements
	if (rank == 0) {
		// read first, then signal next
		check = read_xml (filename, nprocs, rank, sys);
		
		if (nprocs > 1) {
			MPI_Send (&check, 1, MPI_INT, 1, 0, MPI_WORLD_COMM);
		}
	} else {
		// wait mode for signal from rank-1
		MPI_Recv (&isignal, 1, MPI_INT, rank-1, rank-1, MPI_WORLD_COMM, &Stat);
		if (isignal == 0) {
			// recieved signal, read
			check = read_xml (filename, nprocs, rank, sys);
			
			// signal finished to next processor in line
			if (rank+1 < nprocs) {
				MPI_Send (&check, 1, MPI_INT, rank+1, rank, MPI_WORLD_COMM);
			}
		} else {
			check = isignal;
			if (rank < nprocs-1) {
				MPI_Send (&isignal, 1, MPI_INT, rank+1, rank, MPI_WORLD_COMM);
			}
		}
	}
	
	MPI_Barrier (MPI_WORLD_COMM);
	return check;
}

/*!
 Parse an XML file to obtain atom information. Initializes a System object but only stores
 information that belongs to this processor rank according to domain decomposition.
 \param [in] filename Name of file to open and read.
 \param [in] nprocs Number of processors total.
 \param [in] rank Rank of this process.
 \param [in,out] \*sys System object to store its information at.
 */
int read_xml (const string filename, const int nprocs, const int rank, System *sys) {
	FILE *input = mfopen(filename, "r");
	if (input == NULL) {
		sprintf(err_msg, "Could not initialize from %s on rank %d", filename, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		return FILE_ERROR;
	}
	
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
		sprintf(err_msg, "Could not locate natoms in %s on rank %d", filename, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(fp1);
		return FILE_ERROR;
	}
	if (natoms < 0) {
		sprintf(err_msg, "Number of atoms < 0 in %s on rank %d", filename, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(fp1);
		return ILLEGAL_VALUE;
	}
	vector<Atom> new_atoms(natoms);
	
	// read box
	rewind(fp1);
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
		sprintf(err_msg, "Could not locate box specification in %s on rank %d", filename, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		return FILE_ERROR;
	}
	for (int i = 0; i < 3; ++i) {
		if (box[i] <= 0.0) {
			sprintf(err_msg, "Could not read box dimensions (%g,%g,%g) properly from %s on rank %d", box[0], box[1], box[2], filename, rank);
			flag_error (err_msg, __FILE__, __LINE__);
			return ILLEGAL_VALUE;
		}
	}
	sys->set_box(box);
	
	// also read bonds now
	rewind(fp1);
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
		sprintf(err_msg, "Could not locate bond information in %s on rank %d", filename, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(fp1);
		return FILE_ERROR;
	}
	if (nbonds < 0) {
		sprintf(err_msg, "Number of bonds < 0 in %s on rank %d", filename, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(fp1);
		return ILLEGAL_VALUE;
	}
	
	// allocate space to temporarily store atom coordinates, velocities, etc.
	vector <double> atom_coords(3);
	vector <int> lbond (nbonds, -1);
	vector <int> rbond (nbonds, -1);
	
	// read in positions, remember that by default this program normalizes the box corner to 0,0,0 and hoomd is shifted
	rewind(fp1);
	check = 0;
	vector <int> atom_belongs;
	vector <double> tmp_pos(3);
	int processor;
	while (fgets(buff, buffsize, fp1) != NULL) {
		if (strstr(buff, "position") != NULL) {
			check = 1;
			for (int i = 0; i < natoms; ++i) {
				fgets(buff, buffsize, fp1);
				if (sscanf(buff, "%lf %lf %lf", &atom_coords[0], &atom_coords[1], &atom_coords[2]) != 3) {
					sprintf(err_msg, "Could not read atom index %d coordinates from %s on rank %d", i, filename, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(fp1);
					return FILE_ERROR;
				}
				for (int j = 0; j < 3; ++j) {
					// hoomd format places box at origin, we normalize corner to 0,0,0
					atom_coords[j] += box[j]/2.0;
					new_atoms[i].pos[j] = atom_coords[j];
					tmp_pos[j] = atom_coords[j];
				}
				
				// see if this atom belongs on this processor
				processor = get_processor (tmp_pos, sys->widths, &sys->final_breakup);
				if (processor == rank) {
					atom_belongs.push_back(i);
				}
				
				new_atoms[i].sys_index = i;
			}
			break;
		}
	}	
	if (check != 1) {
		sprintf(err_msg, "Could not find atom coordinates properly from %s on rank %d", filename, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(fp1);
		return FILE_ERROR;
	}
	
	// read velocities
	rewind(fp1);
	check = 0;
	while (fgets(buff, buffsize, fp1) != NULL) {
		if (strstr(buff, "velocity") != NULL) {
			check = 1;
			for (int i = 0; i < natoms; ++i) {
				fgets(buff, buffsize, fp1);
				if (sscanf(buff, "%lf %lf %lf", &new_atom[i].vel[0], &new_atom[i].vel[1], &new_atom[i].vel[2]) != 3) {
					sprintf(err_msg, "Could not read atom index %d velocity from %s on rank %d", i, filename, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(fp1);
					return FILE_ERROR;
				}
			}
			break;
		}
	}
	if (check == 0) {
		sprintf(err_msg, "No velocities found in %s on rank %d", filename, rank);
		flag_notify (err_msg, __FILE__, __LINE__);
		fclose(fp1);
		return FILE_ERROR;
	}
	
	// read mass
	rewind(fp1);
	check = 0;
	while (fgets(buff, buffsize, fp1) != NULL) {
		if (strstr(buff, "mass") != NULL) {
			check = 1;
			for (int i = 0; i < natoms; ++i) {
				fgets(buff, buffsize, fp1);
				if (sscanf(buff, "%lf", &new_atoms[i].mass) != 1) {
					sprintf(err_msg, "Could not read atom index %d mass from %s on rank %d", i, filename, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(fp1);
					return ILLEGAL_VALUE;
				}
				if (new_atoms[i].mass < 0.0) {
					sprintf(err_msg, "Atom index %d mass = %g < 0.0 in %s on rank %d", i, new_atoms[i].mass, filename, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(fp1);
					return ILLEGAL_VALUE;
				}
			}
			break;
		}
	}
	if (check == 0) {
		sprintf(err_msg, "No masses found in %s on rank %d", filename, rank);
		flag_notify (err_msg, __FILE__, __LINE__);
		fclose(fp1);
		return FILE_ERROR;
	}
	sprintf(err_msg, "Read atom masses from %s on rank %d", filename, rank);
	flag_notify (err_msg, __FILE__, __LINE__);
	
	// read sizes (diameters)
	rewind(fp1);
	check = 0;
	while (fgets(buff, buffsize, fp1) != NULL) {
		if (strstr(buff, "diameter") != NULL) {
			check = 1;
			for (int i = 0; i < natoms; ++i) {
				fgets(buff, buffsize, fp1);
				if (sscanf(buff, "%lf", &new_atoms[i].diam) != 1) {
					sprintf(err_msg, "Could not read atom index %d diameter from %s on rank %d", i, filename, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(fp1);
					return FILE_ERROR;
				}
				if (new_atoms[i].diam < 0.0) {
					sprintf(err_msg, "Atom index %d diameter = %g < 0.0 in %s on rank %d", i, new_atoms[i].diam, filename, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(fp1);
					return ILLEGAL_VALUE;
				}
			}
			break;
		}
	}
	if (check == 0) {
		sprintf(err_msg, "No atom diameters found in %s on rank %d", filename, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(fp1);
		return FILE_ERROR;
	}

	// read types
	rewind(fp1);
	char aname[ATOM_NAME_LENGTH];
	check = 0;
	while (fgets(buff, buffsize, fp1) != NULL) {
		if (strstr(buff, "type") != NULL) {
			check = 1;
			for (int i = 0; i < natoms; ++i) {
				fgets(buff, buffsize, fp1);
				if (sscanf(buff, "%s", aname) != 1) {
					sprintf(err_msg, "Could not read atom index %d type name from %s on rank %d", i, filename, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(fp1);
					return FILE_ERROR;
				}
				if (sys->add_atom_type(aname) == -1) {
					sprintf(err_msg, "Could not convert atom index %d type name from %s to an internal index on rank %d", i, filename, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(fp1);
					return ILLEGAL_VALUE;
				}
				else {
					new_atoms[i].type = sys->atom_type(aname);
				}
			}
			break;
		}
	}
	if (check == 0) {
		sprintf(err_msg, "No atom types found in %s on rank %d", filename, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(fp1);
		return FILE_ERROR;
	}

	// read bonds
	rewind(fp1);
	char bname[BOND_NAME_LENGTH];
	check = 0;
	while (fgets(buff, buffsize, fp1) != NULL) {
		if (strstr(buff, "bond") != NULL) {
			check = 1;
			for (int i = 0; i < nbonds; ++i) {
				fgets(buff, buffsize, fp1);
				if (sscanf(buff, "%s %d %d", bname, &lbond[i], &rbond[i]) != 3) {
					sprintf(err_msg, "Could not read bond number %d from %s on rank %d", i+1, filename, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(fp1);
					return FILE_ERROR;
				}
				if (lbond[i] < 0 || rbond[i] < 0) {
					sprintf(err_msg, "Bond number %d is between illegal values for atoms (%d,%d) in %s on rank %d", i+1, lbond[i], rbond[i], filename, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(fp1);
					return ILLEGAL_VALUE;
				}
				if (lbond[i] == rbond[i]) {
					sprintf(err_msg, "Atom number %d appears bonded to itself in %s on rank %d", i+1, filename, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(fp1);
					return ILLEGAL_VALUE;
				}
				if (sys->add_bond_type(bname) == -1) {
					sprintf(err_msg, "Could not convert bond index %d type name from %s to an internal index on rank %d", i, filename, rank);
					flag_error (err_msg, __FILE__, __LINE__);
					fclose(fp1);
					return ILLEGAL_VALUE;
				} else {
					sys->add_bond(lbond, rbond, sys->bond_type(bname));
				}	
			}
			break;
		}
	}
	if (check == 0) {
		sprintf(err_msg, "No bonds found in %s on rank %d", filename, rank);
		flag_error (err_msg, __FILE__, __LINE__);
		fclose(fp1);
		return FILE_ERROR;
	}

	// now add the atoms that belong to this domain to the System object
	Atom atom_array[(const int)atom_belongs.size()];
	for (int i = 0; i < atom_belongs.size(); ++i) {
		atom_array[i] = new_atoms[atom_belongs[i]];
	}
	sys->add_atoms(atom_belongs.size(), atom_array);
	
	return 0;
}

/*!
 Print atom information to an xml file. Returns 0 if successful, -1 if failure. This only operates on the node it was called from when MPI is used, 
 the rest pause and are sequentially informed to send information as needed.
 \param [in] filename Name of file to open and read.
 \param [in] nprocs Number of processors total
 \param [in] rank Rank of process calling the routine
 \param [in,out] \*sys System object for the main node to stores its information at.
 */
int print_xml (const string filename, const int nprocs, const int rank, const System *sys) {
	MPI_Status Stat;
	Atom *atom_vec;
	MPI_Barrier (MPI_WORLD_COMM);
	int status = 0;
	
	if (rank == 0) {
		sprintf(err_msg, "Printing configuration to %s", filename.c_str());
		flag_notify (err_msg, __FILE__, __LINE__);
		
		FILE *fp1 = mfopen (filename.c_str(), "w");
		int i_need_to_print = 0, i_finished;
		if (fp1 == NULL) {
			for (int i = 1; i < nprocs; ++i) {
				MPI_Send(&i_need_to_print, 1, MPI_INT, i, 0, MPI_WORLD_COMM);
			}
			return -1;
		}
		
		// get atoms from workers
		int worker_atoms[nprocs-1] = {0};
		i_need_to_print = 1;
		MPI_Request send_reqs[nprocs-1], recv_reqs[nprocs-1], recv_reqs2[nprocs-1];
		MPI_Status worker_stats[nprocs-1], worker_stats2[nprocs-1];
		
		for (int i = 1; i < nprocs; ++i) {
			MPI_Isend (&i_need_to_print, 1, MPI_INT, i, 0, MPI_WORLD_COMM, &send_reqs[i-1]);
			MPI_Irecv (&worker_atoms[i-1], 1, MPI_INT, i, i, MPI_WORLD_COMM, &recv_reqs[i-1]);
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
			MPI_Isend (&i_need_to_print, 1, MPI_INT, i+1, 0, MPI_WORLD_COMM, &send_reqs[i]);
		}
		
		for (int i = 1; i < nprocs; ++i) {
			MPI_Irecv(&system_atoms[i-1], worker_atoms[i-1], MPI_ATOM, i, i, MPI_WORLD_COMM, &recv_reqs2[i-1]);
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
		MPI_Recv (&print_flag, 1, MPI_INT, 0, 0, MPI_WORLD_COMM, &Stat);
		if (print_flag == 1) {
			MPI_Send (&natoms, 1, MPI_INT, rank, 0, MPI_WORLD_COMM);
			MPI_Recv (&dummy, 1, MPI_INT, 0, 0, MPI_WORLD_COMM, &Stat);
			atom_vec = new Atom[natoms]; 
			MPI_Send (atom_vec, natoms, MPI_ATOM, rank, 0, MPI_WORLD_COMM);
		} else {
			status = -1;
		}
	}
				  
    MPI_Barrier (MPI_WORLD_COMM);
	if (rank > 0 && status == 0) {
		delete [] atom_vec;
	}
	MPI_Barrier (MPI_WORLD_COMM);
	return status;
}


