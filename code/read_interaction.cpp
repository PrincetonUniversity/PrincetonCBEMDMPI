#include "read_interaction.h"

int read_interactions (const string filename, System *sys) {
	const int buffsize = 1000;
	char buff[buffsize];
	double rcut_max = -1.0;
	
	char err_msg[MYERR_FLAG_SIZE];
	const char *filename_cstr=filename.c_str();
	FILE *input = mfopen(filename_cstr, "r");
	if (input == NULL) {
		sprintf(err_msg, "Could not initialize from %s", filename_cstr);
		flag_error (err_msg, __FILE__, __LINE__);
		return FILE_ERROR;
	}
	
	vector <pair <string, string> > atom_pairs;
	vector <Interaction> inters_PPOT, inters_BOND;
	vector <string> bond_list;
	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	while (fgets(buff, buffsize, input) != NULL) {
		vector <string> fields;
		split(fields, buff, is_any_of("=,\", "), token_compress_on);
		
		if (fields.size() < 1) {
			sprintf(err_msg, "Insufficient data");
			flag_error (err_msg, __FILE__, __LINE__);
			return FILE_ERROR;
		}

		if (fields[0] == "PPOT") {
			if (fields.size() < 4) {
				sprintf(err_msg, "Insufficient data");
				flag_error (err_msg, __FILE__, __LINE__);
				return FILE_ERROR;
			}
			Interaction interaction;
			
			vector <double> force_args;
			for (int i = 4; i < fields.size(); ++i) {
				force_args.push_back(atof(fields[i].c_str()));
			}
			
			force_energy_ptr fn = get_fn(fields[3], &force_args, &rcut_max);
			if (fn == NULL) {
				sprintf(err_msg, "Invalid function.");
				flag_error (err_msg, __FILE__, __LINE__);
				return ILLEGAL_VALUE;
			}
			
			atom_pairs.push_back(make_pair(fields[1], fields[2]));
			interaction.set_force_energy(fn);
			interaction.set_args(force_args);
			
			inters_PPOT.push_back(interaction);
		}
		else if (fields[0] == "BOND") {
			if (fields.size() < 3) {
				sprintf(err_msg, "Insufficient data");
				flag_error (err_msg, __FILE__, __LINE__);
				return FILE_ERROR;
			}
			Interaction interaction;
			
			// get and check force arguments
			vector <double> force_args;
			for (int i = 3; i < fields.size(); ++i) {
				force_args.push_back(atof(fields[i].c_str()));
			}
			
			force_energy_ptr fn = get_fn(fields[2], &force_args, &rcut_max);
			if (fn == NULL) {
				sprintf(err_msg, "Invalid function.");
				flag_error (err_msg, __FILE__, __LINE__);
				return ILLEGAL_VALUE;
			}
			
			bond_list.push_back(fields[1]);
			interaction.set_args(force_args);
			interaction.set_force_energy(fn);
			inters_BOND.push_back(interaction);
		} 
		else {
			sprintf(err_msg, "Unrecognized keyword %s", fields[0].c_str());
			flag_error (err_msg, __FILE__, __LINE__);
			return FILE_ERROR;
		}
	}
	fclose(input);
	
	// find total number of atoms
	int total_atoms = sys->global_atom_types.size();

	// resize interact
	try {
		sys->interact.resize(total_atoms);
	}
	catch (bad_alloc& ba) {
		sprintf(err_msg, "Out of memory");
		flag_error (err_msg, __FILE__, __LINE__);
		return FILE_ERROR;
	}
	for (int i = 0; i < total_atoms; ++i) {
		try {
			sys->interact[i].resize(total_atoms);
		}
		catch (bad_alloc& ba) {
			sprintf(err_msg, "Out of memory");
			flag_error (err_msg, __FILE__, __LINE__);
			return FILE_ERROR;
		}
	}
	
	// initialize interact with specified pair potentials
	for (int i = 0; i < total_atoms; ++i) {
		for (int j = 0; j < i; ++j) {
			pair <string, string> a_pair;
			a_pair.first = sys->global_atom_types[i];
			a_pair.second = sys->global_atom_types[j];
			
			int index = distance(atom_pairs.begin(),find(atom_pairs.begin(), atom_pairs.end(), a_pair));
			if (index == atom_pairs.size()) {
				a_pair.second = sys->global_atom_types[i];
				a_pair.first = sys->global_atom_types[j];
				index = distance(atom_pairs.begin(),find(atom_pairs.begin(), atom_pairs.end(), a_pair));
				if (index == atom_pairs.size()) {
					sprintf(err_msg, "Not a valid pair");
					flag_error (err_msg, __FILE__, __LINE__);
					return FILE_ERROR;
				}		
			}
					     
			sys->interact[i][j] = inters_PPOT[index];
			sys->interact[j][i] = inters_PPOT[index];
		}
	}

	// go through bonded atoms and change their interaction to bonding potentials in interact
	for (int i = 0; i < sys->nbonds(); ++i) {
		string b_type = sys->bond_name(sys->get_bond_type(i));
		int b_index = sys->bond_type(b_type);
		if (b_index < 0) {
			sprintf(err_msg, "Bond type %s has not been specified in %s", b_type.c_str(), filename_cstr);
			flag_error (err_msg, __FILE__, __LINE__);
			return ILLEGAL_VALUE;
		}

		int index = distance(bond_list.begin(), find(bond_list.begin(), bond_list.end(), b_type));
		if (index >= inters_BOND.size()) {
			sprintf(err_msg, "Bond type %s was found in input coordinates but not in %s", b_type.c_str(), filename_cstr);
			flag_error (err_msg, __FILE__, __LINE__);
			return FILE_ERROR;
		}
		
		// add interaction to interact
		// symmetric
		const pair <int, int> b_pair = sys->get_bond(i);
		sys->interact[b_pair.first][b_pair.second] = inters_BOND[index];
		sys->interact[b_pair.second][b_pair.first] = inters_BOND[index];
	}

	// set max cutoff radius for "skin" calculations with MPI later on
	sys->set_max_rcut(rcut_max);
	
	sprintf(err_msg, "Successfully read energy parameters from %s on rank %d", filename_cstr, rank);
	flag_notify (err_msg, __FILE__, __LINE__);
	return 0;
}

//! Given a name, return the force_energy_ptr associated with it.  Also check the arguments that will be passed to it are in acceptable range.
/*!
 \param [in] name Name of interaction, i.e. "fene" or "slj"
 \param [in] \*args Pointer to vector of arguments that will be passed to this interaction later on
 */
force_energy_ptr get_fn(const string name, vector <double> *args, double *r_cut_max) {
	char err_msg[MYERR_FLAG_SIZE];

	if (name == "fene" || name == "FENE") {
		if (args->size() != 5) {
			sprintf(err_msg, "Fene bond type %s has %d argument(s), expects 5 (epsilon, sigma, delta, k, r0)", name.c_str(), (int) args->size());
			flag_error (err_msg, __FILE__, __LINE__);
			return NULL;
		}
		char* names[] = {"epsilon", "sigma", "delta", "k", "r0"};
		for (int i = 0; i < 5; ++i) {
			if (args->at(i) < 0.0) {
				sprintf(err_msg, "Fene bond type %s argument %s = %g, must be > 0", name.c_str(), names[i], (double) args->at(i));
				flag_error (err_msg, __FILE__, __LINE__);
				return NULL;
			}
		}
		
		// for fene bonds, the equivalent "r_cut" is r0
		if (args->at(4) > *r_cut_max) {
			*r_cut_max = args->at(4);
		}
		return &fene;
	}
	else if (name == "harmonic" || name == "HARM") {
		if (args->size() != 2) {
			sprintf(err_msg, "Harmonic bond type %s has %d argument(s), expects 2 (k, r0)", name.c_str(), (int)args->size());
			flag_error (err_msg, __FILE__, __LINE__);
			return NULL;
		}
		char* names[] = {"k", "r0"};
		for (int i = 0; i < 2; ++i) {
			if (args->at(i) < 0.0) {
				sprintf(err_msg, "Harmonic bond type %s argument %s = %g, must be > 0", name.c_str(), names[i], (double) args->at(i));
				flag_error (err_msg, __FILE__, __LINE__);
				return NULL;
			}
		}
		return &harmonic;
	}
	else if (name == "slj" || name == "SLJ") {
		if (args->size() != 5) {
			sprintf(err_msg, "Shifted Lennard-Jones interaction type %s has %d argument(s), expects 5 (epsilon, sigma, delta, U_shift, r_cut)", name.c_str(), (int)args->size());
			flag_error (err_msg, __FILE__, __LINE__);
			return NULL;
		}
		char* names[] = {"epsilon", "sigma", "delta", "U_shift", "r_cut"};
		for (int i = 0; i < 5; ++i) {
			// skip U_shift test which can be < 0 if we want
			if (i == 3) {
				continue;
			}
			if (args->at(i) < 0.0) {
				sprintf(err_msg, "Shifted Lennard-Jones interaction type %s argument %s = %g, must be > 0", name.c_str(), names[i], (double) args->at(i));
				flag_error (err_msg, __FILE__, __LINE__);
				return NULL;
			}
		}
		
		// test r_cut_max
		if (*r_cut_max < args->at(4)) {
			*r_cut_max = args->at(4);
		}
		
		// square to rcut value here since it is stored internally this way; it will be much faster to have this as such
		args->at(4) = args->at(4)*args->at(4);
		return &slj;
	}
	else {
		sprintf(err_msg, "Undefined interaction type %s", name.c_str());
		flag_error (err_msg, __FILE__, __LINE__);
		return NULL;
	}	
}
