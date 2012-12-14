#include "read_interaction.h"

int read_interactions(const string filename, System *sys) {
	const int buffsize = 1000;
	char buff[buffsize];
	
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
			
			atom_pairs.push_back(make_pair(fields[1], fields[2]));
			force_energy_ptr fn = get_fn(fields[3]);
			if (fn == NULL) {
				sprintf(err_msg, "Invalid function.");
				flag_error (err_msg, __FILE__, __LINE__);
				return FILE_ERROR;
			}
			interaction.set_force_energy(fn);
			
			vector <double> force_args;
			for (int i = 4; i < fields.size(); ++i) {
				force_args.push_back(atof(fields[i].c_str()));
			}
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
			
			bond_list.push_back(fields[1]);
			force_energy_ptr fn = get_fn(fields[2]);
			if (fn == NULL) {
				sprintf(err_msg, "Invalid function.");
				flag_error (err_msg, __FILE__, __LINE__);
				return FILE_ERROR;
			}
			interaction.set_force_energy(fn);
			
			vector <double> force_args;
			for (int i = 3; i < fields.size(); ++i) {
				force_args.push_back(atof(fields[i].c_str()));
			}
			interaction.set_args(force_args);
			
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

	// resize interact_
	try {
		sys->interact_.resize(total_atoms);
	}
	catch (bad_alloc& ba) {
		sprintf(err_msg, "ZOMG no more memory");
		flag_error (err_msg, __FILE__, __LINE__);
		return FILE_ERROR;
	}
	for (int i = 0; i < total_atoms; ++i) {
		try {
			sys->interact_[i].resize(total_atoms);
		}
		catch (bad_alloc& ba) {
			sprintf(err_msg, "ZOMG no more memory");
			flag_error (err_msg, __FILE__, __LINE__);
			return FILE_ERROR;
		}
	}
	
	// initialize interact_ with default pair potentials
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
					     
			sys->interact_[i][j] = inters_PPOT[index];
			sys->interact_[j][i] = inters_PPOT[index];
		}
	}
	
	// go thorugh bonded atoms and change their interaction in interact_
	for (int i = 0; i < sys->nbonds(); ++i) {
		
		string b_type = sys->bond_name(i);
		int index = distance(bond_list.begin(), find(bond_list.begin(), bond_list.end(), b_type));
		
		if (index == atom_pairs.size()) {
			sprintf(err_msg, "Not a valid pair");
			flag_error (err_msg, __FILE__, __LINE__);
			return FILE_ERROR;
		}		
		
		
		// add interaction to interact_
		// symmetric
		const pair <int, int> b_pair = sys->get_bond(i);
		sys->interact_[b_pair.first][b_pair.second] = inters_BOND[index];
		sys->interact_[b_pair.second][b_pair.first] = inters_BOND[index];
	}
	return 0;
}

force_energy_ptr get_fn(const string name) {
	char err_msg[MYERR_FLAG_SIZE];

	if (name == "fene" || name == "FENE") {
		return &fene;
	}
	else if (name == "harmonic" || name == "HARM") {
		return &harmonic;
	}
	else if (name == "slj" || name == "SLJ") {
		return &slj;
	}
	else {
		sprintf(err_msg, "Undefined interaction type");
		flag_error (err_msg, __FILE__, __LINE__);
		return NULL;
	}	
}
