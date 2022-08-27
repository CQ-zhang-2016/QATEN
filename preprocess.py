import pickle, json
import numpy as np
from subprocess import call
from collections import defaultdict as ddict
import os

def createSortedNeighbors(contacts, bonds, max_neighbors=50):
	"""
	Given a list of list contacts of the form [index1, index2, distance, x1, y1, z1, x2, y2, z2]
	generate the k nearest neighbors for each atom based on distance.
	Parameters
	----------
	contacts        : list
	bonds           : list
	max_neighbors   : int
		Limit for the maximum neighbors to be set for each atom.
	"""
	bond_true       = 1
	bond_false      = 0
	neighbor_map    = ddict(list)
	dtype           = [('index2', int), ('distance', float), ('x1', float), ('y1', float), ('z1', float), ('bool_bond', int)]
	idx             = 0

	for contact in contacts:
		if ([contact[0], contact[1]] or [contact[1], contact[0]]) in bonds:
			neighbor_map[contact[0]].append((contact[1], contact[2], contact[3], contact[4], contact[5], bond_true))
			neighbor_map[contact[1]].append((contact[0], contact[2], contact[6], contact[7], contact[8], bond_true))
		else:
			neighbor_map[contact[0]].append((contact[1], contact[2], contact[3], contact[4], contact[5], bond_false))
			neighbor_map[contact[1]].append((contact[0], contact[2], contact[6], contact[7], contact[8], bond_false))
		idx += 1

	for k, v in neighbor_map.items():
		if len(v) < max_neighbors:
			true_nbrs = np.sort(np.array(v, dtype=dtype), order='distance', kind='mergesort').tolist()[0:len(v)]
			true_nbrs.extend([(0, 0, 0, 0, 0, 0) for _ in range(max_neighbors - len(v))])
			neighbor_map[k] = true_nbrs
		else:
			neighbor_map[k] = np.sort(np.array(v, dtype=dtype), order='distance', kind='mergesort').tolist()[0:max_neighbors]

	return neighbor_map


def pdb2pkl(datapath):
	print('start getting features')
	
	for directory in os.listdir(datapath):
		if '.pdb' not in directory: continue
		json_filepath   		= datapath + directory.strip('.pdb') + '.json'
		command         		= './preprocess/get_features -i ' + datapath + directory + ' -j ' + json_filepath
		call(command, shell=True)
		print('Processing ', directory)

		with open(json_filepath, 'r') as file:
			json_data = json.load(file)

		neighbor_map    = createSortedNeighbors(json_data['contacts'], json_data['bonds'])
		amino_atom_idx  = json_data['res_idx']
		atom_fea        = json_data['atoms']
		nbr_fea_idx     = np.array([list(map(lambda x: x[0], neighbor_map[idx])) for idx in range(len(json_data['atoms']))])
		nbr_fea         = np.array([list(map(lambda x: x[1:], neighbor_map[idx])) for idx in range(len(json_data['atoms']))])

		with open(datapath + directory.strip('.pdb')+'.pkl', 'wb') as file:
			pickle.dump(atom_fea, file)
			pickle.dump(nbr_fea, file)
			pickle.dump(nbr_fea_idx, file)
			pickle.dump(amino_atom_idx, file)
			pickle.dump(directory.strip('.pdb'), file)

	print('over!')
	
