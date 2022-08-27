from __future__ import print_function, division
import os, functools, math, csv, random, pickle, json
import numpy as np

import torch
from torch.utils.data import Dataset, DataLoader
from torch.utils.data.dataloader import default_collate
from torch.utils.data.sampler import SubsetRandomSampler, SequentialSampler
from collections import defaultdict as ddict
import json

def get_loader(test_dataset, collate_fn=default_collate, batch_size=64, num_workers=1, pin_memory=False):


	test_indices    = [i for i in range(len(test_dataset))]


	random.shuffle(test_indices)

	test_sampler    = SubsetRandomSampler(test_indices)

	test_loader     = DataLoader(test_dataset, batch_size=batch_size, sampler=test_sampler, num_workers=num_workers, collate_fn=collate_fn, pin_memory=pin_memory)
	return test_loader



def collate_pool(dataset_list):
	N   = max([x[0][0].size(0) for x in dataset_list])  # max atoms
	M   = dataset_list[0][0][1].size(1)                 # num neighbors are same for all so take the first value
	B   = len(dataset_list)                             # Batch size
	h_b = dataset_list[0][0][1].size(2)                 # Edge feature length

	final_protein_atom_fea = torch.zeros(B, N)
	final_nbr_fea          = torch.zeros(B, N, M, h_b)
	final_nbr_fea_idx      = torch.zeros(B, N, M, dtype=torch.long)
	final_atom_amino_idx   = torch.zeros(B, N)
	final_atom_mask        = torch.zeros(B, N)
	amino_base_idx         = 0

	batch_protein_ids, amino_crystal = [], 0
	for i, ((protein_atom_fea, nbr_fea, nbr_fea_idx, atom_amino_idx), protein_id) in enumerate(dataset_list):
		num_nodes                             = protein_atom_fea.size(0)
		final_protein_atom_fea[i][:num_nodes] = protein_atom_fea.squeeze()
		final_nbr_fea[i][:num_nodes]          = nbr_fea
		final_nbr_fea_idx[i][:num_nodes]      = nbr_fea_idx
		final_atom_amino_idx[i][:num_nodes]   = atom_amino_idx + amino_base_idx
		final_atom_amino_idx[i][num_nodes:]   = amino_base_idx
		amino_base_idx                       += torch.max(atom_amino_idx) + 1
		final_atom_mask[i][:num_nodes]        = 1


		batch_protein_ids.append(protein_id)
		amino_crystal += 1

	
	return (final_protein_atom_fea, final_nbr_fea, final_nbr_fea_idx, final_atom_amino_idx, final_atom_mask),\
			(batch_protein_ids)


class GaussianDistance(object):
	"""
	Expands the distance by Gaussian basis.

	Unit: angstrom
	"""
	def __init__(self, dmin, dmax, step, var=None):
		"""
		Parameters
		----------
		dmin: float
			Minimum interatomic distance
		dmax: float
			Maximum interatomic distance
		step: float
			Step size for the Gaussian filter
		"""
		assert dmin < dmax
		assert dmax - dmin > step
		self.filter = np.arange(dmin, dmax+step, step)
		if var is None:
			var = step
		self.var = var

	def expand(self, distances):
		"""
		Apply Gaussian distance filter to a numpy distance array

		Parameters
		----------
		distance: np.array shape n-d array
			A distance matrix of any shape

		Returns
		-------
		expanded_distance: shape (n+1)-d array
			Expanded distance matrix with the last dimension of length
			len(self.filter)
		"""
		return np.exp(-(distances[..., np.newaxis] - self.filter)**2 / self.var**2)


class AtomInitializer(object):
	"""
	Base class for intializing the vector representation for atoms.
	"""
	def __init__(self, atom_types):
		self.atom_types = set(atom_types)
		self._embedding = {}

	def get_atom_fea(self, atom_type):
		assert atom_type in self.atom_types
		return self._embedding[atom_type]

	def load_state_dict(self, state_dict):
		self._embedding = state_dict
		self.atom_types = set(self._embedding.keys())
		self._decodedict = {idx: atom_type for atom_type, idx in self._embedding.items()}

	def state_dict(self):
		return self._embedding

	def decode(self, idx):
		if not hasattr(self, '_decodedict'):
			self._decodedict = {idx: atom_type for atom_type, idx in self._embedding.items()}
		return self._decodedict[idx]


class AtomCustomJSONInitializer(AtomInitializer):
	"""
	Initialize atom feature vectors using a JSON file, which is a python
	dictionary mapping from element number to a list representing the
	feature vector of the element.

	Parameters
	----------
	elem_embedding_file: str
		The path to the .json file
	"""
	def __init__(self, elem_embedding_file):
		with open(elem_embedding_file) as f:
			elem_embedding = json.load(f)
		elem_embedding = {key: value for key, value in elem_embedding.items()}
		atom_types = set(elem_embedding.keys())
		super(AtomCustomJSONInitializer, self).__init__(atom_types)
		counter = 0
		for key, _ in elem_embedding.items():
			self._embedding[key] = counter; counter += 1


class ProteinDataset(Dataset):

	def __init__(self, pkl_dir, atom_init_filename, random_seed=123):
		assert os.path.exists(pkl_dir), '{} does not exist!'.format(pkl_dir)
		self.pkl_dir = pkl_dir
		self.pkl_data = []
		for i in os.listdir(pkl_dir):
			if '.pkl' in i:
				self.pkl_data.append(i)


		random.seed(random_seed)

		assert os.path.exists(atom_init_filename), '{} does not exist!'.format(atom_init_filename)
		self.ari                = AtomCustomJSONInitializer(atom_init_filename)
		self.gdf                = GaussianDistance(dmin=0, dmax=15, step=0.4)

	def __len__(self):
		return len(self.pkl_data)

	def __getitem__(self, idx):
		return self.get_idx(idx)

	def get_idx(self, idx):
		prot_name = self.pkl_data[idx]


		with open(self.pkl_dir + prot_name, 'rb') as f:
			protein_atom_fea    = torch.Tensor(np.vstack([self.ari.get_atom_fea(atom) for atom in pickle.load(f)]))     # Atom features (here one-hot encoding is used)
			nbr_info            = pickle.load(f)                                                                        # Edge features for each atom in the graph
			nbr_fea_idx         = torch.LongTensor(pickle.load(f))                                                      # Edge connections that define the graph

			atom_amino_idx  = torch.LongTensor(pickle.load(f))  # Mapping that denotes which atom corresponds to which amino residue in the protein graph

			protein_id          = pickle.load(f)
			nbr_fea             = torch.Tensor(np.concatenate([self.gdf.expand(nbr_info[:,:,0]), nbr_info[:,:,1:]], axis=2))    # Use Gaussian expansion for edge distance
			

		return (protein_atom_fea, nbr_fea, nbr_fea_idx, atom_amino_idx), protein_id
	
