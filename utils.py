import torch
import config as cfg

class Normalizer(object):
	"""Normalize a Tensor and restore it later. """
	def __init__(self, tensor):
		"""tensor is taken as a sample to calculate the mean and std"""
		self.mean = torch.mean(tensor).type(cfg.FloatTensor)
		self.std = torch.std(tensor).type(cfg.FloatTensor)

	def norm(self, tensor):
		if self.mean != self.mean or self.std != self.std:
			return tensor
		return (tensor - self.mean) / self.std

	def denorm(self, normed_tensor):
		if self.mean != self.mean or self.std != self.std:
			return normed_tensor
		return normed_tensor * self.std + self.mean

	def state_dict(self):
		return {'mean': self.mean,
				'std': self.std}

	def load_state_dict(self, state_dict):
		self.mean = state_dict['mean']
		self.std = state_dict['std']


def count_parameters(model):
	return sum(p.numel() for p in model.parameters() if p.requires_grad)


def randomSeed(random_seed):
	"""Given a random seed, this will help reproduce results across runs"""
	if random_seed is not None:
		torch.manual_seed(random_seed)

