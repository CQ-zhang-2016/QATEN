import sys, time, csv, os, random, math, argparse
import numpy as np
import torch
from collections import OrderedDict
from arguments import buildParser
import json
from model import QATEN
from data import ProteinDataset, collate_pool, get_loader
from utils import Normalizer, count_parameters, randomSeed
import config as cfg
import time
from preprocess import pdb2pkl

def main():
	global args

	parser  = buildParser()
	args    = parser.parse_args()
	
	try:os.mkdir(args.save_dir)
	except:pass
	
	pdb2pkl(args.data_dir)

	print('Torch Device being used: ', cfg.device)
	test_dataset = ProteinDataset(args.data_dir, args.atom_init)
	print('Dataset length: ', len(test_dataset))

	# load all model args from pretrained model
	model_checkpoint    = torch.load(args.pretrained, map_location=lambda storage, loc: storage)
	model_args          = argparse.Namespace(**model_checkpoint['args'])
	print(model_args)
	print("=> loaded model params '{}'".format(args.pretrained))

	# build model
	kwargs = {
		'atom_init'     : model_args.atom_init,       # Atom Init filename
		'h_a'           : model_args.h_a,             # Dim of the hidden atom embedding learnt
		'h_g'           : model_args.h_g,             # Dim of the hidden graph embedding after pooling
		'head_num'		: model_args.head_num,
		'n_conv'        : model_args.n_conv,          # Number of GCN layers

		'random_seed'   : model_args.seed,            # Seed to fix the simulation
		'lr'            : model_args.lr,              # Learning rate for optimizer
	}

	structures, _    = test_dataset[0]
	h_b                 = structures[1].shape[-1]
	kwargs['h_b']       = h_b                       # Dim of the bond embedding initialization


	model = QATEN(**kwargs)
	model.cuda()

	print('Trainable Model Parameters: ', count_parameters(model))
	# Create dataloader to iterate through the dataset in batches
	test_loader = get_loader(test_dataset, 
								collate_fn    = collate_pool,
								num_workers   = args.workers,
								batch_size    = args.batch_size,
								pin_memory    = False)

	

	normalizer_global    = Normalizer(torch.tensor([0.0]))
	normalizer_local    = Normalizer(torch.tensor([0.0]))

	model.load_state_dict(model_checkpoint['state_dict'])
	normalizer_global.load_state_dict(model_checkpoint['normalizer_global'])
	normalizer_local.load_state_dict(model_checkpoint['normalizer_local'])
	[_, _] = trainModel(test_loader, model, normalizer_global, normalizer_local, args.save_dir)




def trainModel(data_loader, model, normalizer_global, normalizer_local, save_dir):


	for protein_batch_iter, (input_data, batch_data) in enumerate(data_loader):

		input_var = getInputs(input_data)
		# evaluate one iteration
		with torch.no_grad():
			# Switch to evaluation mode
			model.eval()
			predicted = model(input_var)

		test_pred_global    = normalizer_global.denorm(predicted[0].data).squeeze().tolist()
		test_pred_local    = normalizer_local.denorm(predicted[1].data).squeeze().tolist()
		res = '{"global": ' + str(test_pred_global)+', \n' + '"local": ' + str(test_pred_local)+'}'

		with open(save_dir + '/' + batch_data[0] + '_QATEN_result.json','w')as f:
			f.write(res)

	return None, None


def getInputs(inputs):

	input_var               = [inputs[0].cuda(), inputs[1].cuda(), inputs[2].cuda(), inputs[3].cuda(), inputs[4].cuda()]
	return input_var


if __name__ == '__main__':
	start = time.time()
	main()
	print('Time taken: ', time.time() - start)
