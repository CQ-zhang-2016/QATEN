import configargparse

def buildParser():
	parser = configargparse.ArgParser(default_config_files=['settings.conf'])

	# Data source
	parser.add('--data_dir',    default='./data/',              	help='Source directory for pdb files')
	parser.add('--save_dir',    default='./result/',      			help='Destination directory for results')
	parser.add('--atom_init',   default='protein_atom_init.json',   help='atom_init filename')
	parser.add('--pretrained',  default='QATEN.pth.tar',			help='Path to pretrained model')

	# Training setup
	parser.add('--batch_size',  default=1,                          help='Batch size for training',    type=int)
	parser.add('--workers',     default=0,                 			help='Number of workers for data loading',  type=int)

	return parser
