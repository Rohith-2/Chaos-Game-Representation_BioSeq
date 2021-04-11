from collections import defaultdict
import matplotlib.pyplot as plt

from Bio import SeqIO

try:
	from pycgr import helper
except ModuleNotFoundError:
	import helper
# defining cgr graph
# CGR_CENTER = (0.5, 0.5)
CGR_X_MAX = 1
CGR_Y_MAX = 1
CGR_X_MIN = 0
CGR_Y_MIN = 0
CGR_A = (CGR_X_MIN, CGR_Y_MIN)
CGR_T = (CGR_X_MAX, CGR_Y_MIN)
CGR_G = (CGR_X_MAX, CGR_Y_MAX)
CGR_C = (CGR_X_MIN, CGR_Y_MAX)
CGR_CENTER = ((CGR_X_MAX - CGR_Y_MIN) / 2, (CGR_Y_MAX - CGR_Y_MIN) / 2)

# Add color code for each element


def empty_dict():
	"""
	None type return vessel for defaultdict
	:return:
	"""
	return None


CGR_DICT = defaultdict(
	empty_dict,
	[
		('A', CGR_A),  # Adenine
		('T', CGR_T),  # Thymine
		('G', CGR_G),  # Guanine
		('C', CGR_C),  # Cytosine
		('U', CGR_T),  # Uracil demethylated form of thymine
		('a', CGR_A),  # Adenine
		('t', CGR_T),  # Thymine
		('g', CGR_G),  # Guanine
		('c', CGR_C),  # Cytosine
		('u', CGR_T)  # Uracil/Thymine
		]
)


def fasta_reader(fasta):
	"""Return a generator with sequence description and sequence
	:param fasta: str filename
	"""
	# TODO: modify it to be capable of reading genebank etc
	flist = SeqIO.parse(fasta, "fasta")
	for i in flist:
		yield i.description, i.seq


def mk_cgr(seq):
	"""Generate cgr
	:param seq: list of nucleotide
	:return cgr: [['nt', (x, y)]] List[List[Tuple(float, float)]]
	"""
	cgr = []
	cgr_marker = CGR_CENTER[:
		]    # The center of square which serves as first marker
	for s in seq:
		cgr_corner = CGR_DICT[s]
		if cgr_corner:
			cgr_marker = (
				(cgr_corner[0] + cgr_marker[0]) / 2,
				(cgr_corner[1] + cgr_marker[1]) / 2
			)
			cgr.append([s, cgr_marker])
		else:
			sys.stderr.write("Bad Nucleotide: " + s + " \n")

	return cgr


def mk_plot(cgr, name, figid):
	"""Plotting the cgr
		:param cgr: [(A, (0.1, 0.1))]
		:param name: str
		:param figid: int
		:return dict: {'fignum': figid, 'title': name, 'fname': helper.slugify(name)}
	"""
	x_axis = [i[1][0] for i in cgr]
	y_axis = [i[1][1] for i in cgr]
	plt.figure(figid)
	plt.title("Chaos Game Representation\n" + name, wrap=True)
	# diagonal and vertical cross
	# plt.plot([x1, x2], [y1, y2])
	# plt.plot([0.5,0.5], [0,1], 'k-')
	plt.plot([CGR_CENTER[0], CGR_CENTER[0]], [0, CGR_Y_MAX], 'k-')

	# plt.plot([0,1], [0.5,0.5], 'k-')
	plt.plot([CGR_Y_MIN, CGR_X_MAX], [CGR_CENTER[1], CGR_CENTER[1]], 'k-')
	plt.scatter(x_axis, y_axis, alpha=0.5, marker='.')

	return {'fignum': figid, 'title': name, 'fname': helper.slugify(name)}


def write_figure(fig, output_dir, dpi=300):
	"""Write plot to png
	:param fig:  {'fignum':figid, 'title':name, 'fname':helper.slugify(name)}
	:param dpi: int dpi of output
	:param output_dir: str
	Usage:
		figures = [mk_plot(cgr) for cgr in all_cgr]
		for fig in figures:
			write_figure(fig, "/var/tmp/")
		The figid in the mk_plot's return dict must be present in plt.get_fignums()
	"""
	all_figid = plt.get_fignums()
	if fig['fignum'] not in all_figid:
		raise ValueError("Figure %i not present in figlist" % fig['fignum'])
	plt.figure(fig['fignum'])
	target_name = os.path.join(
		output_dir,
		helper.slugify(fig['fname']) + ".png"
	)
	plt.savefig(target_name, dpi=dpi)


def get_args():
	"""Get args"""
	import argparse
	parser = argparse.ArgumentParser(
		description="Chaos Game Representation",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)

	parser.add_argument(
		"--dest-dir",
		"-d",
		dest="save_dir",
		action='store',
		default="/var/tmp",
		help="Where do you want to save the resulting CGRs?"
	)
	parser.add_argument(
		"--show",
		"-p",
		action='store_true',
		default=False,
		help="Do you want to display resulting CGRs?"
	)
	parser.add_argument(
		"--save",
		"-s",
		action='store_true',
		default=False,
		help="Do you want to save resulting CGRs?"
	)
	parser.add_argument(
		"--dpi",
		"-i",
		dest="dpi",
		action='store',
		type=int,
		default=300,
		help="dpi for generated image"
	)
	parser.add_argument('files', nargs='*')
	args = parser.parse_args()
	if not os.path.isdir(args.save_dir):
		raise FileNotFoundError(
			'The requested save_dir %s does not exist!!!' % args.save_dir
		)
	if not args.files:
		raise TypeError('Not Enough Argument; expected at least 1')
	for i in args.files:
		if not os.path.isfile(i):
			raise TypeError(
				'File %s does not exist or is not the right type' % i
			)
	return args


if __name__ == '__main__':
	fig_id = 1
	args = get_args()
	my_plots = []    # new plot
	mycgr = []    # hold new cgrs; pointless for now
	for i in args.files:
		fasta_seq = fasta_reader(i)
		for name, seq in fasta_seq:
			cgr = mk_cgr(seq)
			# TODO:Add facility to write cgr to file
			my_plots.append(mk_plot(cgr, name, fig_id))
			fig_id += 1
	if args.save:
		for i in my_plots:
			write_figure(i, args.save_dir, dpi=args.dpi)

	if args.show:
		plt.show()