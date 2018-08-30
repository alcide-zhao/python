# -*- coding: utf-8 -*-
import matplotlib.pyplot as pltt
import os
import subprocess
import tempfile
import sys


reload(sys)
sys.setdefaultencoding('utf8')


"""
The following functions can be used by scripts to get the sizes of
the various elements of the figures.
"""
def plot_setup():
	params = {
		'figure.dpi': 600,
		'text.usetex': False,
		'axes.labelsize': 10,
		'axes.titlesize': 9,
		'axes.linewidth': 0.6,
		'text.fontsize': 8,
		'axes.annotate.size':8,
		'legend.fontsize': 8,
		'legend.markerscale': 1.0,
		'legend.handlelength': 3.0,
		'legend.columnspacing': 1.0,
		'legend.markerfirst':True,
		'xtick.labelsize': 8,
		'ytick.labelsize': 8,
		'ps.useafm': True,
		# 'pdf.fonttype': 42
		'ps.fonttype' : 42,
		'font.family' : 'sans-serif',
		'font.size': 12,
		'font.sans-serif': 'Helvetica'
		}
	pltt.rcParams.update(params)


def save_fig(fig, file_name, fmt, dpi=600, tight=True):
	"""
	Save a Matplotlib figure as EPS/PNG/PDF to the given path and trim it.
	"""
	if not fmt:
		fmt = file_name.strip().split('.')[-1]
	if fmt not in ['eps', 'png', 'pdf']:
		raise ValueError('unsupported format: %s' % (fmt,))
	extension = '.%s' % (fmt,);print extension
	if not file_name.endswith(extension):
		file_name += extension
	file_name = os.path.abspath(file_name);
    # with tempfile.NamedTemporaryFile().encode('utf-8') as tmp_file:
        # tmp_name = tmp_file.name + extension; print tmp_name
	tmp_name = 'cache'+ extension;
    # save figure
	if tight:
		fig.savefig(tmp_name, dpi=dpi, bbox_inches='tight')
	else:
		fig.savefig(tmp_name, dpi=dpi)

    # trim it
	if fmt == 'eps':
		subprocess.call('epstool --bbox --copy %s %s' %
                        (tmp_name, file_name), shell=True)
	elif fmt == 'png':
		subprocess.call('convert %s -trim %s' %
                        (tmp_name, file_name), shell=True)
	elif fmt == 'pdf':
		subprocess.call('pdfcrop %s %s' % (tmp_name, file_name), shell=True)