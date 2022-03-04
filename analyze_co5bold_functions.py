# functions and tools for loading and analyzing output
# from co5bold.




"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.io.idl import readsav

a=readsav('st29gm04n002_seq.sav')
seq = a['seq']

# print(a.keys())  # get restored variables

# loaded sav file is as a dict file, with each restored variable as a numpy.recarray file
# (i.e. an ndarray that allows field access using attributes) so get the tag names
# (tag_names in IDL) via the dtype.names
tag_names = a['seq'].dtype.names
# print(tag_names)

# then can access tag data via
# B = a['seq']['tag_name'][0]

# misc
# can always check type(XX) and if array, check dims via .shape

# e.g. plot the radial velocity
B = a['seq']['v3'][0]
B = (1e-5 * seq.rhov3 / seq.rho)[0]
xx = a['seq']['time'][0] / (365*24*3600)
xx = xx - xx[0]
yy = seq.xc3_radial[0] / 6.98000e+10
hh = B

plt.pcolormesh(xx, yy, hh.transpose(), cmap=plt.cm.get_cmap('bwr'), shading='auto')
cbar=plt.colorbar()
cbar.set_label('Radial velocity ' + r'$[kms^{-1}]$', rotation=270, labelpad=25)
plt.xlabel(r'Time $[yr]$')
plt.ylabel(r'$R \: \: [R_{\odot}]$')
plt.tick_params(direction='in', top ='on', right = 'on')
plt.show()
"""