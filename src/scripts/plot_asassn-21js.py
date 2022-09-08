import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import unique,vstack,Table
import paths

from photomutils import *

fin='obs_ASASSN-21js_ASASSN.ecsv'
t = ascii.read(paths.data / fin)

obj='ASASSN-21js'

# get a list of the unique bandpasses
t_by_filter = t.group_by('Filter')
print('all observed photometric bands:')
print(t_by_filter.groups.keys)

fig, (ax) = plt.subplots(1,1,figsize=(12,6))
ax.set_ylabel('Normalised flux')
ax.set_xlabel('Epoch [MJD]')
ax.set_title('data from {}'.format(fin))

for key, group in zip(t_by_filter.groups.keys, t_by_filter.groups):
    ax.errorbar(group['MJD'],group['fnorm'],yerr=group['fnormerr'],label=key['Filter'],fmt='.',alpha=0.8)

ax.legend()
plt.show()
fig.savefig(paths.figures / 'asassn-21js.pdf')
