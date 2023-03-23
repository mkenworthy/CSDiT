import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import unique,vstack,Table
import paths

from photomutils import *

obj='ASASSN-21sa'

fin=f'obs_{obj}_ASASSN.ecsv'
t = ascii.read(paths.data / fin)


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

ty = dict(color='k', fontsize=16, fontweight='bold', va='bottom')
ax.text(0.05, 0.05, obj, transform=ax.transAxes, **ty)

ax.legend()
#plt.draw()
fig.savefig(paths.figures / f'{obj}.pdf', bbox_inches='tight')
