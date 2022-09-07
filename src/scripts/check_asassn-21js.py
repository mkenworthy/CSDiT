import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import unique,vstack
import paths

from photomutils import *

fin='ASASSN-21js/light_curve_f9818a9a-2dfc-4c33-ad97-826ae7b78a33.csv'
t = ascii.read(paths.data / fin)

#      HJD           UT Date       Camera FWHM Limit   mag   mag_err flux(mJy) flux_err Filter
# ------------- ------------------ ------ ---- ------ ------ ------- --------- -------- ------
# 2457420.65322 2016-02-02.1500246     be 1.46 17.458  13.45   0.005    15.995     0.08      V

t['MJD'] = t['HJD']-2400000.5

# reject noisy points
t = t[(t['flux(mJy)']<30)]
print('rejected ASASSN points with high fluxes in both bands')

fig, (ax) = plt.subplots(1,1,figsize=(12,6))
ax.errorbar(t['MJD'],t['flux(mJy)'],yerr=t['flux_err'],fmt='.')
ax.set_ylabel('Flux [mJy]')
ax.set_xlabel('Epoch [MJD]')
ax.set_title('data from {}'.format(fin))
fig.savefig('_check_asassn0.pdf')


# get a list of the unique bandpasses
t_by_filter = t.group_by('Filter')
print('all observed photometric bands:')
print(t_by_filter.groups.keys)

fig, (ax) = plt.subplots(1,1,figsize=(12,6))
ax.set_ylabel('Flux [mJy]')
ax.set_xlabel('Epoch [MJD]')
ax.set_title('data from {}'.format(fin))


for key, group in zip(t_by_filter.groups.keys, t_by_filter.groups):
    ax.errorbar(group['MJD'],group['flux(mJy)'],yerr=group['flux_err'],label=key['Filter'],fmt='.')

    print('')

ax.legend()
fig.savefig('_check_asassn1.pdf')

(tV, tg) = t_by_filter.groups

(cum_reqV, xV, cumV) = hist_errors (tV['flux_err'])
(cum_reqg, xg, cumg) = hist_errors (tg['flux_err'])
print(cum_reqV,cum_reqg)

fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,6))
# cumulative plots....
ax1.plot(cumV,xV)
ax2.plot(cumg,xg)
ax1.set_xlabel('V Flux err')
ax2.set_xlabel('g Flux err')
ax1.set_ylabel('Cumulative fraction')
ax2.set_ylabel('Cumulative fraction')
ax1.set_title('data from {}'.format(fin))


# check errors in both filters
ax3.hist(tV['flux_err'],bins=50,range=(0,2*cum_reqV))
ax4.hist(tg['flux_err'],bins=50,range=(0,2*cum_reqg))
ymin,ymax=ax3.get_ylim()
ax3.vlines(cum_reqV,ymin=ymin,ymax=ymax,linestyle='dashed')
ymin,ymax=ax4.get_ylim()
ax4.vlines(cum_reqg,ymin=ymin,ymax=ymax,linestyle='dashed')

ax3.set_ylabel('N')
ax3.set_xlabel('V Flux err ')
ax4.set_xlabel('g Flux err ')
fig.savefig('_check_asassn_2err.pdf')

# reject low flux points
tV = tV[(tV['flux_err']<cum_reqV)]
tg = tg[(tg['flux_err']<cum_reqg)]
print('rejecting noisy points in g and V data separately in ASASSN')

# tmin, tmax = 58290, 58360
#
# (V_flux_norm, V_flux_norm_err)= \
#     mean_rms_region(tV['MJD'],tV['flux(mJy)'],tV['flux_err'], tmin, tmax)
#
# (g_flux_norm, g_flux_norm_err)= \
#     mean_rms_region(tg['MJD'],tg['flux(mJy)'],tg['flux_err'], tmin, tmax)
#
# print(V_flux_norm, g_flux_norm)
#
# tV['fnorm'] = tV['flux(mJy)']/V_flux_norm
# tg['fnorm'] = tg['flux(mJy)']/g_flux_norm
#
# tV['fnormerr'] = tV['flux_err']/V_flux_norm
# tg['fnormerr'] = tg['flux_err']/g_flux_norm

fig, (ax) = plt.subplots(1,1,figsize=(12,6))

# 1.973 2.200 at delta 0.227 over 202 days
tedges = np.arange(57420.20-0.5,58336.02+0.5,1.0002237)

(tc_binV, fc_binV, fsigc_binV, fc_npoiV) = \
    rebinsig(tV['MJD'],tV['flux(mJy)'], tedges, ax=ax, err=tV['flux_err'])

tedges = np.arange(58282.163-0.5,59822.02+0.5,1.0002237)

(tc_bing, fc_bing, fsigc_bing, fc_npoig) = \
    rebinsig(tg['MJD'],tg['flux(mJy)'], tedges, ax=ax, err=tg['flux_err'])

ax.set_ylabel('Flux [mJy]')
ax.set_xlabel('Epoch [MJD]')
ax.set_title('data from {}'.format(fin))
ax.errorbar(tV['MJD'],tV['flux(mJy)'],yerr=tV['flux_err'],label='V',fmt='.')
ax.errorbar(tg['MJD'],tg['flux(mJy)'],yerr=tg['flux_err'],label='g',fmt='.')

ax.errorbar(tc_binV,fc_binV,yerr=fsigc_binV,label='V2',fmt='.')

ax.errorbar(tc_bing,fc_bing,yerr=fsigc_bing,label='g2',fmt='.')

fig.savefig('_check_asassn3.pdf')

fig, (ax) = plt.subplots(1,1,figsize=(12,6))

ax.set_ylabel('Flux [mJy]')
ax.set_xlabel('Epoch [MJD]')
ax.set_title('data from {}'.format(fin))

ax.errorbar(tc_binV,fc_binV,yerr=fsigc_binV/np.sqrt(fc_npoiV),label='V2',fmt='.')
ax.errorbar(tc_bing,fc_bing,yerr=fsigc_bing/np.sqrt(fc_npoig),label='g2',fmt='.')
fig.savefig('_check_asassn4.pdf')

#
# fig, (ax) = plt.subplots(1,1,figsize=(12,6))
# ax.set_ylabel('Flux [normalised]')
# ax.set_xlabel('Epoch [MJD]')
# ax.set_title('Normalised ASASSN data {}'.format(fin))
#
# ax.scatter(tV['MJD'],tV['fnorm'],label='V',
#     marker='o',edgecolors='none',alpha=0.1)
# ax.scatter(tg['MJD'],tg['fnorm'],label='g',
#     marker='o',edgecolors='none',alpha=0.1)
# fig.savefig('_check_asassn4.pdf')


#ax.errorbar(tc_bing,fc_bing,yerr=fsigc_bing,label='g2',fmt='.')

binned_g_err = fsigc_bing/np.sqrt(np.sqrt(fc_npoig))

(cum_bing, xbing, cumbing) = hist_errors (binned_g_err,frac_keep=0.90)
print(cum_bing)


binned_V_err = fsigc_binV/np.sqrt(np.sqrt(fc_npoiV))

(cum_binV, xbinV, cumbinV) = hist_errors (binned_V_err,frac_keep=0.90)
print(cum_bing,cum_binV)


fig, (ax1,ax2) = plt.subplots(1,2,figsize=(12,6))
# cumulative plots....
ax1.plot(xbing,cumbing)
ax1.set_xlabel('g Flux err')
ax1.set_ylabel('Cumulative fraction')
ax2.plot(xbinV,cumbinV)
ax2.set_xlabel('V Flux err')
ax2.set_ylabel('Cumulative fraction')
ax2.set_title('data from {}'.format(fin))

# now filter binned g band flux on errors
mg = (binned_g_err<cum_bing)
mV = (binned_V_err<cum_binV)
fig, (ax) = plt.subplots(1,1,figsize=(12,6))
ax.set_ylabel('Flux [mJy]')
ax.set_xlabel('Epoch [MJD]')
ax.errorbar(tc_bing[mg],fc_bing[mg],yerr=binned_g_err[mg],label='g Band',fmt='.')
ax.errorbar(tc_binV[mV],fc_binV[mV],yerr=binned_V_err[mV],label='V band',fmt='.')

ax.set_title('ASASSN-21js light curve')
ax.legend()
fig.savefig('_check_asassn6.pdf')


plt.show()


tn = vstack([tV,tg])
tn['Survey'] = "ASASSN"

tn.write(paths.data / 'obs_ASASSN-21js_ASASSN.ecsv',
    format='ascii.ecsv',overwrite=True)
