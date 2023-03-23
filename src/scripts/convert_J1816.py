import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import unique,vstack,Table
import paths

from photomutils import *

fin='J1816/light_curve_b1c2ca34-b951-4bf6-94aa-443a8dfe053c.csv'
t = ascii.read(paths.data / fin)

#      HJD           UT Date       Camera FWHM Limit   mag   mag_err flux(mJy) flux_err Filter
# ------------- ------------------ ------ ---- ------ ------ ------- --------- -------- ------
# 2457420.65322 2016-02-02.1500246     be 1.46 17.458  13.45   0.005    15.995     0.08      V

t['MJD'] = t['HJD']-2400000.5
obj='J1816'
flux_high = 75
# reject noisy points
t = t[(t['flux(mJy)']<flux_high)]
print(f'rejected ASASSN points with high fluxes (> {flux_high} mJy) in both bands')

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


# building histograms of errors and doing a percentage cut in the cumulative distributions
(cum_reqV, xV, cumV) = hist_errors (tV['flux_err'])
(cum_reqg, xg, cumg) = hist_errors (tg['flux_err'])
print(f'cutoff for error in V:{cum_reqV:5.2f} and g:{cum_reqg:5.2f}')

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

# reject noisy points
tV = tV[(tV['flux_err']<cum_reqV)]
tg = tg[(tg['flux_err']<cum_reqg)]
print(f'rejecting noisy points in g (>{cum_reqg:5.3f}) and V (>{cum_reqV:5.3f}) data separately in ASASSN')

print('binning photometry that is taken within the same night in g and V')
fig, (ax) = plt.subplots(1,1,figsize=(12,6))
deltaV = np.ediff1d(tV['MJD'])
ax.hist(deltaV,bins=1000,range=(-1,3))

# get indices of where delta > 0.5, that's a night boundary
big_deltaV = (deltaV>0.5)
(nz) = np.flatnonzero(big_deltaV)
time_lower = tV['MJD'][nz]
time_upper = tV['MJD'][nz+1]
time_mean = (time_lower+time_upper)/2.

tmin = np.min(tV['MJD'])
tmax = np.max(tV['MJD'])
tedges = np.concatenate((np.array([tmin-0.5]),time_mean,np.array([tmax+0.5])))

fig, (ax) = plt.subplots(1,1,figsize=(12,6))

(tc_binV, fc_binV, fsigc_binV, fc_npoiV) = \
    rebinsig(tV['MJD'],tV['flux(mJy)'], tedges, ax=ax, err=tV['flux_err'])

# get indices of where delta > 0.5, that's a night boundary
deltag = np.ediff1d(tg['MJD'])
big_deltag = (deltag>0.5)
(nz) = np.flatnonzero(big_deltag)
time_lower = tg['MJD'][nz]
time_upper = tg['MJD'][nz+1]
time_mean = (time_lower+time_upper)/2.

tmin = np.min(tg['MJD'])
tmax = np.max(tg['MJD'])
tedges = np.concatenate((np.array([tmin-0.5]),time_mean,np.array([tmax+0.5])))

(tc_bing, fc_bing, fsigc_bing, fc_npoig) = \
    rebinsig(tg['MJD'],tg['flux(mJy)'], tedges, ax=ax, err=tg['flux_err'])

ax.set_ylabel('Flux [mJy]')
ax.set_xlabel('Epoch [MJD]')
ax.set_title('data from {}'.format(fin))
ax.errorbar(tV['MJD'],tV['flux(mJy)'],yerr=tV['flux_err'],label='V',fmt='.')
ax.errorbar(tg['MJD'],tg['flux(mJy)'],yerr=tg['flux_err'],label='g',fmt='.')

ax.errorbar(tc_binV,fc_binV,yerr=fsigc_binV,label='V binned',fmt='.')

ax.errorbar(tc_bing,fc_bing,yerr=fsigc_bing,label='g binned',fmt='.')

fig.savefig('_check_asassn3.pdf')

fig, (ax) = plt.subplots(1,1,figsize=(12,6))

ax.set_ylabel('Flux [mJy]')
ax.set_xlabel('Epoch [MJD]')
ax.set_title('data from {}'.format(fin))

ax.errorbar(tc_binV,fc_binV,yerr=fsigc_binV/np.sqrt(fc_npoiV),label='V binned',fmt='.')
ax.errorbar(tc_bing,fc_bing,yerr=fsigc_bing/np.sqrt(fc_npoig),label='g binned',fmt='.')
fig.savefig('_check_asassn4.pdf')

binned_g_err = fsigc_bing/np.sqrt(np.sqrt(fc_npoig))

(cum_bing, xbing, cumbing) = hist_errors (binned_g_err,frac_keep=0.95)

binned_V_err = fsigc_binV/np.sqrt(np.sqrt(fc_npoiV))

(cum_binV, xbinV, cumbinV) = hist_errors (binned_V_err,frac_keep=0.95)
print(f'rejecting noisy points in binned g (>{cum_bing:5.3f}) and binned V (>{cum_binV:5.3f}) data separately in ASASSN')


fig, (ax1,ax2) = plt.subplots(1,2,figsize=(12,6))
# cumulative plots....
ax1.plot(xbing,cumbing)
ax1.set_xlabel('g Flux err')
ax1.set_ylabel('Cumulative fraction')
ax2.plot(xbinV,cumbinV)

ax2.set_xlabel('V Flux err')
ax2.set_ylabel('Cumulative fraction')
ax2.set_title('data from {}'.format(fin))

# now filter binned g and V band flux on errors
mg = (binned_g_err<cum_bing)
mV = (binned_V_err<cum_binV)

fig, (ax) = plt.subplots(1,1,figsize=(12,6))
ax.set_ylabel('Flux [mJy]')
ax.set_xlabel('Epoch [MJD]')
ax.errorbar(tc_bing[mg],fc_bing[mg],yerr=binned_g_err[mg],label='g Band',fmt='.')
ax.errorbar(tc_binV[mV],fc_binV[mV],yerr=binned_V_err[mV],label='V band',fmt='.')

ax.set_title('{obj} binned and filtered light curve')
ax.legend()
fig.savefig('_check_asassn6.pdf')


fig, (ax) = plt.subplots(1,1,figsize=(12,6))



tminV, tmaxV = 58150,58450

(V_flux_norm, V_flux_norm_err)= \
    mean_rms_region(tc_binV[mV], fc_binV[mV], binned_V_err[mV], tminV, tmaxV)

tming, tmaxg = 58150,58450

(g_flux_norm, g_flux_norm_err)= \
    mean_rms_region(tc_bing[mg], fc_bing[mg], binned_g_err[mg], tming, tmaxg)

print(f'Flux V band normalised {V_flux_norm:5.2}\pm{V_flux_norm_err:5.2f}')
print(f'Flux g band normalised {g_flux_norm:5.2}\pm{g_flux_norm_err:5.2f}')

tVout = Table([tc_binV[mV],fc_binV[mV]/V_flux_norm,
    binned_V_err[mV]/V_flux_norm],
    names=('MJD','fnorm','fnormerr'))
tVout['Filter'] = 'V'

tgout = Table([tc_bing[mg],fc_bing[mg]/g_flux_norm,
    binned_g_err[mg]/g_flux_norm],
    names=('MJD','fnorm','fnormerr'))
tgout['Filter'] = 'g'

fig, (ax) = plt.subplots(1,1,figsize=(12,6))
ax.set_ylabel('Flux [normalised]')
ax.set_xlabel('Epoch [MJD]')
ax.set_title('Normalised ASASSN data {}'.format(fin))

ax.errorbar(tgout['MJD'],tgout['fnorm'],yerr=tgout['fnormerr'],label='g Band',fmt='.')
ax.errorbar(tVout['MJD'],tVout['fnorm'],yerr=tVout['fnormerr'],label='V band',fmt='.')

ax.legend()

fig.savefig('_check_asassn9.pdf')

tn = vstack([tVout,tgout])
tn['Survey'] = "ASASSN"
tn['Source'] = obj
#plt.show()
tn.write(paths.data / 'obs_J1816_ASASSN.ecsv',
    format='ascii.ecsv',overwrite=True)
