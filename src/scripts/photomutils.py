import numpy as np

def hist_errors (e,frac_keep=0.95):
    '''hist_errors -  show histogram and CMD for photometric errors and return a cropped sample'''
    x = np.arange(100)
    cum = np.percentile(e,x)

    cum_req = np.percentile(e,frac_keep*100)

    return cum_req, x, cum

def mean_rms_region (t,f,e, tmin, tmax):
    '''mean_rms_region - determine the mean and r.m.s. of photometric points between tmin and tmax'''
    m = (t>tmin) * (t<tmax)

    if (np.sum(m) < 1):
        print(f'Warning! No points selected between {tmin} and {tmax}')

    f_mean = np.average(f[m], weights=(1./(np.power(e[m],2))))

    f_std = np.std(f[m])

    return f_mean, f_std



def rebinsig(t, f, tedges, ax=None, err=None):
    from astropy.stats import sigma_clip
    import matplotlib.patches as patches

    tc_bin = np.array([])
    fc_bin = np.array([])
    fsigc_bin = np.array([])
    fc_npoi = np.array([])

    t_lower = tedges[:-1]
    t_upper = tedges[1:]
    t_midpoint = (t_lower+t_upper)/2.

    for t_mid, t_low, t_high in zip(t_midpoint, t_lower, t_upper):
        # select points from flux that have
        selpoints = ((t > t_low) * (t < t_high))
        fluxsel = f[selpoints]
        timesel = t[selpoints]

        # if there's just one point, just copy the values over
        if (fluxsel.size == 1):
        #    print('time {} has {} points'.format(t_mid, fluxsel.size))
            tc_bin = np.append(tc_bin, timesel)
            fc_bin = np.append(fc_bin, fluxsel)
            if err is not None:
                fsigc_bin = np.append(fsigc_bin, err[selpoints].value)
            else:
                fsigc_bin = np.append(fsigc_bin, 0.)

            fc_npoi = np.append(fc_npoi, 1.)

       # if there are two points, don't run any sigma clipping
        if (fluxsel.size == 2):
        #    print('time {} has {} points'.format(t_mid, fluxsel.size))
            tc_bin = np.append(tc_bin, np.mean(timesel))
            fc_bin = np.append(fc_bin, np.mean(fluxsel))
            if err is not None:
                fsigc_bin = np.append(fsigc_bin, np.sqrt(np.sum(np.power(err[selpoints],2))))
            else:
                fsigc_bin = np.append(fsigc_bin, np.ma.std(fluxsel))

            fc_npoi = np.append(fc_npoi, 2.)

        # more than two points, go sigma clipping
        if (fluxsel.size > 2):
        #    print('time {} has {} points'.format(t_mid, fluxsel.size))

            meanflux = sigma_clip(fluxsel, sigma=3, maxiters=1, masked=True)
            meane = np.ma.mean(meanflux)
            # NOTE!!! should select the unmasked times and take the average of that, maybe???
            if np.isfinite(meane):
                tc_bin = np.append(tc_bin, np.ma.mean(timesel))
                fc_bin = np.append(fc_bin, np.ma.mean(meanflux))
                fsigc_bin = np.append(fsigc_bin, np.ma.std(meanflux))
                fc_npoi = np.append(fc_npoi, np.sum(selpoints))

        if (ax is not None) and (fluxsel.size > 0):
            f_low = np.min(fluxsel)
            f_hig = np.max(fluxsel)
            rect = patches.Rectangle((t_low,f_low), t_high-t_low, f_hig-f_low, edgecolor='r',facecolor='none')
            ax.add_patch(rect)

    return (tc_bin, fc_bin, fsigc_bin, fc_npoi)
