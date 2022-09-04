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
