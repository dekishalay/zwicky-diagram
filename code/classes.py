""" Events / parameter space for different class of transients """
import os
import sncosmo
import extinction
import numpy as np
import pandas as pd
from astropy.table import Table

# Define cosmology
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=69.6, Om0=0.286)

# Import plot functions
from matplotlib.patches import Rectangle


def snia(ax, band = 'p48g', timescale='efoldDecline',
         remove_nearby = True):
    """
    Please contact me if you have any questions: yyao@astro.caltech.edu
    
    Sample from SNe Ia young sample observed in the ZTF extragalactic high-cadence
    partnership experiment. We use 121 ``normal'' SNe Ia from 2018 (including 
    normal, 91T-like, and 91bg-like SNe).
    Light curves presented in Yao et al. (2019)
    
    Parameters:
    -----------
    band: 
        'p48g': use g band ZTF data
        'p48r': use r band ZTF data
        'swope2::v', 'f555w', ...: use the SALT2 (Guy et al 2007) tool to estimate light curves in other bands
    
    timescale: 
        'efoldDecline': plot the rest-frame duration for the flux to drop from peak to 1/e*peak
        [other options TBD]
        
    remove_nearby [bool]: 
        If True, remove two SNe with z < 0.02 (ZTF18aasdted and ZTF18abgmcmv),
        since their absolute peak luminosity are largely affected by the local 
        group peculiar velocity.
    """
    from snia_funcs import add_ZTFfilters, mylinear_fit
    
    cwd = os.getcwd()
    datafile = cwd[:-4]+'data/snia/' + band + '_' + timescale + '.csv'
    tbfile = cwd[:-4]+'data/snia/Yao2019_catalog.csv'
    try:
        os.stat(datafile)
    except:
        tb = pd.read_csv(tbfile)
        add_ZTFfilters()
        zp_temp = 26.325
        
        Mpeaks = np.zeros(len(tb))
        Mpeaks_unc = np.zeros(len(tb))
        Tdeclines = np.zeros(len(tb))
        
        for i in range(len(tb)):
            name = tb['name'].values[i]
            lcfile = cwd[:-4]+'data/snia/lightcurve/' + name + '.csv'
            z = tb['z'].values[i]
            z_unc = tb['z_unc'].values[i]
            ebv = tb['E_B_V_SandF'].values[i]
            
            mytb = pd.read_csv(lcfile)
            lcband = mytb['band'].values
            Fratio = mytb['Fratio'].values
            Fratio_unc = mytb['Fratio_unc'].values
            zp = np.ones(len(mytb)) * zp_temp
            zpsys = ['ab']*len(mytb)
            time = mytb['jdobs'].values
            lc = Table(data = [time, lcband, Fratio, Fratio_unc, zp, zpsys],
                         names = ['time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'])

            dust = sncosmo.CCM89Dust()
            model = sncosmo.Model(source='salt2', effects=[dust],
                                  effect_names=['mw'], effect_frames=['obs'])
            model.set(z=z, mwebv= ebv, mwr_v=3.1)
            res, fitted_model = sncosmo.fit_lc(lc, model, 
                                               ['t0', 'x0', 'x1', 'c'],
                                               bounds={'x1':(-10.,10.), 'c':(-5.,5.)})
            
            t0jd_B = res.parameters[1] # B band maximum light epoch
            t0jd_B_unc = res['errors']['t0']
            
            ta = np.round(t0jd_B-25, 2)+0.01
            tgrid = np.linspace(ta, ta+75, 100*100+1)
            fitted_band = fitted_model.bandflux(band, tgrid, zp=zp_temp, zpsys='ab')
    
            # estimate peak epoch and apparent magnitude
            ixx = np.argsort(fitted_band)
            t_max = tgrid[ixx[-1]]
            fratio_max = fitted_band[ixx[-1]]
            appmag_max = -2.5 * np.log10(fratio_max) 
            appmag_around = fitted_model.bandmag(band, 'ab', [t_max - t0jd_B_unc, 
                                                              t_max + t0jd_B_unc]) - zp_temp
            appmag_max_unc = (np.sum(appmag_around)-2*appmag_max)/2.
            
            # estimate extinction corrected peak epoch apparent magnitude
            wave_eff = sncosmo.get_bandpass(band).wave_eff
            Aband = extinction.ccm89(np.array([wave_eff]), 3.1*ebv, 3.1)[0]
            appmag0_max = appmag_max - Aband
            
            # estimate extinction corrected peak epoch absolute magnitude
            D = cosmo.luminosity_distance([z])[0].value * 1e+6 # in pc
            D_unc = (cosmo.luminosity_distance([z+z_unc])[0].value - cosmo.luminosity_distance([z-z_unc])[0].value) / 2. * 1e+6
            dis_mod = 5*np.log10(D / 10) 
            dis_mod_unc = 5 / np.log(10) / D * D_unc
            absmag0_max = appmag0_max - dis_mod
            absmag0_max_unc = np.sqrt(appmag_max_unc**2 + dis_mod_unc**2)
            Mpeaks[i] = absmag0_max
            Mpeaks_unc[i] = absmag0_max_unc
            
            # estimate the characteristic timescale
            if timescale=='efoldDecline':
                fratio_end = fratio_max/np.e
                ind = np.argsort(abs(fitted_band - fratio_end))
                ix_end = ind[ind > ixx[-1]][0]
                t_end = tgrid[ix_end]
                time_duration = t_end - t_max
                Tdeclines[i] = time_duration
        
        tbnew = Table(data = [tb['name'].values,tb['z'].values, Mpeaks, Mpeaks_unc, Tdeclines],
                      names = ['name', 'z', 'mag', 'mag_unc', 'timescale'])
        tbnew.write(datafile, overwrite=True)
        
    data = pd.read_csv(datafile)
    data = data[data.z > 0.02]
    x = data['timescale'].values / (1 + data['z'].values)
    y = data['mag'].values
    y_err = data['mag_unc'].values
    
    # plt.figure(figsize=(12, 8))
    # ax = plt.subplot(111)
    
    # plot individual data points
    ax.errorbar(x, y, y_err, fmt='.k')
    
    # plot the grey contour
    slope, e_slope, intercept = mylinear_fit(x, y, y_err, npar=2)
    nstd = 5. # to draw 5-sigma intervals
    slope_up = slope + nstd * e_slope
    slope_dw = slope - nstd * e_slope
    x_fit = np.linspace(min(x)-0.5, max(x)+0.5, 100)
    # fit = slope * x_fit+ intercept
    fit_up = slope_up * x_fit+ intercept
    fit_dw = slope_dw * x_fit+ intercept
    ax.fill_between(x_fit, fit_up, fit_dw, color = 'grey', alpha=.25, label="5-sigma interval")
           


def core_collapse(ax):
    # Yashvi
    pass


def gap(ax):
    # Kishalay
    pass


def relativistic(ax):
    # Anna
    pass


def slsne(ax):
    # Lin
    pass


def novae(ax):
    # TBD
    pass


def lrne(ax):
    # Nadia
    pass


