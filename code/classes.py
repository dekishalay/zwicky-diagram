""" Events / parameter space for different class of transients """
import os
import sncosmo
import extinction
import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.cosmology import Planck15
from astropy.io import ascii
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/fast_transient_search/code")
from read_kann_lc import load_lc
sys.path.append("/Users/annaho/Dropbox/Projects/Research/Koala/code")

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
        
        tbnew = Table(data = [tb['name'].values,tb['z'].values, 
            Mpeaks, Mpeaks_unc, Tdeclines],
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
    #ax.errorbar(x, y, y_err, fmt='.k', ms=0.1) #, lw=0.1)
    ax.scatter(x, y, marker='.', s=1, c='k', zorder=1)
    
    # plot the grey contour
    slope, e_slope, intercept = mylinear_fit(x, y, y_err, npar=2)
    nstd = 5. # to draw 5-sigma intervals
    slope_up = slope + nstd * e_slope
    slope_dw = slope - nstd * e_slope
    x_fit = np.linspace(min(x)-0.5, max(x)+0.5, 100)
    # fit = slope * x_fit+ intercept
    fit_up = slope_up * x_fit+ intercept
    fit_dw = slope_dw * x_fit+ intercept
    ax.fill_between(
            x_fit, fit_up, fit_dw, facecolor = 'grey', 
            label="5-sigma interval",
            edgecolor='k', lw=1.0, zorder=0)
           

def fix_df(df):
    """ Helper function provided by Yashvi Sharma """
    newcol = []
    for i in range(len(df)):
        lcstr = df.loc[i]['lc']
        data = []
        lcstr = lcstr.split('\n')
        for line in lcstr[1:]:
            data.append(line.split())
        lcdf = pd.DataFrame(columns=['index','filter','jd','mag','magerr'],data=data)
        newcol.append(lcdf)
    df2 = df.drop(columns='lc')
    df2['lc'] = newcol
    return df2


def core_collapse(ax):
    """ Light curves of 230 CC SNe from RCF provided by Yashvi Sharma """
    # Load the table of events
    df = pd.read_csv("../data/ccsne/ccsne_final.csv")
    keep = df['redshift'] != 'None'

    nearby = df['redshift'][keep].astype(float)<0.05
    sampled = df['numpoints_lc'][keep]>50
    choose = np.logical_and(nearby, sampled)
    z = df['redshift'][keep].values[choose].astype(float)

    # Load the array of light curves
    lcs = np.load(
            "../data/ccsne/ccsne_lcs.npy", 
            allow_pickle=True, encoding='latin1')
    lcs = lcs[keep][choose]

    # Plot one of them
    for ii in np.arange(sum(choose)):
        lc = lcs[ii]
        filt = lc['filter']
        jd = lc['jd']
        mag = lc['mag']
        gband = filt == 'g'
        if sum(gband) > 30:
            x = jd[gband]-jd[gband].values[0]
            y = mag[gband]
            peakmag = np.min(mag)
            thalf = np.interp(peakmag+0.75, y, x)
            absmag = peakmag-Planck15.distmod(z=z[ii]).value
            ax.scatter(thalf, absmag, c='orange', marker='s', zorder=0)


def fbots(ax):
    """ Will include an inset showing fast-rising transients """
    # Koala
    x, y = 6/1.2714, -21.2
    ax.scatter(x, y, marker='D', c='k')
    #ax.text(x, y, "ZTF18abvkwla", fontsize=10, horizontalalignment='right')

    # PS1-10bjp
    x, y = (1+7.7)/1.113, -18.2
    ax.scatter(x, y, marker='D', c='k')

    # PS1-11qr
    x, y = (3+8.7)/1.324, -19.3
    ax.scatter(x, y, marker='D', c='k')

    # PS1-11qr
    x, y = (2+5)/1.405, -19.1
    ax.scatter(x, y, marker='D', c='k')

    # PS1-12brf
    x, y = (9)/1.275, -18.3
    ax.scatter(x, y, marker='D', c='k')


def gap(ax):
    """ Data points from prev version of diagram """
    dat = np.loadtxt("../data/gap.txt", delimiter=',')
    ax.scatter(dat[:,0], dat[:,1], marker='X', c='k')


def relativistic(ax):
    """ 
    Light curves of optical GRB afterglows provided by Anna Ho
    ah@astro.caltech.edu
    """

    datadir = "/Users/annaho/Dropbox/Projects/Research/fast_transient_search/data/lc"
    ptf11agg = ascii.read(datadir + "/lc_11agg.txt")
    agg_mjd = ptf11agg['MJD']
    agg_dt = agg_mjd-agg_mjd[0]
    agg_mag = ptf11agg['mag']
    agg_emag = ptf11agg['emag']
    z = 1.5
    peak = min(agg_mag)
    thalf = np.interp(peak+0.75, agg_mag, agg_dt)/(1+z)
    mag = peak-Planck15.distmod(z=z).value
    ax.text(thalf, mag*1.01, "PTF11agg", fontsize=10,
            horizontalalignment='left', verticalalignment='bottom')
    ax.scatter(thalf, mag, marker='o', c='k', zorder=1)

    # Temporary values from the previous version of the plot
    x,y = 0.010179047221370183, -28.38709677419355
    ax.scatter(x, y, marker='o', c='k', zorder=1)
    ax.text(x, y*1.01, "iPTF14cva")

    x,y = 0.011941841178367309, -27.81783681214422
    ax.scatter(x, y, marker='o', c='k', zorder=1)
    ax.text(x, y*1.01, "iPTF14cyb")

    x,y = 0.0420995502849078, -27.390891840607214
    ax.scatter(x, y, marker='o', c='k', zorder=1)
    ax.text(
            x, y, "iPTF14yb", 
            verticalalignment='top', horizontalalignment='right')

    #x,y = 0.06003714594928769, -27.01138519924099
    #ax.scatter(x, y, marker='o', c='k')
    #ax.text(x, y*1.01, "iPTF13ekl")

    x,y = 0.14841699095798874, -26.631878557874767
    ax.scatter(x, y, marker='o', c='k', zorder=1)
    ax.text(x, y/1.01, "ATLAS17aeu",
            verticalalignment='top', horizontalalignment='center')

    x,y = 0.28618535042609117, -26.29981024667932
    ax.scatter(x, y, marker='o', c='k', zorder=1)
    ax.text(x, y*1.01, "iPTF14aue")

    x,y = 0.47880000093436637, -25.872865275142317
    ax.scatter(x, y, marker='o', c='k', zorder=1)
    ax.text(x, y, "iPTF13dsw", verticalalignment='top',
            horizontalalignment='left')

    # Add Kann sample if I can
    lc = load_lc()
    for key,value in lc.items():
        try:
            z = float(value['z'])
            t = value['t'].astype(float)
            mag = value['mag'].astype(float)
            thalf = np.interp(min(mag)+0.75, mag, t)/(1+z)
            Mpeak = min(mag)-Planck15.distmod(z=z).value
            ax.scatter(thalf, Mpeak, marker='o', c='lightgrey', zorder=0)
        except:
            pass


def slsne(ax):
    """ Data points from Dan Perley """
    dat = np.loadtxt("../data/slsne/values.txt", delimiter=',')
    ax.scatter(dat[:,0], dat[:,1], marker='v', c='purple')


def novae(ax):
    # TBD
    pass


def ilrt(ax):
    """ Data points from prev version of diagram """
    dat = np.loadtxt("../data/ilrt.txt", delimiter=',')
    ax.scatter(dat[:,0], dat[:,1], marker='v', c='red')


def lrne(ax):
    """ Data points from prev version of diagram """
    dat = np.loadtxt("../data/lrn.txt", delimiter=',')
    ax.scatter(dat[:,0], dat[:,1], marker='v', c='red')


