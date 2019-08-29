""" 
For all sources in RCF,
calculate the rise time and fade time
(from half-max to max)
and save to a big table
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import subprocess
import pandas as pd
import glob
from penquins import Kowalski
from astropy.time import Time


def logon():
    """ Log onto Kowalski """
    username = 'ah'
    password = 'TetraodonInsists'
    s = Kowalski(
        protocol='https', host='kowalski.caltech.edu', port=443,
        verbose=False, username=username, password=password)
    return s


def download_obs_table(datestr):
    """ download the obs table from night of scanning """
    # Import the log from that night, to see which fields were observed
    # as part of these programs
    path = 'exp.%s.tbl' %datestr
    check = glob.glob(path)
    if len(check) == 0:
        subprocess.call(['wget', '--user', 'ztf', '--password', 'discover', 
            'http://www.astro.caltech.edu/~tb/ztfops/sky/%s/%s' %(datestr,path)])


def read_obs_table(datestr):
    # Skip the lines that contain the repeated header
    skip = []
    path = 'exp.%s.tbl' %datestr
    with open(path, "r") as inputf:
        all_lines = np.array(inputf.readlines())
        nlines = len(all_lines)
        to_skip = np.array(['UT_START' in val for val in all_lines])
        ind_to_skip = np.arange(0,nlines)[to_skip] 

    obs = pd.read_table(
            path, skiprows=ind_to_skip, delim_whitespace=True, header=None, 
            names=['UT_START', 'sun', 'exp', 'fil', 'type', 'field', 'pid', 
                'ra', 'dec', 'slew', 'wait', 'fileroot', 'progpi', 'qcomment'])

    return obs


def count_obs(obs, program):
    """ Count the number of observations taken
    in a given field, in each filter """
    choose = obs['qcomment'] == program
    fields = np.unique(obs['field'][choose])
    counts = np.zeros((len(fields), 2))

    for ii,fieldid in enumerate(fields):
        filts = [1,2]
        counts[ii,:] = np.array(
            [np.sum(np.logical_and(obs['field'][choose]==fieldid, obs['fil'][choose]==f)) for f in filts])
        
    return counts


def diagnostics(obs, datestr):
    """ For the night requested, output a histogram
    of the number of fields and how many visits were obtained
    per filter and per field """
    counts_partnership = count_obs(obs, 'high_cadence')
    counts_caltech = count_obs(obs, 'Caltech_1DC')

    # Make a histogram reporting the results
    fig = plt.figure(figsize=(8,2), dpi=120)

    # Partnership results
    ax1 = fig.add_subplot(121)
    x = np.arange(5)
    crit1 = sum(
        np.logical_and(counts_partnership[:,0]==2, counts_partnership[:,1]==1))
    crit2 = sum(np.sum(counts_partnership, axis=1)==2)
    crit3 = sum(np.sum(counts_partnership, axis=1)==1)
    crit4 = sum(
        np.logical_and(counts_partnership[:,0]==3, counts_partnership[:,1]==0))
    crit5 = sum(
        np.logical_and(counts_partnership[:,0]==3, counts_partnership[:,1]==1))
    ax1.bar(x, [crit4,crit5,crit1,crit2,crit3])
    ax1.set_xticklabels(['','3g1r','3g0r','2g,1r', '2g or 1g1r', '1g or 1r'])
    ax1.set_ylabel("Number of fields")
    ax1.set_title("Partnership (%s fields)" %counts_partnership.shape[0])

    # Caltech results
    ax2 = fig.add_subplot(122)
    x = np.arange(2)
    crit1 = sum(np.logical_and(counts_caltech[:,0]==1, counts_caltech[:,1]==1))
    crit2 = sum(np.sum(counts_caltech, axis=1)==1)
    ax2.bar(x, [crit1,crit2])
    ax2.set_xticklabels(['','1g,1r','', '1g or 1r'])
    ax2.set_ylabel("Number of fields")
    ax2.set_title("Caltech 1DC (%s fields)" %counts_caltech.shape[0])

    plt.savefig("diagnostics_%s.png" %datestr, dpi=500)

    return obs


def initial_search(s, obs):
    """ Run an initial broad query """

    # Set search window
    obst = obs['UT_START'].values
    start = Time(obst[0], format='isot').jd
    end = Time(obst[-1], format='isot').jd

    choose = np.logical_or(
            obs['qcomment'] == 'high_cadence', 
            obs['qcomment'] == 'Caltech_1DC')
    fields = [int(val) for val in np.unique(obs['field'][choose])]

    # Note: ndethist includes the current detection. 
    # It's never 0. For some reason using it > 1 causes the query to take
    # forever, so I'm just going to filter on that afterwards.

    q = {"query_type": "find",
         "query": {
             "catalog": "ZTF_alerts",
             "filter": {
                     'candidate.jd': {'$gt': start, '$lt': end},
                     'candidate.field': {'$in': fields},
                     'candidate.ndethist': {'$gt': 0},
                     'candidate.magpsf': {'$lt': 20},
                     'candidate.isdiffpos': {'$in': ['1', 't']},
                     '$expr': {'$lt': [{'$subtract': ['$candidate.jd', '$candidate.jdstarthist']}, 30]},
                     'candidate.programid': {'$gt': 1},
                     'candidate.ssdistnr': {'$lt': -1},
                     'candidate.drb': {'$gt': 0.5},
             },
             "projection": {
                     "_id": 0,
                     "objectId": 1,
             }
         }
         }

    query_result = s.query(query=q)
    out = query_result['result_data']['query_result']

    names = np.unique([val['objectId'] for val in out])
    print("There are %s unique cands from this initial filter." %len(names))

    return names


def star_check(s, name):
    """ Remove stars from the list """
    query = """db['ZTF_alerts'].find_one(
            {'objectId': {'$eq': '%s'}},
            {'candidate.distpsnr1', 'candidate.sgscore1', 'candidate.srmag1', 'candidate.sgmag1', 'candidate.simag1', 'candidate.szmag1'}
        )""" %name
    query_result = s.query(
            {'query_type': 'general_search', 'query': query})
    out = query_result['result_data']['query_result']
    dist = out['candidate']['distpsnr1']
    sg = out['candidate']['sgscore1']
    rmag = out['candidate']['srmag1']
    gmag = out['candidate']['sgmag1']
    imag = out['candidate']['simag1']
    zmag = out['candidate']['szmag1']

    pointunderneath = False
    if np.logical_and.reduce((dist>0, dist<2, sg>0.76)):
        pointunderneath = True
    if np.logical_and.reduce((sg==0.5, dist>0, dist<0.5)):
        if np.logical_or.reduce((
            np.logical_and(rmag>0,rmag<17),
            np.logical_and(gmag>0,gmag<17),
            np.logical_and(imag>0,imag<17),
            np.logical_and(zmag>0,zmag<17),
            )):
            pointunderneath = True

    brightstar = False
    if np.logical_and.reduce((rmag > 0, rmag < 12, sg > 0.49, dist < 20)):
        brightstar = True
    if np.logical_and.reduce((rmag > 0, rmag < 15, sg > 0.8, dist < 20)):
        brightstar = True
    if np.logical_and.reduce((gmag > 0, gmag < 12, sg > 0.49, dist < 20)):
        brightstar = True
    if np.logical_and.reduce((gmag > 0, gmag < 15, sg > 0.8, dist < 20)):
        brightstar = True
    if np.logical_and.reduce((imag > 0, imag < 12, sg > 0.49, dist < 20)):
        brightstar = True
    if np.logical_and.reduce((imag > 0, imag < 15, sg > 0.8, dist < 20)):
        brightstar = True
    if np.logical_and.reduce((zmag > 0, zmag < 12, sg > 0.49, dist < 20)):
        brightstar = True
    if np.logical_and.reduce((zmag > 0, zmag < 15, sg > 0.8, dist < 20)):
        brightstar = True

    return np.logical_or(pointunderneath, brightstar)


def get_dets(s, name):
    q = {"query_type": "find",
         "query": {
             "catalog": "ZTF_alerts",
             "filter": {
                     'objectId': {'$eq': name},
                     'candidate.isdiffpos': {'$in': ['1', 't']},
             },
             "projection": {
                     "_id": 0,
                     "prv_candidates": 1,
                     "candidate.jd": 1,
                     "candidate.magpsf": 1,
                     "candidate.sigmapsf": 1,
                     "candidate.fid": 1,
                     "candidate.pid": 1,
             }
         }  
         }  
    query_result = s.query(query=q)
    out = query_result['result_data']['query_result']
    return out


def get_lc(s, name):
    """ Retrieve LC for object """
    out = get_dets(s, name)
    jd = []
    mag = []
    emag = []
    filt = []
    program = []
    
    for det in out:
        cand = det['candidate']

        current_jd = cand['jd']
        current_mag = cand['magpsf']
        current_emag = cand['sigmapsf']
        current_filter = cand['fid']
        current_prog = cand['pid']
        
        if current_jd not in jd:
            jd.append(current_jd)
            mag.append(current_mag)
            emag.append(current_emag)
            filt.append(current_filter)
            program.append(current_prog)
        
        if det['prv_candidates'] != None:
            for prv_cand in det['prv_candidates']:
                if prv_cand['magpsf'] != None:
                    if prv_cand['jd'] not in jd:
                        jd.append(prv_cand['jd'])
                        mag.append(prv_cand['magpsf'])
                        emag.append(prv_cand['sigmapsf'])
                        filt.append(prv_cand['fid'])
                        program.append(prv_cand['pid'])

    jd = np.array(jd)
    mag = np.array(mag)
    emag = np.array(emag)
    filt = np.array(filt)
    program = np.array(program)

    # Sort in order of jd
    order = np.argsort(jd)
    jd = jd[order]
    mag = mag[order]
    emag = emag[order]
    filt = filt[order]
    program = program[order]
    return jd,mag,emag,filt,program


def plot_band(ax, dt,mag,emag,mask,col):
    indpeak = np.argmin(mag[mask])
    tpeak = dt[mask][indpeak]
    ax.errorbar(dt[mask], mag[mask], emag[mask], c=col, fmt='o')
    ax.axhline(y=mag[mask][indpeak], c=col, lw=0.5)
    ax.axvline(x=tpeak, c=col, lw=0.5)
    ax.axhline(y=mag[mask][indpeak]+0.75, c=col, lw=0.5)


def plot_lc(ax,name,dt,mag,emag,filt):
    if sum(filt==1) > 0:
        plot_band(ax, dt, mag, emag, filt==1, 'green')
    if sum(filt==2) > 0:
        plot_band(ax, dt, mag, emag, filt==2, 'red')
    if sum(filt==3) > 0:
        plot_band(ax, dt, mag, emag, filt==3, 'grey')

    ax.set_xlabel("Days", fontsize=14)
    ax.set_ylabel("App Mag", fontsize=14)
    ax.invert_yaxis()
    ax.tick_params(axis='both', labelsize=14)


def plotlims_band(ax, dt,mag,mask,col):
    ax.scatter(
            dt[mask], mag[mask], edgecolor=col, facecolor='white', marker='v')


def plot_lims(ax, limdt,limmag,limfilt):
    plotlims_band(ax, limdt, limmag, limfilt==1, 'green')
    plotlims_band(ax, limdt, limmag, limfilt==2, 'red')


def get_nondet(s, name):
    q = {"query_type": "find",
         "query": {
             "catalog": "ZTF_alerts",
             "filter": {
                     'objectId': {'$eq': name}
             },
             "projection": {
                     "_id": 0,
                     "prv_candidates": 1,
             }
         }  
         }  
    query_result = s.query(query=q)
    out = query_result['result_data']['query_result']
    
    limjds = []
    limmags = []
    limfilts = []
    for res in out:
        if res['prv_candidates'] is not None:
            for det in res['prv_candidates']:
                if det['magpsf'] is None:
                    limjds.append(det['jd'])
                    limmags.append(det['diffmaglim'])
                    limfilts.append(det['fid'])
    
    limjds = np.array(limjds)
    limmags = np.array(limmags)
    limfilts = np.array(limfilts)    
    
    order = np.argsort(np.array(limjds))
    limjds = limjds[order]
    limmags = limmags[order]
    limfilts = limfilts[order]
    
    return limjds,limmags,limfilts


def calc_trise(dt,mag,mask):
    if len(dt[mask]) > 0:
        indpeak = np.argmin(mag[mask])
        tpeak = dt[mask][indpeak]
        mpeak = mag[mask][indpeak]
        dt_fit = dt[mask][0:indpeak+1]
        mag_fit = mag[mask][0:indpeak+1]
        order = np.argsort(mag_fit)
        dt_fit = dt_fit[order]
        mag_fit = mag_fit[order]
        tval = np.interp(mpeak+0.75, mag_fit, dt_fit, left=-999, right=-999)
        return tpeak-tval
    else:
        return 1000
    
    
def calc_tfade(dt,mag,mask):
    if len(dt[mask]) > 0:
        indpeak = np.argmin(mag[mask])
        tpeak = dt[mask][indpeak]
        mpeak = mag[mask][indpeak]
        dt_fit = dt[indpeak:]
        mag_fit = mag[indpeak:]
        order = np.argsort(mag_fit)
        dt_fit = dt_fit[order]
        mag_fit = mag_fit[order]
        tval = np.interp(mpeak+0.75, mag_fit, dt_fit, left=-999, right=-999)
        return tpeak-tval
    else:
        return 1000
    
    
def calc_dmdt(dt,mag,emag):
    dm = np.array(mag[1:]-mag[0:-1])
    dt = np.array(dt[1:]-dt[0:-1])
    dmdt = np.abs(dm/dt)
    sig = np.array(np.abs(dm)/np.sqrt(emag[1:]**2 + emag[0:-1]**2))
    # only pay attention to dm/dt values that are significant
    choose = sig > 1
    if sum(choose) > 0:
        return np.max(dmdt[choose])
    return 0


def plot_source_summary(s, name):
    # Get LC -- dets and non-dets separately
    jd,mag,emag,fid,program = get_lc(s, name)
    if jd[-1]-jd[0] > 0.02:
        limjds,limmags,limfilts = get_nondet(s, name)

        # Decide on a t0 for the plot
        t0 = jd[np.argmin(mag)]

        # Plot the LC
        fig = plt.figure(figsize=(8,10))
        gs = gridspec.GridSpec(4, 2, figure=fig, hspace=0.4, wspace=0.4)

        # Plot the FULL LC
        ax = plt.subplot(gs[0,:]) # first row, all columns
        plot_lc(ax,name,jd-t0,mag,emag,fid)
        plot_lims(ax,limjds-t0,limmags,limfilts)  
        ax.set_title("%s" %(name), fontsize=14)

        # Plot a 20-day zoom-in of the LC
        ax = plt.subplot(gs[1,:]) 
        plot_lc(ax,name,jd-t0,mag,emag,fid)
        plot_lims(ax,limjds-t0,limmags,limfilts)  
        # somewhat of a zoom-in
        ax.set_xlim(-10,10)

        # Calculate some basic numbers for each filter
        for fidx in [1, 2]:
            choose = fid == fidx
            if sum(choose) > 1:
                jd_fit = jd[choose]
                dt_fit = jd_fit-t0
                mag_fit = mag[choose]
                emag_fit = emag[choose]
                ndet = len(jd_fit)

                choose_lim = limfilts == fidx
                limdt_fit = limjds[choose_lim]-t0
                limmags_fit = limmags[choose_lim]
                limfilts_fit = limfilts[choose_lim]

                ax = plt.subplot(gs[fidx,0])
                plot_lc(
                        ax,name,dt_fit,mag_fit,emag_fit,
                        np.array([fidx]*ndet))
                plot_lims(ax,limdt_fit,limmags_fit,limfilts_fit)  
                ax.set_xlim(-20,20)

                # Calculate the duration
                src_duration = dt_fit[-1]-dt_fit[0]
                ax.text(
                        0.9, 0.1, 
                        "Duration [d]: %s" %np.round(src_duration, 2),
                        transform=ax.transAxes, horizontalalignment='right')

                # Calculate the max dmag/dt of the detections
                dmdt_det = calc_dmdt(dt_fit, mag_fit, emag_fit)

                # If it exists,
                # calculate the dmag/dt of the last non-detection
                # to the first detection
                choose_lim = limfilts == fidx
                before_first_det = limjds[choose_lim] < jd_fit[0]
                if sum(before_first_det) > 0:
                    lastlim_jd = limjds[choose_lim][before_first_det][-1] 
                    lastlim_dt = lastlim_jd - t0
                    lastlim_mag = limmags[choose_lim][before_first_det][-1] 
                    dmdt_lim = calc_dmdt(
                            np.array([lastlim_dt, dt_fit[0]]), 
                            np.array([lastlim_mag, mag_fit[0]]),
                            np.array([0, emag_fit[0]]))

                    # Keep whichever is bigger,
                    # assuming that the limiting mag is actually
                    # fainter than the detection
                    if np.logical_and(dmdt_lim>dmdt_det,
                            lastlim_mag > mag_fit[0]):
                        dmdt = dmdt_lim
                    else:
                        dmdt = dmdt_det
                else:
                    dmdt = dmdt_det
                ax.text(
                        0.90, 0.05, "dm/dt [mag/d]: %s" %np.round(dmdt, 2),
                        transform=ax.transAxes, horizontalalignment='right')

                # Check whether there is a peak
                indpeak = np.argmin(mag_fit)
                mpeak = mag_fit[indpeak]
                empeak = emag_fit[indpeak]
                dmrise = mag_fit[0]-mpeak
                dmfade = mag_fit[-1]-mpeak
                rise_sig = dmrise/np.sqrt(empeak**2+emag_fit[0]**2)
                fade_sig = dmfade/np.sqrt(empeak**2+emag_fit[-1]**2)

                has_peak = np.logical_and.reduce((
                        mag_fit[0] > mpeak, mag_fit[-1] > mpeak,
                        rise_sig > 1, fade_sig > 1))

                # If the source has peaked, check the rise time
                ax = plt.subplot(gs[fidx,1])
                tpeak = dt_fit[indpeak]
                dt_rise = dt_fit[0:indpeak+1]
                mag_rise = mag_fit[0:indpeak+1]
                emag_rise = emag_fit[0:indpeak+1]
                ndet = len(dt_rise)
                plot_lc(
                        ax,name,dt_rise,mag_rise,emag_rise,
                        np.array([fidx]*ndet))
                ax.set_xlim(min(dt_rise), max(dt_rise))
                if has_peak:
                    order = np.argsort(mag_rise)
                    dt_rise_sorted = dt_rise[order]
                    mag_rise_sorted = mag_rise[order]
                    tval = np.interp(
                            mpeak+0.75, mag_rise_sorted, dt_rise_sorted, 
                            left=-999, right=-999)
                    rise_time = tpeak-tval
                    if rise_time > 998:
                        dur_rise = dt_rise[-1]-dt_rise[0]
                        ax.text(
                                0.90, 0.05, 
                                "$t_{1/2}$ [d] > %s" %np.round(dur_rise,2),
                                transform=ax.transAxes, 
                                horizontalalignment='right')
                    else:
                        ax.text(
                                0.90, 0.05, 
                                "$t_{1/2}$ [d]: %s" %np.round(rise_time,2),
                                transform=ax.transAxes, 
                                horizontalalignment='right')
                else:
                    ax.text(
                            0.90, 0.05, 
                            "No significant peak", 
                            transform=ax.transAxes, 
                            horizontalalignment='right')


        # Save figure
        # plt.show()
        plt.savefig("%s.png" %(name))
        plt.close()

