from ztfquery import query,marshal
from daily_scan import *

def dump_rcf_table():
    m = marshal.MarshalAccess()
    m.load_target_sources("Redshift Completeness Factor")
    names = np.array(m.target_sources['name'])
    classes = np.array(m.target_sources['classification'])
    redshift = np.array(m.target_sources['redshift'])


s = logon()
dat = np.loadtxt("rcf.txt", delimiter=',', dtype=str)
names = dat[:,0]
classes = dat[:,1]
redshift = dat[:,2]


# Munging
cl = np.copy(classes)

# ignore sources with no known redshift
no_z = dat[:,2] == 'None'
cl[no_z] = ''

choose = np.array(['None' in val for val in classes])
cl[choose] = ''
choose = np.array(['AGN' in val for val in classes])
cl[choose] = ''
choose = np.array(['CV' in val for val in classes])
cl[choose] = ''
choose = np.array(['Duplicate' in val for val in classes])
cl[choose] = ''
choose = np.array(['Galaxy' in val for val in classes])
cl[choose] = ''
choose = np.array(['Q' in val for val in classes])
cl[choose] = ''
choose = np.array(['?' in val for val in classes])
cl[choose] = ''
choose = np.array(['Ic-BL' in val for val in classes])
cl[choose] = 'Ibc'
choose = np.array([val=='Ic' for val in classes])
cl[choose] = 'Ibc'
choose = np.array([val=='Ib' for val in classes])
cl[choose] = 'Ibc'
choose = np.array([val=='Ib/c' for val in classes])
cl[choose] = 'Ibc'

choose = np.array([val=='SN II' for val in classes])
cl[choose] = 'Type II'
choose = np.array([val=='SNII' for val in classes])
cl[choose] = 'Type II'
choose = np.array(['SN IIL' in val for val in classes])
cl[choose] = 'Type II'
choose = np.array(['SN IIP' in val for val in classes])
cl[choose] = 'Type II'

bad = np.array(['?' in val for val in classes])
cl[bad] = ''

search_terms = ['Ia', 'Type II', 'IIb', 'IIn', 'SLSN', 'Ibc', 'Gap', 'ILRT']

for st in search_terms[4:5]:
    print("Running for %s" %s)
    keep = np.array([st in val for val in cl])
    cands = names[keep]
    z = redshift[keep]

    ncands = len(cands)
    print("%s total objects" %ncands)
    print(cands)
    rises_all = np.zeros(ncands)
    fades_all = np.zeros(ncands)

    names_out = []
    z_out = []
    m_out = []
    rise_out = []
    fade_out = []

    for ii,cand in enumerate(cands):
        print(ii)
        jd,mag,emag,fid,program = get_lc(s, cand)
        choose = fid == 1 # g-band
        if sum(choose) > 5:
            jd_g = jd[choose]
            dt_g = jd_g-jd_g[0]
            mag_g = mag[choose]
            emag_g = emag[choose]
            
            # Check whether there is a peak
            indpeak = np.argmin(mag_g)
            mpeak = mag_g[indpeak]
            empeak = emag_g[indpeak]
            dmrise = mag_g[0]-mpeak
            dmfade = mag_g[-1]-mpeak
            rise_sig = dmrise/np.sqrt(empeak**2+emag_g[0]**2)
            fade_sig = dmfade/np.sqrt(empeak**2+emag_g[-1]**2)

            has_peak = np.logical_and.reduce((
                    mag_g[0] > mpeak + 0.75, mag_g[-1] > mpeak + 0.75,
                    rise_sig > 1, fade_sig > 1))

            # If peak, calculate rise time and fade time
            if has_peak:
                print("has peak")
                tpeak = dt_g[indpeak]

                # Rise Time
                dt_rise = dt_g[0:indpeak+1]
                mag_rise = mag_g[0:indpeak+1]
                emag_rise = emag_g[0:indpeak+1]
                order = np.argsort(mag_rise)
                dt_rise_sorted = dt_rise[order]
                mag_rise_sorted = mag_rise[order]
                tval = np.interp(
                        mpeak+0.75, mag_rise_sorted, dt_rise_sorted)
                rise_time = tpeak-tval

                # Fade Time
                dt_fade = dt_g[indpeak:]
                mag_fade = mag_g[indpeak:]
                emag_fade = emag_g[indpeak:]
                order = np.argsort(mag_fade)
                dt_fade_sorted = dt_fade[order]
                mag_fade_sorted = mag_fade[order]
                tval = np.interp(
                        mpeak+0.75, mag_fade_sorted, dt_fade_sorted)
                fade_time = tval-tpeak

                names_out.append(cand)
                z_out.append(z[ii])
                m_out.append(mpeak)
                rise_out.append(rise_time)
                fade_out.append(fade_time)

    np.savetxt(
        "%s.txt" %st, 
        np.array([names_out, z_out, m_out, rise_out, fade_out]).T, fmt='%s')
