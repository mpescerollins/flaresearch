import matplotlib.pyplot as plt
import os, sys
import numpy as np
import scipy as sp

np.seterr(all="ignore")
from threeML import *
from threeML.io.package_data import get_path_of_data_file
from threeML.utils.data_download.Fermi_LAT.download_LAT_data import LAT_dataset

from astropy.io import fits as pyfits


from GtBurst import IRFS
irfs = IRFS.IRFS.keys()

sys.path.append('/home/pesce/work/Solar/solarflares/scripts/')
import MET2date,Date2MET
import SolarTemplates
import sunpy
#from sunpy.map import Map
from sunpy.coordinates import frames
import astropy
from astropy.coordinates import SkyCoord,get_sun
import math
import datetime
import ephem
def d(x):    return (math.degrees(float(repr(x))))

from threeML.utils import bayesian_blocks

from astropy.io import fits as pyfits
import pickle as pickle

from GtBurst import dataHandling, angularDistance

def str2datetime(datestr):
    """
    Convert a date string in the format 'yyyy-mm-ddTHH:MM:SS' to a datetime object.
    
    :param datestr: Date string in the format 'yyyy-mm-ddTHH:MM:SS'.
    :return: datetime object.
    """
    DateTime=datestr.split('T')
    ymd=DateTime[0].split('-')
    hms=DateTime[1].split(':')
    yyyy= int(ymd[0])
    mm  = int(ymd[1])
    dd  = int(ymd[2])
    hr  = int(hms[0])
    mi  = int(hms[1])
    se  = int(hms[2].split('.')[0])
    return datetime.datetime(yyyy,mm,dd,hr,mi,se)
   
def build_source_name(SolarFlareDate):
    SourceString = 'SOL%02d%02d%02d%02d'%(SolarFlareDate.year,SolarFlareDate.month,SolarFlareDate.day,SolarFlareDate.hour)
    return SourceString

def sun_in_field_of_view(SolarFlareDate, RA, DEC, flarename):
    """
    Check if the Sun is in the Fermi LAT field of view at the given time.
    
    :param SolarFlareDate: datetime object of the solar flare.
    :return: True if the Sun is in the field of view, False otherwise.
    """
    ra_sun  = RA
    dec_sun = DEC
    ft2filepath = 'FermiData/bn%s/gll_ft2_tr_bn%s_v00.fit' % (flarename,flarename)
    ft2_file = pyfits.open(ft2filepath)
    ra_scz = ft2_file['SC_DATA'].data.field("RA_SCZ")
    dec_scz = ft2_file['SC_DATA'].data.field("DEC_SCZ")
    ra_zenith = ft2_file['SC_DATA'].data.field("RA_ZENITH")
    dec_zenith = ft2_file['SC_DATA'].data.field("DEC_ZENITH")
    time = ft2_file['SC_DATA'].data.field("START")
    triggerTime = Date2MET.computeMET_datetime(SolarFlareDate)
    time = np.array([x - triggerTime for x in time])
    ft2_file.close()

    time_mask = time >= 0
    time = time[time_mask]
    ra_scz = ra_scz[time_mask]
    dec_scz = dec_scz[time_mask]
    theta = [angularDistance.getAngularDistance(x[0], x[1], ra_sun, dec_sun) for x in zip(ra_scz, dec_scz)]
    theta = np.array(theta)
    masked_theta = theta[theta < 70.0]
    # Check if the Sun is in the field of view at least for 3 time bins
    if len(masked_theta) > 2:
        print("Sun is in the field of view for %d time bins." % len(masked_theta))
        return True
    else:
        print("Sun is not in the field of view.")
        return False
    
    

def set_up_likelihood_analysis(T0, dt_minutes=10):
    """
    Set up the likelihood analysis for the solar flare.
    T0 is the trigger time from GOES or STIX. it will be in the form of yyyy-mm-ddTHH:MM:SS
    """
    print("Setting up likelihood analysis...")
    #grab the T0 in the form of  yyyy-mm-ddTHH:MM:SS and convert it to a datetime object  
    
    SolarFlareDate = str2datetime(T0)

    FinalTime      = SolarFlareDate + datetime.timedelta(minutes=dt_minutes)

    MET            = Date2MET.computeMET_datetime(SolarFlareDate)


    sun  = ephem.Sun()
    sun.compute(SolarFlareDate,epoch='2000')
    RA  = d(sun.a_ra)
    DEC = d(sun.a_dec)
    print ('Date   = ',SolarFlareDate,FinalTime)
    print ('MET    = ',MET)
    print ('RA,DEC =',d(sun.a_ra),d(sun.a_dec))

    roi       = 10
    zmax      = 100.
    thetamax  = 180.0
    emin      = 60.0
    emax      = 10000.0
    #irfs      = 'p8_source'
    irfs      = 'p8_transient020'
    strategy  = 'time'
    tstart    = 0
    tstop     = (int((FinalTime - SolarFlareDate).total_seconds()))
    flemin    = 100
    flemax    = 1000


    myLATdataset = LAT_dataset()

    myLATdataset.make_LAT_dataset(
        ra                    = RA,
        dec                   = DEC,
        radius                = 12,
        trigger_time          = MET,
        tstart                = tstart,
        tstop                 = tstop,
        data_type             = "Extended",
        destination_directory = 'FermiData',
        Emin= emin,    
        Emax= emax
        )

    print("Checking if the Sun is in the field of view...")
    sun_in_fov = sun_in_field_of_view(SolarFlareDate, RA, DEC, myLATdataset.grb_name)
    if not sun_in_fov:
        print("Sun is not in the field of view. Skipping analysis.")
        return None, RA, DEC, SolarFlareDate, FinalTime, myLATdataset.grb_name
        #exit()


    myLATdataset.extract_events(roi, zmax, irfs, thetamax, strategy=strategy)
   

    T0     = tstart
    T1     = tstop
    tstarts='%.3f' % tstart
    tstops ='%.3f' % tstop
    print("Time for time-integrated analysis: %s,%s" %(tstarts,tstops))
    tstarts = tstarts.replace('-','\\-')
    tstops  = tstops.replace('-','\\-')

    flarename = myLATdataset.grb_name

    analysis_builder = TransientLATDataBuilder(myLATdataset.grb_name,
                                           outfile=myLATdataset.grb_name,
                                           roi=roi,
                                           tstarts=tstarts,
                                           tstops=tstops,
                                           irf=irfs,
                                           zmax=zmax,
                                           source_model = 'PowerLaw',
                                           galactic_model='template',
                                           particle_model='isotr template',
                                           datarepository='FermiData',
                                           flemin=flemin,
                                           flemax=flemax,
                                           emin=emin,
                                           emax=emax,
                                           strategy=strategy,
                                           thetamax=thetamax)


    return analysis_builder, RA, DEC, SolarFlareDate, FinalTime, flarename
    


def doAnalysis(T0):   
    """
    Run a simpl likelihood analysis over a given time inteval using the simple power-law model.
    
    :param T0: T0 from GOES or STIX
    :param T1: T0 + duration fixed by user
    :param RA: Right Ascension of the Sun at T0
    :param DEC: Declination of the Sun at T0
    :param analysis_builder: The TransientLATDataBuilder object
    :return: JointLikelihood object after fitting.
    """
    #This is needed so the script doesn't fail when looking for the previeous intervals
    os.system('rm -rf interval*')

    analysis_builder, RA, DEC, SolarFlareDate, FinalTime, flarename = set_up_likelihood_analysis(T0)
    if analysis_builder is None:
        return -1.0

    print("----------doAnalysis from %s to %s using a power-law ------" %(SolarFlareDate,FinalTime))
    
    df=analysis_builder.display(get=True)
    LAT_observations = analysis_builder.run(include_previous_intervals = True)

    source_name = build_source_name(SolarFlareDate)
    
    LAT_plugins={}     
    
    LAT_name = 'LAT_%s' % (source_name)
    t0 = 0 
    t1 = int((FinalTime - SolarFlareDate).total_seconds())
    LAT_model_name = ('LAT%dX%d' % (t0,t1)).replace('-','n')
    

    for l in LAT_observations:
        LAT_plugins[LAT_name] = l.to_LATLike()

    datalist = DataList(LAT_plugins[LAT_name])
    

    sun = PointSource('sun',ra=RA,dec=DEC,spectral_shape=Powerlaw_flux())

    model = Model(sun)

    model.sun.spectrum.main.Powerlaw_flux.F       =  1e-5
    model.sun.spectrum.main.Powerlaw_flux.index   =  -2.1
    model.sun.spectrum.main.Powerlaw_flux.a       =  1e5
    model.sun.spectrum.main.Powerlaw_flux.b       =  1e8


    model.display(complete=True)
   
    jl=JointLikelihood(model,datalist,verbose=False)
    model[LAT_model_name+'_GalacticTemplate_Value'].value=1.0
    model[LAT_model_name+'_GalacticTemplate_Value'].fix=True
    model[LAT_model_name+'_GalacticTemplate_Value'].fix=True
    jl.set_minimizer('MINUIT')

   
    dfs = jl.fit(compute_covariance=True)
    
    display_spectrum_model_counts(jl, step=False,figsize=(10,10), ratio_residuals=True)
    TS = {}
    TS[LAT_name] = jl.compute_TS('sun',dfs[1])
    ts  = TS[LAT_name]['TS'][0]
    if ts >= 20:
        print("Significant detection with TS=%.2f" % ts)
    else:
        print("No significant detection. TS=%.2f" % ts)

    return ts





if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--goes_t0", type=str, help="GOES flare trigger time (MET)")
    args = parser.parse_args()

    t0 = args.goes_t0 

    ts = doAnalysis(t0)
    outstr = "likelihood_results.txt" 
    with open(outstr,'w') as f:
        f.write("#GOES_T0\t TS\n")
        f.write("%s,%.2f\n" % (t0, ts))
    print("Results written to %s" % outstr)
    print("Done.")