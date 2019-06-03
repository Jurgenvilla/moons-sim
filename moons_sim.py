from PyAstronomy import pyasl
import numpy as np
from astropy.io import fits
import sys
from scipy import ndimage
import matplotlib
import matplotlib.pyplot as plt
import scipy.io as sc
import math
import warnings

def get_input():
    uservals={}
    uservals['telescope']='VLT'
    uservals['instrument']='MOONS'
    uservals['template']=str(sys.argv[1])
    uservals['outfits']=str(sys.argv[2])
    uservals['moons_mode']=str(sys.argv[3])
    uservals['band']=str(sys.argv[4])
    uservals['ab']=float(sys.argv[5])
    uservals['dit']=float(sys.argv[6])
    uservals['N_dit']=float(sys.argv[7])
    uservals['seeing']=float(sys.argv[8])
    uservals['airmass']=str(sys.argv[9])
    uservals['airmass_fl']=float(sys.argv[9])
    uservals['atm_dif_ref']=float(sys.argv[10])
    uservals['sky_residual']=float(sys.argv[11])
    uservals['sky_template']=sys.argv[12]
    uservals['telluric']=int(sys.argv[13])
    uservals['set_line_profile']='NO' # do not change...to be applied later
    # detailed check for valid inputs to be added here
    return uservals

def setup_moons(uservals):
    instrumentconfig={}
    if (uservals['moons_mode'] == "HR"):
        if ( uservals['band'] == "H" ) :
            print("Adopting moons mode: ", uservals['moons_mode'], uservals['band'])
            instrumentconfig['wlr']=[1.521,1.641]
            instrumentconfig['RON'] = 3.
            DK = 30.
            instrumentconfig['gain']=0.5 #ADU/e- for 4RGs consistent with Gianluca's model
            instrumentconfig['saturation_level']=30000.0 #ADU
            instrumentconfig['pix_size']=15.0
            instrumentconfig['spec_sampling']=2.7
            instrumentconfig['resolution']=19400.00
            print("Setting wavelength range to ", str(instrumentconfig['wlr']))
            instrumentconfig['template_wl_norm']=1.630
            QE_file='Inst_setup/QE_4RG.txt'
            instrumentconfig['QE_wav'], instrumentconfig['QE_eff']  = np.loadtxt(QE_file, unpack = True)
        if ( uservals['band'] == "YJ" ) :
            print("Adopting moons mode: ", uservals['moons_mode'], uservals['band'])
            instrumentconfig['wlr']=[0.934,1.350]
            instrumentconfig['RON'] = 3.
            DK = 30.
            instrumentconfig['gain']=0.5 #ADU/e- for 4RGs consistent with Gianluca's model
            instrumentconfig['saturation_level']=30000.0 #ADU
            instrumentconfig['pix_size']=15.0
            instrumentconfig['spec_sampling']=2.7
            instrumentconfig['resolution']=4000.00
            print("Setting wavelength range to ", str(instrumentconfig['wlr']))
            instrumentconfig['template_wl_norm']=1.22
            QE_file='Inst_setup/QE_4RG.txt'
            instrumentconfig['QE_wav'], instrumentconfig['QE_eff']  = np.loadtxt(QE_file, unpack = True)
        if ( uservals['band'] == "RI" ) :
            print("Adopting moons mode: ", uservals['moons_mode'], uservals['band'])
            instrumentconfig['wlr']=[0.765,0.898]
            DK = 2.
            instrumentconfig['gain']=0.6 #ADU/e- for LBNL consistent with Gianluca's model
            instrumentconfig['saturation_level']=58000.0 #ADU
            instrumentconfig['pix_size']=15.0
            instrumentconfig['spec_sampling']=2.7
            instrumentconfig['resolution']=9400.00
            print("Setting wavelength range to ", str(instrumentconfig['wlr']))
            instrumentconfig['template_wl_norm']=0.797
            #setting detector_model_RI to 'LBNL':
            QE_file='Inst_setup/QE_LBNL.txt'
            instrumentconfig['QE_wav'], instrumentconfig['QE_eff']  = np.loadtxt(QE_file, unpack = True)
            instrumentconfig['RON']=3.0
    if (uservals['moons_mode'] == "LR"):
        if ( uservals['band'] == "H" ) :
            print("Adopting moons mode: ", uservals['moons_mode'], uservals['band'])
            instrumentconfig['wlr']=[1.452,1.800]
            instrumentconfig['RON'] = 3.
            DK = 30.
            instrumentconfig['gain']=0.5 #ADU/e- for 4RGs consistent with Gianluca's model
            instrumentconfig['saturation_level']=30000.0 #ADU
            instrumentconfig['pix_size']=15.0
            instrumentconfig['spec_sampling']=2.7
            instrumentconfig['resolution']=6400.00
            print("Setting wavelength range to ", str(instrumentconfig['wlr']))
            instrumentconfig['template_wl_norm']=1.63
            QE_file='Inst_setup/QE_4RG.txt'
            instrumentconfig['QE_wav'], instrumentconfig['QE_eff']  = np.loadtxt(QE_file, unpack = True)
        if ( uservals['band'] == "YJ" ) :
            print("Adopting moons mode: ", uservals['moons_mode'], uservals['band'])
            instrumentconfig['wlr']=[0.934,1.350]
            instrumentconfig['RON'] = 3.
            DK = 30.
            instrumentconfig['gain']=0.5 #ADU/e- for 4RGs consistent with Gianluca's model
            instrumentconfig['saturation_level']=30000.0 #ADU
            instrumentconfig['pix_size']=15.0
            instrumentconfig['spec_sampling']=2.7
            instrumentconfig['resolution']=4500.00
            print("Setting wavelength range to ", str(instrumentconfig['wlr']))
            instrumentconfig['template_wl_norm']=1.22
            QE_file='Inst_setup/QE_4RG.txt'
            instrumentconfig['QE_wav'], instrumentconfig['QE_eff']  = np.loadtxt(QE_file, unpack = True)
        if ( uservals['band'] == "RI" ) :
            print("Adopting moons mode: ", uservals['moons_mode'], uservals['band'])
            instrumentconfig['wlr']=[0.647,0.934]
            instrumentconfig['pix_size']=15.0
            instrumentconfig['spec_sampling']=2.7
            instrumentconfig['resolution']=5400.00
            DK = 2.
            instrumentconfig['RON']=3.
            instrumentconfig['saturation_level']=58000.0 #ADU
            instrumentconfig['gain']=0.6 #ADU/e- for LBNL consistent with Gianluca's model
            print("Setting wavelength range to ", str(instrumentconfig['wlr']))
            instrumentconfig['template_wl_norm']=0.797
            #setting detector_model_RI to 'LBNL':
            QE_file='Inst_setup/QE_LBNL.txt'
            instrumentconfig['QE_wav'], instrumentconfig['QE_eff']  = np.loadtxt(QE_file, unpack = True)
    instrumentconfig['DKsec'] = DK/3600.
    instrumentconfig['sky_aperture']=1.1
    print('')
    return instrumentconfig

def set_telescope(uservals):
    telescopeconfig={}
    if (uservals['telescope']=='ELT'):
        telescopeconfig['t_aperture']=39.0
    if (uservals['telescope']=='VLT'):
        telescopeconfig['t_aperture']=8.1
    return telescopeconfig

def set_detector(telescopeconfig,uservals,instrumentconfig):
    #equivalent pixel sizes in arcsec
    detectorconfig={}
    xanpix=0.072*(instrumentconfig['pix_size']/15.0)*(1.1/1.0)*(39.0/telescopeconfig['t_aperture'])
    yanpix=0.072*(instrumentconfig['pix_size']/15.0)*(1.1/1.0)*(39.0/telescopeconfig['t_aperture'])
    detectorconfig['ypix_fwhm']=instrumentconfig['sky_aperture']/yanpix #math.pi/4.0*sky_aperture*sky_aperture/(xanpix*yanpix)
    print(' ')
    print('Spectral sampling for current configuration: ',str(round(instrumentconfig['spec_sampling'],1)))
    print(' ')
    print('Spectral resolving power for current configuration: ',str(round(instrumentconfig['resolution'],1)))
    print(' ')
    detectorconfig['disp']=(instrumentconfig['wlr'][1]+instrumentconfig['wlr'][0])/2.0*1.0e4/instrumentconfig['resolution']/instrumentconfig['spec_sampling']
    detectorconfig['npix']=int(1.0e4*(instrumentconfig['wlr'][1]-instrumentconfig['wlr'][0])/detectorconfig['disp']) # Total number of pixels by fixing spectral range and sampling requirements in current baseline design
    print('Spectral dispersion for current configuration: ',str(round(detectorconfig['disp'],1)))
    print(' ')
    return detectorconfig

def get_template(uservals):
    template_name=uservals['template']
    try:
        print(' ')
        print("Reading FITS template spectrum: %s " %template_name)
        fits.getdata(str(template_name))
    except FileNotFoundError:
        print("FITS file not found or not valid input file")
        exit()
    hdu=fits.open(str(template_name))
    spec=hdu[0].data
    spec_header=hdu[0].header
    naxis=spec_header['naxis1']
    crval=spec_header['crval1']
    cdelt=spec_header['cdelt1']
    tunit_wav=spec_header['TUNIT1']
    template_data={}
    template_data['header']=spec_header
    temp_pix_array=np.arange(0,naxis,1)
    # This are the wave and flux arrays of the templates
    unitsok=0
    if (tunit_wav=='Angstroms'):
        wl_tpl=crval+temp_pix_array*cdelt
        fl_tpl=spec
        template_data['wave']=wl_tpl # in Angstroms
        template_data['cdelt']=cdelt
        template_data['flux']=fl_tpl # in ergs/s/cm2/A
        unitsok=1
    if (tunit_wav=='nm'): # specific for templates from Jorge. To be removed later on....
        z_val=1.5
        wl_tpl=(crval+temp_pix_array*cdelt)*10.0 #in A
        fl_tpl=spec
        template_data['wave']=wl_tpl*(1.+float(z_val)) #in A
        template_data['cdelt']=cdelt*10.0 # in A
        template_data['flux']=fl_tpl # in ergs/s/cm2/A
        unitsok=1
    if (unitsok==0):
        print('ERROR: Wavelength unit expect to be Angstroms or nm')
        exit()
    return template_data

def get_diffraction(Lambda,airmass,atm_ref_wav):
    Lambda1=0.5
    Lambda2=1.9
    Lambda_Step=0.10
    Lambda0=float(atm_ref_wav)
    TC=11.5			#Temperature [C]
    RH=14.5			#Relative Humidity [%]
    P=743.0			#Pressure [mbar]
    Z=np.arccos(1.0/float(airmass))*57.2958
    ZD=Z*0.0174533
    T=TC+273.16
    PS=-10474.0+116.43*T-0.43284*T**2+0.00053840*T**3
    P2=RH/100.0*PS
    P1=P-P2
    D1=P1/T*(1.0+P1*(57.90*1.0e-8-(9.3250*1.0e-4/T)+(0.25844/T**2)))
    D2=P2/T*(1.0+P2*(1.0+3.7e-4*P2)*(-2.37321e-3+(2.23366/T)-(710.792/T**2)+(7.75141e4/T**3)))
    S0=1.0/Lambda0
    S=1.0/Lambda
    N0_1=1.0E-8*((2371.34+683939.7/(130.0-S0**2)+4547.3/(38.9-S0**2))*D1+(6487.31+58.058*S0**2-0.71150*S0**4+0.08851*S0**6)*D2)
    N_1=1.0E-8*((2371.34+683939.7/(130.0-S**2)+4547.3/(38.9-S**2))*D1+(6487.31+58.058*S**2-0.71150*S**4+0.08851*S**6)*D2)
    DR=np.tan(ZD)*(N0_1-N_1)*206264.8
    return DR

def create_fits(sim_spectrum, uservals,template_data,instrumentconfig): #sim_spectrum, outfits,template_data,gain):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        hdu = fits.PrimaryHDU()
        hdu.data=sim_spectrum['flux']*instrumentconfig['gain']
        hdu.header['CDELT1']=sim_spectrum['cdelt']
        hdu.header['CRVAL1']=sim_spectrum['crval']
        hdu.header['TUNIT1']='Angstroms'
        hdu.header['MTYPE']=template_data['header']['MTYPE']
        hdu.header['MNAME']=template_data['header']['MNAME']
        hdu.header['TUNIT2']='ADU'
        hdu.header['GAIN']=instrumentconfig['gain']
        hdu.header.comments['GAIN']='ADU/e-'
        hdu.writeto(uservals['outfits'],clobber=True)

def get_line_profile(filename, wave_spec, flux_spec, disp):
    profile_angs, profile_line  = np.loadtxt(filename, unpack=True)
    k_x=np.arange(0,len(profile_angs),1)
    k_x_wlr=profile_angs[0]+k_x*(profile_angs.max()-profile_angs.min())/len(k_x)
    k_x_sim = np.interp(k_x_wlr,profile_angs,profile_line)
    flux_conv=np.convolve(flux_spec,k_x_sim,mode="same") #Convolve with line profile kernel, resampled to simulation dispersion steps
    return flux_conv

def get_sky_model(sky_template,airmass_fl):
    if sky_template=='eso_skycalc':
        print('Using ESO SkyCalc template for Mean Sky conditions')
        available_airmass=np.array([1.0,1.2,1.4,1.6,1.8,2.0])
        closest_airmass=available_airmass[np.argmin(np.abs(available_airmass-airmass_fl))]
        print('Skymodel selected: SkyTemplate_ESO_a'+np.str(closest_airmass)+'.fits')
        spec_hdu=fits.open('Skymodel/SkyTemplate_ESO_a'+np.str(closest_airmass)+'.fits')
        spec=spec_hdu[0].data
        trans_spec=spec_hdu[1].data
        spec_header=spec_hdu[0].header
        trans_header=spec_hdu[1].header
        naxis_sf=spec_header['naxis1']
        crval_sf=spec_header['crval1']
        cdelt_sf=spec_header['cdelt1']
        temp_pix_array_sf=np.arange(0,naxis_sf,1)
        naxis_st=spec_header['naxis1']
        crval_st=spec_header['crval1']
        cdelt_st=spec_header['cdelt1']
        temp_pix_array_st=np.arange(0,naxis_st,1)
        # This are the wave and flux arrays of the sky templates
        rwl0=crval_sf+temp_pix_array_sf*cdelt_sf
        rfn0=spec
        atmwl=crval_st+temp_pix_array_st*cdelt_st
        atmtr=trans_spec
        oh_f={}
        oh_f['rwl0']=rwl0
        oh_f['rfn0']=rfn0
        oh_f['atmwl']=atmwl
        oh_f['atmtr']=atmtr
        oh_f['rdwl']=cdelt_sf
    else:
        skyfile=sky_template
        try:
            fits.getdata(str(skyfile))
        except FileNotFoundError:
            print("FITS file not found or not valid input file")
            exit()
        spec_hdu=fits.open(str(skyfile))
        spec=spec_hdu[1].data
        trans_spec=spec_hdu[2].data
        spec_header=spec_hdu[1].header
        trans_header=spec_hdu[2].header
        naxis_sf=spec_header['naxis1']
        crval_sf=spec_header['crval1']
        cdelt_sf=spec_header['cdelt1']
        temp_pix_array_sf=np.arange(0,naxis_sf,1)
        naxis_st=spec_header['naxis1']
        crval_st=spec_header['crval1']
        cdelt_st=spec_header['cdelt1']
        temp_pix_array_st=np.arange(0,naxis_st,1)
        # This are the wave and flux arrays of the sky templates
        rwl0=crval_sf+temp_pix_array_sf*cdelt_sf
        rfn0=spec
        atmwl=crval_st+temp_pix_array_st*cdelt_st
        atmtr=trans_spec
        oh_f={}
        oh_f['rwl0']=rwl0
        oh_f['rfn0']=rfn0
        oh_f['atmwl']=atmwl
        oh_f['atmtr']=atmtr
        oh_f['rdwl']=cdelt_sf
    return oh_f

def get_efficiency(outputwl,min_wl,max_wl):
    # Add Telescope transmission
    if (uservals['telescope']=='VLT'):
        telescope_file='Inst_setup/telescope_eff.txt'

    telescope_wav_micron, telescope_eff  = np.loadtxt(telescope_file, unpack = True)
    telescope_wav=telescope_wav_micron*1.0e04
    # Resample telescope transmission on wlgrid
    tel_eff = np.interp(outputwl,telescope_wav,telescope_eff)
    # Add transmission curve from overall instrument+telescope efficiency
    #data_dir='data_dir/Inst_setup'
    if (uservals['instrument']=='MOONS'):
        trans_file='Inst_setup/throughput.sav'
        efficiency=sc.readsav(trans_file,verbose=False)
        eff_lr0=efficiency.lreff
        eff_hr0=efficiency.hreff
        eff_wl0=efficiency.weff*10.0
        effok = np.where((eff_wl0 > min_wl) & (eff_wl0 < max_wl))[0] #Selecting instrument transmission in Setup
        eff_lr=eff_lr0[effok]
        eff_hr=eff_hr0[effok]
        eff_wl=eff_wl0[effok]
        # Resample instrument efficiency on wlgrid and detector QE (updated to detector FDR)
        if (uservals['moons_mode'] == "HR"):
            eff = np.interp(outputwl,eff_wl,eff_hr)
            QE_detector=np.interp(outputwl,instrumentconfig['QE_wav']*1.e4,instrumentconfig['QE_eff'])
            eff = QE_detector/100.0*eff*tel_eff
        if (uservals['moons_mode'] == "LR"):
            eff = np.interp(outputwl,eff_wl,eff_lr)
            QE_detector=np.interp(outputwl,instrumentconfig['QE_wav']*1.e4,instrumentconfig['QE_eff'])
            eff = QE_detector/100.0*eff*tel_eff
    print('Mean total efficiency (Telescope+Instrument+Detector): ',str(round(np.mean(eff),2)))
    print(' ')
    return eff

def seeing_to_ImageQ(seeing,cen_wav,airmass_fl):
    r_0=0.100/seeing*(cen_wav/500.0)**(1.2)*airmass_fl**(-0.6)
    F_kolb=-0.981644
    fwhm_atm=seeing*airmass_fl**(0.6)*(cen_wav/500.0)**(-0.2)*np.sqrt(1.0+F_kolb*2.183*(r_0/46.0)**(0.356))
    fwhm_tel=0.000212*(cen_wav/8.1)
    fwhm_iq=np.sqrt(fwhm_tel**2+fwhm_atm**2)
    print(" ")
    print('Expected Image Quality in selected band:', round(fwhm_iq,2))
    print(" ")
    return fwhm_iq

def get_saturation(count_level,saturation_level):
    check_level=np.where(count_level>saturation_level)
    if (np.size(check_level)>0):
        return True
    else:
        return False
def get_PeakIntensity(spec,npix_y,DK_perpix):
    spec2d_central=np.zeros([np.int(npix_y)*2,np.size(spec)])
    central_pix=np.int(npix_y-1)
    spec2d_central[central_pix]=spec
    spec2d=ndimage.gaussian_filter1d(spec2d_central,npix_y/2.0/2.355,axis=0)
    peak_intensity=spec2d[central_pix]+DK_perpix
    return peak_intensity

def make_simulation(template_data, uservals, detectorconfig, telescopeconfig,instrumentconfig):#template_data,resolution,wlr,t_aperture,sky_aperture,template_wl_norm,ab,airmass,npix,Instrument,aper_sampl,dit,eff_opt,seeing,atm_dif_ref):

    ####################################################################
    #Set dispersion axis
    cen_wav=(instrumentconfig['wlr'][1]+instrumentconfig['wlr'][0])/2.0*1.0e3
    wav_range_length=(instrumentconfig['wlr'][1]-instrumentconfig['wlr'][0])*1.0e3
    pix_arr=np.arange(0,detectorconfig['npix'],1)
    outputwl = instrumentconfig['wlr'][0]*1.0e4+pix_arr*detectorconfig['disp']
    ####################################################################

    ####################################################################
    #obtain Sky spectrum from template (either provided or ESO sky_calc)
    oh_f=get_sky_model(uservals['sky_template'], uservals['airmass_fl'])
    rdwl=oh_f['rdwl']
    rwl0=oh_f['rwl0'] # wavelength of sky model in Angstroms
    rfn0=oh_f['rfn0']*math.pi*(telescopeconfig['t_aperture']*1.e2/2.0)**2*rdwl #Photons/s/arcsec2 per pixel of the OH grid
    atmtr0=oh_f['atmtr'] # atmospheric transmission
    atmwl0=oh_f['atmwl']
    tol = (instrumentconfig['wlr'][1]-instrumentconfig['wlr'][0])*0.05
    min_wl=(instrumentconfig['wlr'][0]-tol)*1.0e4
    max_wl=(instrumentconfig['wlr'][1]+tol)*1.0e4
    iok_sky = np.where((rwl0 > min_wl) & (rwl0 < max_wl))[0]
    iok_atm = np.where((atmwl0 > min_wl) & (atmwl0 < max_wl))[0]
    rfn = rfn0[iok_sky] # flux of OH spectrum
    rwl = rwl0[iok_sky] # wave of OH spectrum
    atmwl=atmwl0[iok_atm]
    atmtr=atmtr0[iok_atm]
    rfn_sky = (rfn)*math.pi*(instrumentconfig['sky_aperture']/2.0)**2 # Photons/s/pix of the OH spectrum ALONE

    ####################################################################

    ####################################################################
    #Prepare the input template
    fl_tpl_in=template_data['flux']
    wl_tpl_in=template_data['wave']
    setup_range_tpl=np.where((wl_tpl_in>min_wl) & (wl_tpl_in<max_wl))[0]
    fl_tpl=fl_tpl_in[setup_range_tpl]
    wl_tpl=wl_tpl_in[setup_range_tpl]
    # normalising the template source spectrum to the magnitude
    i = np.where(wl_tpl >= instrumentconfig['template_wl_norm']*1.0e4)[0]
    templ_fl_ref = fl_tpl[i[0]]
    fl_tpl_normalized_to_one = fl_tpl/templ_fl_ref
    # template photon flux per A in cgs assuming a telescope effective diameter
    templ_phot_fl = fl_tpl_normalized_to_one*3.63e3*10.0**(-(uservals['ab']+57.5)/2.5)*math.pi*(telescopeconfig['t_aperture']*1.0e2/2.0)**2/(wl_tpl*6.626e-27)
    # template photon flux per PER SPECTRAL PIXEL OF THE TEMPLATE
    templ_phot_fl_pix = templ_phot_fl*template_data['cdelt']#*rdwl
    npix_y=detectorconfig['ypix_fwhm']*2.0
    ####################################################################

    ####################################################################
    # Modelling Fibre Injection
    #Obtain Image Quality in corresponding band from seeing provided (seeing defined in zenith at 500nm)
    fwhm_iq=seeing_to_ImageQ(uservals['seeing'],cen_wav,uservals['airmass_fl'])
    uservals['seeing']=fwhm_iq
    fib_frac=1.0 # allowing for tuning of fraction loss manually, keep as 1.0 for no additional loss
    #Calculate atmospheric difraction effect
    Lambda=np.arange(0.5,1.9,0.1)
    atm_diff=get_diffraction(Lambda, uservals['airmass_fl'], uservals['atm_dif_ref'])
    atm_diff_oh=np.interp(wl_tpl,Lambda*1.0e4,atm_diff)
    seeing_arr=np.sqrt(uservals['seeing']**2+atm_diff_oh**2) #Co-added effect of seeing and atm_diff offset.
    seeing_cor_flux= (1.0-1.0/np.exp( (np.sqrt(math.log(2.0)) * instrumentconfig['sky_aperture']/seeing_arr)**2))*templ_phot_fl_pix
    sp_conv_src = fib_frac*seeing_cor_flux
    sp_conv_sky=rfn_sky
    ####################################################################
    #resampling to detector pixels, conserving total flux
    #transmission:
    atminterp_res, fwhm = pyasl.instrBroadGaussFast(atmwl, atmtr, instrumentconfig['resolution'],edgeHandling="firstlast", fullout=True)
    atminterp = np.interp(outputwl,atmwl,atminterp_res)# resample atmospheric transmission to detector pixels
    #
    #source:
    sp_conv_src_res, fwhm = pyasl.instrBroadGaussFast(wl_tpl, sp_conv_src, instrumentconfig['resolution'],edgeHandling="firstlast", fullout=True)
    sp_det_src = np.interp(outputwl,wl_tpl,sp_conv_src_res)
    inband = np.where((wl_tpl >= outputwl[0]) & (wl_tpl <= outputwl[detectorconfig['npix']-1]))[0]
    sp_conv_src_sum=sp_conv_src_res[inband].sum()
    sp_det_src_sum=sp_det_src.sum()
    renorm=sp_conv_src_sum/sp_det_src_sum
    sp_det_src_rn = sp_det_src*renorm
    if (uservals['set_line_profile']=='YES'):
        line_profile_file='IP_HIRES_m68_wave17834.1A.txt'
        print(' ')
        print('Applying line profile convolution')
        print('Extracting line profile from file: ',line_profile_file)
        sp_det_src_rn=get_line_profile(line_profile_file,outputwl,sp_det_src_rn,detectorconfig['disp'])
    else:
        print('No LSF provided, adopting Gaussian kernel convolution')
        print(' ')
    #
    #sky emission:
    sp_conv_sky_res, fwhm = pyasl.instrBroadGaussFast(rwl, sp_conv_sky, instrumentconfig['resolution'],edgeHandling="firstlast", fullout=True)
    sp_det_sky = np.interp(outputwl,rwl,sp_conv_sky_res)
    inband = np.where((rwl >= outputwl[0]) & (rwl <= outputwl[detectorconfig['npix']-1]))[0]
    sp_conv_sky_sum=sp_conv_sky_res[inband].sum()
    sp_det_sky_sum=sp_det_sky.sum()
    renorm_sky=sp_conv_sky_sum/sp_det_sky_sum
    sp_det_sky_rn = sp_det_sky*renorm_sky # for sky only

    if (uservals['set_line_profile']=='YES'):
        line_profile_file='IP_HIRES_m68_wave17834.1A.txt'
        sp_det_sky_rn=get_line_profile(line_profile_file,outputwl,sp_det_sky_rn,detectorconfig['disp'])
    ####################################################################

    ####################################################################
    #Modelling overall efficiency
    eff=get_efficiency(outputwl,min_wl,max_wl)
    #Calculate telescope and instrument emissivity (not enabled as is not critical for H band)
    ThBK=283.00
    EBK=0.00 #36.0*selector #selector matching temperature of the telescope
    ThBK_ins=283.00
    EBK_ins=0.10 #36.0*selector #selector matching temperature of the telescope
    t_em=1.4*10.0**12*EBK*np.exp(-14388.0/(outputwl/1.0e4*ThBK))/((outputwl/1.0e4)**3/instrumentconfig['resolution'])*detectorconfig['disp']
    ins_em=1.4*10.0**12*EBK_ins*np.exp(-14388.0/(outputwl/1.0e4*ThBK_ins))/((outputwl/1.0e4)**3/instrumentconfig['resolution'])*detectorconfig['disp']
    NBK_tel=(t_em)*math.pi*(instrumentconfig['sky_aperture']/2.0)**2
    NBK_ins=(ins_em)*math.pi*(instrumentconfig['sky_aperture']/2.0)**2
    NBK=NBK_tel+NBK_ins
    #print("Thermal background emissivity Telescope [e-/s]: ",str(round(np.max(NBK_tel*eff),3)))
    #print("Thermal background emissivity Instrument [e-/s]: ",str(round(np.max(NBK_ins*eff),3)))
    #sp_det_sky_rn=sp_det_sky_rn+NBK commented out for the moment.
    ####################################################################
    #Modelling detector influx
    sp_det = sp_det_src_rn+sp_det_sky_rn # total detector flux from both sky and source
    sp_eff=sp_det*eff*atminterp
    sp_eff_src = sp_det_src_rn*eff*atminterp
    sp_eff_sky = sp_det_sky_rn*eff*atminterp
    # number of pixels COLLAPSED ALONG Y DIRECTION for MOS
    #sp_dk = sp_eff+DKsec*npix_y # from total insident flux
#    sp_dk_sky = sp_eff_sky+DKsec*npix_y #from only sky
    # scale by integration time
    spec_total = sp_eff*uservals['dit']
    spec_2d_peak_intensity=get_PeakIntensity(spec_total,npix_y,instrumentconfig['DKsec']*uservals['dit'])
    # the following is for the sky only
    spec_sky = sp_eff_sky*uservals['dit']
    # the following is for the source alone
    spec_source = sp_eff_src*uservals['dit']
    #calculate total noise
    #####################################################
    # Add stray light contribution
    stray=1.0 # 1% difusse stray light contribution as per latest optical modelling
    total_stray=spec_sky.sum()*uservals['N_dit']*500.0/detectorconfig['npix']**2
    spec_stray=total_stray*float(stray)/100.0*npix_y
    #detector NOISE:
    noisedet=np.sqrt(npix_y*(uservals['N_dit']*instrumentconfig['RON']**2+instrumentconfig['DKsec']*uservals['dit']*uservals['N_dit']))
    #background NOISE (including stray):
    noiseback=np.sqrt(spec_sky*uservals['N_dit']+spec_stray)
    #Add residual sky subtraction (optional, if you want then change skyres to percentage):
    if ((uservals['sky_residual'] >= 0) & (uservals['sky_residual'] <= 100)):
        skyres=uservals['sky_residual']
    else:
        if (uservals['sky_residual'] == -1):
            print('Simulation with sky-subtraction OFF')
            print(' ')
            skyres=0
        else:
            print(' Not a valid sky residual value (0-100). Adopting 0.00 percent')
            print(' ')
            skyres=0
    noiseskyres=float(skyres)/100.0*spec_sky
    #total NOISE: comment out when necessary
    noisetot=np.sqrt(noiseback**2+noisedet**2+spec_source*uservals['N_dit'])
    #noise_sim_withsky=np.sqrt(noiseback**2+noisedet**2+spec_source*uservals['N_dit'])
    outnoise=noisetot+noiseskyres
    #outnoise = np.sqrt(outspec+outspec_sky+RON**2*npix_y*2.)
    sn_cont_all = spec_source/outnoise*uservals['N_dit']
    sn_cont = np.median(sn_cont_all)
    cent_range=np.where((outputwl >= (cen_wav-wav_range_length*0.1)*10.0) & (outputwl <= (cen_wav+wav_range_length*0.1)*10.0))[0]
    sn_central=sn_cont_all[cent_range].max()
    print("**** S/N at central wavelength = %.2f ****"%sn_central)
    print(" ")
    count_level=spec_2d_peak_intensity*instrumentconfig['gain']
    saturated=get_saturation(count_level,instrumentconfig['saturation_level'])
    if saturated:
        print('WARNING!!! Counts above saturation/linearity regime!!!')
    #Get figure with summary of results
    #SNR figure
    f=plt.figure(figsize=(10,8),dpi=100)
    ax1=f.add_subplot(221)
    ax1.plot(outputwl/1.0e4, sn_cont_all,label='SNR per pixel')
    ax1.axis([instrumentconfig['wlr'][0], instrumentconfig['wlr'][1], sn_cont_all.min(), sn_cont_all.max()])
    ax1.set_xlabel('Wavelength [um]')
    ax1.set_ylabel('SNR (/pix)')
    plt.legend(loc='upper right',prop={'size':8}, numpoints=1)
    #Sky spectrum
    ax1=f.add_subplot(222)
    ax1.plot(outputwl/1.0e4, spec_2d_peak_intensity,label='Peak Intensity')
    ax1.axis([instrumentconfig['wlr'][0], instrumentconfig['wlr'][1], spec_2d_peak_intensity.min(), spec_2d_peak_intensity.max()])
    if saturated:
        ax1.plot(outputwl/1.0e4,spec_2d_peak_intensity*0.0+instrumentconfig['saturation_level']/instrumentconfig['gain'],color='red',label='Saturation')
    plt.legend(loc='upper right',prop={'size':8}, numpoints=1)
    ax1.set_xlabel('Wavelength [um]')
    ax1.set_ylabel('Counts (e-)')
    #atmospheric transmission
    ax1=f.add_subplot(223)
    ax1.plot(outputwl/1.0e4, atminterp,label='Atmospheric transmission')
    ax1.axis([instrumentconfig['wlr'][0], instrumentconfig['wlr'][1], atminterp.min(), atminterp.max()])
    ax1.set_xlabel('Wavelength [um]')
    ax1.set_ylabel('Transmission fraction')
    plt.legend(loc='lower left',prop={'size':8}, numpoints=1)
    #Simulated spectrum
    ax1=f.add_subplot(224)
    normnoise=np.random.random((detectorconfig['npix']))*2.0-1.0
    res_noise=normnoise*outnoise
    nores_noise=normnoise*noisetot
    sim_spectrum={}
    if uservals['sky_residual']==-1:
        sim_spectrum['flux']=spec_total*uservals['N_dit']+nores_noise
    else:
        sim_spectrum['flux']=spec_source*uservals['N_dit']+res_noise
    if uservals['telluric']==1:
        print('Telluric correction applied (idealistic)')
        print(' ')
        sim_spectrum['flux']=sim_spectrum['flux']/atminterp
    else:
        sim_spectrum['flux']=sim_spectrum['flux']
    sim_spectrum['wave']=outputwl
    sim_spectrum['cdelt']=detectorconfig['disp']
    sim_spectrum['crval']=outputwl[0]
    ax1.plot(outputwl/1.0e4, sim_spectrum['flux'],label='Sim spectrum')
    ax1.plot(outputwl/1.0e4, spec_source*uservals['N_dit'],label='Object (no noise)',alpha=0.6)
    ax1.axis([instrumentconfig['wlr'][0], instrumentconfig['wlr'][1], sim_spectrum['flux'].min(), sim_spectrum['flux'].max()])
    ax1.set_xlabel('Wavelength [um]')
    ax1.set_ylabel('Counts (e-)')
    plt.legend(loc='upper right',prop={'size':8}, numpoints=1)
    plt.tight_layout()
    f.savefig("SIM_results.pdf", bbox_inches='tight')
    plt.close()
    return sim_spectrum

if __name__ == "__main__":
    uservals=get_input()
    if uservals['instrument']=='MOONS':
        instrumentconfig=setup_moons(uservals)
    telescopeconfig=set_telescope(uservals)
    detectorconfig=set_detector(telescopeconfig,uservals,instrumentconfig)
    template_data=get_template(uservals)
    sim_spectrum=make_simulation(template_data, uservals, detectorconfig, telescopeconfig, instrumentconfig) # resolution,wlr,t_aperture,sky_aperture,template_wl_norm,ab,airmass,npix,Instrument,ypix_fwhm,dit,eff_opt,seeing,atm_dif_ref)
    #cdelt1=detectorconfig['dist']#get_disp(resolution,spec_sampling,wlr)
    create_fits(sim_spectrum, uservals,template_data,instrumentconfig)
    print('Output fits file created: %s'%uservals['outfits'])
