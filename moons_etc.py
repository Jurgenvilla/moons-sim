# import modules that I'm using
import matplotlib
matplotlib.use('TKAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as pltlib
import tkinter
from tkinter import filedialog
from tkinter import *
import numpy as np
import scipy.io as sc
from scipy import ndimage
import math
import sys
import os.path
from scipy.ndimage.filters import gaussian_filter

import csv
#import matplotlib.pyplot as pltlib
# lmfit is imported becuase parameters are allowed to depend on each other along with bounds, etc.
from lmfit import minimize, Parameters, Minimizer

#detector_model_RI=str(sys.argv[1])
#fields = 'moons_mode', 'exptime', 'ab_mag', 'template', 'fwhm', 'inwl', 'infl', 'seeing', 'z', 'galsize', 'stray', 'skyres'
field_names = 'MOONS mode', 'Channel', 'Atm. Diff. Wavelength [um]','Magnitude [AB]','Extended', 'Template', 'Em. line FWHM [km/s]','Em. rest wavelength [um]', 'Em. Flux 10^-16 erg/s/cm2', 'Redshift', 'Seeing [arcsec]', 'Airmass', 'Stray light [%]', 'Sky residual [%]','NDIT','DIT [s]'
def_vals= 'mode','channel','1.2','20.0','extended','template','','','','0.0','0.8','1.2','1.0','0.0','1','300'
simmode='batch'

moons_mode_input=str(sys.argv[2])
channel_input=str(sys.argv[3])
atm_dif_ref_input=float(sys.argv[4])
magnitude_input=float(sys.argv[5])
extended_input=float(sys.argv[6])
template_input=str(sys.argv[7])
emission_lines_input=float(sys.argv[8])
emission_lines_wav_input=float(sys.argv[9])
emission_lines_flux_input=float(sys.argv[10])
redshift_input=float(sys.argv[11])
seeing_input=float(sys.argv[12])
airmass_input=float(sys.argv[13])
stray_input=float(sys.argv[14])
sky_res_input=float(sys.argv[15])
ndit_input=float(sys.argv[16])
dit_input=float(sys.argv[17])
var_input=[moons_mode_input,channel_input,atm_dif_ref_input,magnitude_input,extended_input, template_input, emission_lines_input, emission_lines_wav_input,emission_lines_flux_input,redshift_input,seeing_input,airmass_input,stray_input,sky_res_input,ndit_input,dit_input]

def get_band(ent):
    return ent.get()

def get_diffraction(self,Lambda,airmass,atm_ref_wav):
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

def getEntries_batch(self,fields):
    entries = []
    bands= ['RI', 'YJ', 'H']
    Temps= ['','MARCS T:3500K logg:1.5','MARCS T:3500K logg:4.5','MARCS T:4000K logg:1.5','MARCS T:4000K logg:4.5','MARCS T:4500K logg:1.5','MARCS T:4500K logg:4.5','MARCS T:5000K logg:1.5','MARCS T:5000K logg:4.5','MARCS T:5500K logg:1.5','MARCS T:5500K logg:4.5','SSP age:11Gyr solar met','SSP age:2.5Gyr solar met','SSP age:1.4Gyr solar met','SSP age:100Myr solar met','Constant flux']
    respow= ['High resolution','Low resolution']
    airm= ['1.0','1.2','1.4','1.6','1.8']
    inde=0
    for field in field_names:
        var=var_input[inde]
        entries.append((field,var))
        inde=inde+1
    return entries

def initialize_batch(self):
    if (simmode=='batch'):
        ents=self.getEntries_batch(field_names)
        try:
            self.fetch(ents)
        except:
            print('ERROR in input variables for batch mode. Check Manual')

def refreshFigure_batch(self,sn_cont_res_e,npix_e,outnoise_e,outputwl_e,sp_src_e,sn_cont_e,sn_cent,sky_spec,ins_through,atm_transm,res_power,sim_spectrum):
    x_print=outputwl_e/1.0e4
    y_print=sn_cont_e
    sens_out_file="Sensitivity_table.txt"
    name_x=['# Lambda']
    name_y=['SNR']
    with open(sens_out_file,'w+') as f:
        writer = csv.writer(f,delimiter='\t')
        writer.writerows(zip(name_x,name_y))
        writer.writerows(zip(x_print,y_print))

def CheckTemplateLocal_bash(self,template):
    localfile_exists=os.path.isfile(template)
    return localfile_exists

def CheckEntries(self,entries):
    print(" ")
    if (entries[0][1].get() == '' or entries[0][1].get() == "Select MOONS mode"):
        print("ERROR: Please select one the MOONS mode options")
        return False
    if (entries[1][1].get() == '' or entries[1][1].get() == "Select band"):
        print(" ")
        print("Please select one the Channel/band options")
        print(" ")
        return False
    if (entries[5][1].get() == "Select Template" or entries[5][1].get() ==''):
        try:
            template_local_path
            return True
        except NameError:
            print("Please select one the available templates or load a local file")
            return False
    if (entries[11][1].get() == '' or entries[11][1].get() == "Select airmass"):
        print("Please select the airmass of observation")
        return False
    if (entries[15][1].get() == ''):
        print("Exposure time input is not valid. Please check")
        return False
    try:
        ab_mag=float(entries[3][1].get())
        check=True
    except ValueError:
        print("Magnitude input is not valid. Please check")
        return False
    try:
        seeing=float(entries[10][1].get())
        check=True
    except ValueError:
        print("Seeing input is not valid. Please check")
        return False
    if (entries[6][1].get() != '' and entries[7][1].get() != '' and entries[8][1].get() != ''):
        try:
            fwhm=np.array([float(entries[6][1].get())])
            check=True
        except ValueError:
            print("FWHM input is not valid. Please check")
            return False
        try:
            inwl=np.array([float(entries[7][1].get())])
            check=True
        except ValueError:
            print("Emission line wavelength input is not valid. Please check")
            return False
        try:
            infl=np.array([float(entries[8][1].get())])
            check=True
        except ValueError:
            print("Emission line flux input is not valid. Please check")
            return False
    if (entries[2][1].get() != ''):
        try:
            diff_wave=float(entries[2][1].get())
            check=True
        except ValueError:
            print("Differential Atmospheric dispersion wavelength input is not valid. Please check or leave blank")
            return False
    if (entries[9][1].get() != ''):
        try:
            z_val=float(entries[9][1].get())
            check=True
        except ValueError:
            print("Redshift input is not valid. Please check or leave blank")
            return False
    if (entries[4][1].get() != ''):
        try:
            galsize=float(entries[4][1].get())
            check=True
        except ValueError:
            print("Source size (in arcsec) value is not valid, please check or leave blank ")
            return False
    if (entries[12][1].get() != ''):
        try:
            stray=float(entries[12][1].get())
            check=True
        except ValueError:
            print("Stray light contamination (percentage) is not a valid input. Please check")
            return False
    if (entries[14][1].get() != ''):
        try:
            N_dit=float(entries[14][1].get())
            check=True
        except ValueError:
            print("Number of DIT exposures input is not a valid value. Please check")
            return False
    try:
        exptime=float(entries[15][1].get())
        check=True
    except ValueError:
        print("Exposure time input is not a valid value. Please check")
        return False
    if check:
        return True
    else:
        return False
    print(" ")

def fetch(self,entries):
    checkVals=True
    if checkVals:
        fs=1.05
        fs_fieldstop=1.2
        fr=1.02778 #1.04
        #eff=0.20
        npix=4096
        npix_1=npix-1
        # Get entries from GUI. When not provided, assuming default.
        # MOONS MODE selection
        if (simmode=='batch'):
            ab_mag=float(entries[3][1])
            band=str(entries[1][1])
            airmass=str(entries[11][1])
            moons_mode=str(entries[0][1])
            exptime=float(entries[15][1])
            template=str(entries[5][1])
            seeing=float(entries[10][1])
            airmass_fl=float(entries[11][1])
        print(".............................................")
        if (simmode=='batch'):
            if (moons_mode == "HR"):
                moons_mode="High resolution"
            if (moons_mode == "LR"):
                moons_mode="Low resolution"
        if (moons_mode == "High resolution"):
            if ( band == "H" ) :
                print("Adopting moons mode: ", moons_mode, band)
                wlr=[1.521,1.641]
                ab = float(ab_mag)# + 0.04*1.2
                RON = 3.
                DK = 4#30.
                print("Setting wavelength range to ", str(wlr))
                template_wl_norm=1.630
                QE_file=self.resource_path('QE_4RG.txt')
                QE_wav, QE_eff  = np.loadtxt(QE_file, unpack = True)
            if ( band == "YJ" ) :
                print("Adopting moons mode: ", moons_mode, band)
                wlr=[0.934,1.350]
                ab=float(ab_mag) #+0.05*1.2
                RON = 3.
                DK = 4#30.
                print("Setting wavelength range to ", str(wlr))
                template_wl_norm=1.22
                QE_file=self.resource_path('QE_4RG.txt')
                QE_wav, QE_eff  = np.loadtxt(QE_file, unpack = True)
            if ( band == "RI" ) :
                print("Adopting moons mode: ", moons_mode, band)
                wlr=[0.765,0.898]
                ab = float(ab_mag) # + 0.06 * 1.2
                DK = 2.
                print("Setting wavelength range to ", str(wlr))
                template_wl_norm=0.797
                #setting detector_model_RI to 'LBNL':
                QE_file=self.resource_path('QE_LBNL.txt')
                QE_wav, QE_eff  = np.loadtxt(QE_file, unpack = True)
                RON=3.0
        if (moons_mode == "Low resolution"):
            if ( band == "H" ) :
                print("Adopting moons mode: ", moons_mode, band)
                wlr=[1.452,1.800]
                ab=float(ab_mag) #+0.04*1.2
                RON = 3.
                DK = 4#30.
                print("Setting wavelength range to ", str(wlr))
                template_wl_norm=1.63
                QE_file=self.resource_path('QE_4RG.txt')
                QE_wav, QE_eff  = np.loadtxt(QE_file, unpack = True)
            if ( band == "YJ" ) :
                print("Adopting moons mode: ", moons_mode, band)
                wlr=[0.934,1.350]
                ab=float(ab_mag) #+0.05*1.2
                RON = 3.
                DK = 4#30.
                print("Setting wavelength range to ", str(wlr))
                template_wl_norm=1.22
                QE_file=self.resource_path('QE_4RG.txt')
                QE_wav, QE_eff  = np.loadtxt(QE_file, unpack = True)
            if ( band == "RI" ) :
                print("Adopting moons mode: ", moons_mode, band)
                wlr=[0.647,0.934]
                ab = float(ab_mag) #+ 0.06 * 1.2
                DK = 2.
                print("Setting wavelength range to ", str(wlr))
                template_wl_norm=0.797
                #setting detector_model_RI to 'LBNL':
                QE_file=self.resource_path('QE_LBNL.txt')
                QE_wav, QE_eff  = np.loadtxt(QE_file, unpack = True)
                RON=3.0

        print(".............................................")
        if (simmode=='batch'):
            print('Checking for Source Template locally')
            foundfilelocally=self.CheckTemplateLocal_bash(template)
            if foundfilelocally:
                print('Using user provided template: %s'%template)
                templateislocal=True
            else:
                print('Template file not found!!!')
        cen_wav=(wlr[1]+wlr[0])/2.0*1.0e3
        wav_range_length=(wlr[1]-wlr[0])*1.0e3
        atm_dif_ref=float(entries[2][1])
        Lambda=np.arange(0.5,1.9,0.1)
        atm_diff=self.get_diffraction(Lambda, airmass_fl, atm_dif_ref)
        r_0=0.100/seeing*(cen_wav/500.0)**(1.2)*airmass_fl**(-0.6)
        F_kolb=-0.981644
        fwhm_atm=seeing*airmass_fl**(0.6)*(cen_wav/500.0)**(-0.2)*np.sqrt(1.0+F_kolb*2.183*(r_0/46.0)**(0.356))
        fwhm_tel=0.000212*(cen_wav/8.1)
        fwhm_iq=np.sqrt(fwhm_tel**2+fwhm_atm**2)
        seeing=fwhm_iq
        print(" ")
        print('Expected Image Quality in selected band:', round(seeing,2))
        print(" ")
        if (simmode=='batch'):
            z_val=float(entries[9][1])
            galsize=float(entries[4][1])
            skyres=float(entries[13][1])
            N_dit=float(entries[14][1])
            stray=float(entries[12][1])
        if (simmode=='batch'):
            if (float(entries[6][1]) == 0.0 or float(entries[7][1]) == 0.0 or float(entries[8][1]) == 0.0):
                emission_lines=False
                print(" ")
                print("No emission lines (if needed then FWHM, wavelength, and flux are mandatory)")
                print(" ")
            else:
                fwhm=np.array([float(entries[6][1])])
                inwl=np.array([float(entries[7][1])])
                infl=np.array([float(entries[8][1])])
                inwl=inwl*(1.0+float(z_val))
                if (inwl>wlr[0] and inwl<wlr[1]):
                    emission_lines=True
                else:
                    emission_lines=False
                    print("WARNING: Provided rest frame wavelength of emission line does not fall in the selected wavelength setup at input redshift")
        dit=float(exptime)
        DKsec = DK/3600.
        sky_pix = 0.37/fr # pixel size on sky
        r_fiber_pix = (fs/2.)/sky_pix # aperture of the fiber (radius) in pix
        disp = 1.0e4*(wlr[1]-wlr[0])/(float(npix)-1.0) # spectral dispersion in A/pix
        print("")
        print("Dispersion: ",str(round(disp,2)),' A/pix')
        print(" ")
    ##################################################################
    # Extracting the relevant OH template
        #data_dir='data_dir/Skymodel'
        if (disp > 0.65):
    	    #oh_f=idlsave.read('oh_th_spec_reb_lowres.sav')
            #oh_f=sc.readsav('oh_th_spec_reb_lowres.sav',verbose=False)
            skyfile=self.resource_path("skytable_lr_a"+airmass+".sav")
            oh_f=sc.readsav(skyfile,verbose=False)
        else:
    	    #oh_f=sc.readsav('Skymodel/skytable_hr_a'+airmass+'.sav',verbose=False)
            skyfile=self.resource_path("skytable_hr_a"+airmass+".sav")
            oh_f=sc.readsav(skyfile,verbose=False)
        rwl0=oh_f.rwl0*1.0e4 # wavelength of sky model in Angstroms
        rdwl=oh_f.rdwl
        #rfn0=oh_f.rfn0*math.pi*(8.1e2/2.0)**2/1.0e4 # OH lines + sky continuum from ESO skycalc model scaled to VLT aperture e-/s/m2/arcsec2/A
        rfn0=oh_f.rfn0*math.pi*(8.1/2.0)**2*rdwl/1.0e4# * 0.8
        atmtr=oh_f.atmtr # atmospheric transmission
        atmwl=oh_f.atmwl*1.0e4
        tol = (wlr[1]-wlr[0])*0.05
        min_wl=(wlr[0]-tol)*1.0e4
        max_wl=(wlr[1]+tol)*1.0e4
        iok = np.where((rwl0 > min_wl) & (rwl0 < max_wl))[0]
        rfn = rfn0[iok] # flux of OH spectrum
        rwl = rwl0[iok] # wave of OH spectrum
        # normalize the OH to fiber size input
        rfn_sky = rfn*math.pi*(fs_fieldstop/2.0)**2 # OH spectrum ALONE
    ##################################################################
    # Set up emission lines
        if emission_lines:
            def gaussian(x, mu, sig):
                return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
            infles=infl*38.3*inwl/1.65 # emission line in e/sec
            sp_src=np.zeros(np.size(rwl))
            inwlA=inwl*1.0e4
            for iline in range(0,np.size(inwl)):
    	        i=np.where(rwl >= inwlA[iline])[0]
    	        sp_src[min(i)]=infles[iline]
            fwhm_A = inwlA*fwhm/2.99792e5 # FWHM in A
            fwhm_rpix = fwhm_A#/rdwl
            sigma_rpix = fwhm_rpix/2.35
            #print(sigma_rpix)
            #print(np.mean(sigma_rpix))
            sp_src_sm = infles/(np.sqrt(2*math.pi*sigma_rpix**2))*gaussian(rwl, inwlA, sigma_rpix)
            #sp_src_sm = gaussian_filter(sp_src,np.mean(sigma_rpix))
        else:
            sp_src_sm=0.0
    ##################################################################
    # Extracting the relevant SCIENCE template from Library
        #template_lib=idlsave.read(template)
        if templateislocal:
            f2=open(template,"r")
            lines=f2.readlines()
            f2.close()
            wl_tpl_1=[]
            fl_tpl_1=[]
            for line in lines:
                p=line.split()
                wl_tpl_1.append(float(p[0]))
                fl_tpl_1.append(float(p[1]))
            wl_tpl=np.array(wl_tpl_1)
            fl_tpl=np.array(fl_tpl_1)
        else:
            template_lib=sc.readsav(template,verbose=False)
            wl_tpl=template_lib.wl_tpl # wavelenth of the source spectrum
            fl_tpl=template_lib.fl_tpl # flux of the source spectrum
        wl_tpl = wl_tpl*(1.+float(z_val)) # redshift the spectrum to provided z
    ##################################################################
    # normalising the template source spectrum to the magnitude
        i = np.where(wl_tpl >= template_wl_norm*1.0e4)[0]
        templ_fl_ref = fl_tpl[i[0]]
        fl_tpl_normalized_to_one = fl_tpl/templ_fl_ref
    ##################################################################
    # resampling into the OH spectrum. I do it in this direction to allow for any input from the user as science source template, leaving the OH sampling as the fixed one.
        fl_tpl_res=np.interp(rwl,wl_tpl,fl_tpl_normalized_to_one)
        wl_tpl=rwl
    # template photon flux per A in cgs assuming a telescope effective diameter of 8.2m
        templ_phot_fl = fl_tpl_res*3.63e3*10.0**(-(ab+57.5)/2.5)*math.pi*(8.2e2/2.0)**2/(rwl*6.626e-27)
    # template photon flux per PER SPECTRAL PIXEL OF THE OH GRID
        templ_phot_fl_pix = templ_phot_fl*rdwl
        templ_phot_fl_pix = templ_phot_fl_pix + sp_src_sm # Source flux in OH sampling + emission line in OH sampling
        #sp_src_sm=sp_src_sm+fib_templ  # Source flux in OH sampling + emission line in OH sampling
    ##################################################################
    #including fibre losses from seeing and target size for extended sources
    # Also add losses from Atmospheric Refraction by first resampling to OH spectrum and then applying "seeing effect"
        if (atm_dif_ref > 0.5):
            atm_diff_oh=np.interp(rwl,Lambda*1.0e4,atm_diff)
        if (galsize > 0):
            seeing=np.sqrt(float(seeing)**2)#+2*float(galsize/1.68)**2)
            eff_area=math.pi*(fs_fieldstop/2.0)**2
            seeing_cor_flux=templ_phot_fl_pix*eff_area # effectively independent of Seeing!!
        if (galsize == 0):
            if (float(z_val) > 0.0):
                seeing_arr=np.sqrt(seeing**2+atm_diff_oh**2)#+2*(0.3/1.68)**2) #assuming a hl=0.3 arcsec
            if (float(z_val) == 0.0):
                seeing_arr=np.sqrt(seeing**2+atm_diff_oh**2)#+float(galsize)**2)
            seeing_cor_flux= (1.0-1.0/np.exp( (np.sqrt(math.log(2.0)) * fs_fieldstop/seeing_arr)**2))*templ_phot_fl_pix
        templ_phot_fl_pix = seeing_cor_flux
        fib_frac=1.0 # allowing for tuning of fraction loss manually, keep as 1.0 for no additional loss
        fib_templ = fib_frac*templ_phot_fl_pix
    ##################################################################
        sp_src_sm=fib_templ  # Source flux in OH sampling + emission line in OH sampling including fibre losses
        oh_src_sp=rfn_sky+sp_src_sm # OH spectrum + Source (continuum+emission) spectrum
        detpx_rpx = disp/rdwl   # ratio between detector pixel and pixels in OH spectrum
        r_fib_rpix = r_fiber_pix*detpx_rpx  # Fibre radius in pixels of the OH spectrum

            # defining the fibre kernel function
#            k_x=np.arange(-50.0,51,1.0)
#            k_y= k_x-k_x # a zero array same as k_x
#            ik=np.where(np.absolute(k_x) <= r_fib_rpix)[0]
#            k_y[ik]=2.0*np.sqrt(r_fib_rpix**2-k_x[ik]**2)/(math.pi*r_fib_rpix**2) # Fibre Kernel
            # retrieving spectral resolving power
        rs_fwhm_rpix = 2.0*np.sqrt(r_fib_rpix**2-0.5**2) #Maximum FWHM is at 1 pixel in x direction
        rs_fwhm_A = rs_fwhm_rpix*rdwl
        samp = rs_fwhm_A/disp
        rs = 1.0e4*((wlr[1]+wlr[0])/2.0)/rs_fwhm_A
        sigma_spec=cen_wav*10.0/rs/samp # sigma for convolution
        sp_conv=ndimage.gaussian_filter1d(oh_src_sp, sigma_spec)
        sp_conv_src=ndimage.gaussian_filter1d(sp_src_sm, sigma_spec)
        sp_conv_sky=ndimage.gaussian_filter1d(rfn_sky, sigma_spec)

#            sp_conv=np.convolve(oh_src_sp,k_y,mode="same") # convolve Source (cont+emission)+OH spectrum
#            sp_conv_src=np.convolve(sp_src_sm,k_y,mode="same") # convolve only Source+emission spectrum
#            sp_conv_sky=np.convolve(rfn_sky,k_y,mode="same") # convolve only OH spectrum

        # Resample source to detector pixels, conserving total photons
        pix_arr=np.arange(0,npix,1)
        outputwl = wlr[0]*1.0e4+pix_arr*disp
        sp_det_src = np.interp(outputwl,rwl,sp_conv_src)
        inband = np.where((rwl >= outputwl[0]) & (rwl <= outputwl[npix_1]))[0]
        sp_conv_src_sum=sp_conv_src[inband].sum()
        sp_det_src_sum=sp_det_src.sum()
        renorm=sp_conv_src_sum/sp_det_src_sum
        sp_det_src_rn = sp_det_src*renorm # for source only
    # Resample sky to detector pixels, conserving total photons
        sp_det_sky = np.interp(outputwl,rwl,sp_conv_sky)
        sp_conv_sky_sum=sp_conv_sky[inband].sum()
        sp_det_sky_sum=sp_det_sky.sum()
        renorm_sky=sp_conv_sky_sum/sp_det_sky_sum
        sp_det_sky_rn = sp_det_sky*renorm_sky # for sky only
        sp_det = sp_det_src_rn+sp_det_sky_rn # total detector flux
        sp_det_c = sp_det
        sp_det_src_c = sp_det_src_rn
            # Add background sky and wings of OH lines (from Roberto's adoption via Tino's suggestions)
        wing_fact = 0.5
        sp_wings = sp_det_c #gaussian_filter(sp_det_c,wing_fact)
        sp_wings_src = sp_det_src_c #gaussian_filter(sp_det_src_c,wing_fact) # this was commented OUT in v2.0
        sp_wings_sky = sp_det_sky_rn #gaussian_filter(sp_det_sky_rn,wing_fact)
            #here is to add an additional background constant level. Not used
        between_oh_back = 0 #1.5 # to be refined according to MOON phase
        add_back= 0 #between_oh_back*disp*math.pi*(fs/2)**2
        #I do not use these values as the skycalc model has background included.
        sp_totback = sp_wings+add_back
        sp_totback_sky = sp_wings_sky+add_back
        ################################################
        # Add transmission filter from atmosphere
        #atmos=sc.readsav('atmtrans.sav',verbose=False)
        #atmtr=atmos.atmtr
        #atmwl=atmos.atmwl
        # Add Telescope transmission
        telescope_file=self.resource_path('telescope_eff.txt')
        telescope_wav_micron, telescope_eff  = np.loadtxt(telescope_file, unpack = True)
        telescope_wav=telescope_wav_micron*1.0e04
        # Resample atmospheric transmission on wlgrid
        atminterp = np.interp(outputwl,atmwl,atmtr)
        # Resample telescope transmission on wlgrid
        tel_eff = np.interp(outputwl,telescope_wav,telescope_eff)
        # Add transmission curve from overall instrument+telescope efficiency
        #data_dir='data_dir/Inst_setup'
        trans_file=self.resource_path("throughput.sav")
        efficiency=sc.readsav(trans_file,verbose=False)
        eff_lr0=efficiency.lreff
        eff_hr0=efficiency.hreff
        eff_wl0=efficiency.weff*10.0
        effok = np.where((eff_wl0 > min_wl) & (eff_wl0 < max_wl))[0] #Selecting instrument transmission in Setup
        eff_lr=eff_lr0[effok]
        eff_hr=eff_hr0[effok]
        eff_wl=eff_wl0[effok]
        # Resample instrument efficiency on wlgrid and detector QE (updated to detector FDR)
        if (moons_mode == "High resolution"):
            eff = np.interp(outputwl,eff_wl,eff_hr)
#               if (band == 'RI') :
            QE_detector=np.interp(outputwl,QE_wav*1.e4,QE_eff)
            eff = QE_detector/100.0*eff*tel_eff
        if (moons_mode == "Low resolution"):
            eff = np.interp(outputwl,eff_wl,eff_lr)
#               if (band == 'RI') :
            QE_detector=np.interp(outputwl,QE_wav*1.e4,QE_eff)
            eff = QE_detector/100.0*eff*tel_eff
        # and at wavelengths of emission lines (this should be improved to take the average transmission across the line)
        if emission_lines:
            atmlines = np.interp(inwl,atmwl,atmtr)
        # scale by efficiency and atmosphere
        sp_eff = sp_totback*eff*atminterp
        # the following is for the source alone
        sp_eff_src = sp_wings_src*eff*atminterp
        # the following keeps track of the individual line fluxes
        if emission_lines:
            flux_line_eff = infles*eff*atmlines
        # the following is for the sky alone
        sp_eff_sky = sp_totback_sky*eff*atminterp
        # add dark current
        # number of pixels COLLAPSED ALONG Y DIRECTION
        npix_y = r_fiber_pix*2.0
        sp_dk = sp_eff_src+DKsec*npix_y # changed from sp_eff to sp_eff_src to only count on Source Signal without sky
        sp_dk_sky = sp_eff_sky+DKsec*npix_y
        # scale by integration time
        outspec = sp_dk*dit
        # the following is for the sky only
        outspec_sky = sp_dk_sky*dit
        # the following is for the source alone
        sp_src = sp_eff_src*dit
        # the following keeps track of the individual line fluxes
        if emission_lines:
            flux_line_dit = flux_line_eff*dit
        #####################################################
        # Add stray light contribution
        total_stray=outspec_sky.sum()*N_dit*500.0/npix**2
        outspec_stray=total_stray*float(stray)/100.0*npix_y

            #outnoise = np.sqrt(outspec+outspec_sky+outspec_stray+RON**2*npix_y*2.) # if no stray then outspec_stray=0.0
        # sky subtraction residual IN PERCENTAGE
            #outnoise=outnoise+float(skyres)/100.0*outspec_sky

        #####################################################
        #detector NOISE:
        noisedet=np.sqrt(npix_y*(N_dit*RON**2+DKsec*dit*N_dit))
        #background NOISE:
        noiseback=np.sqrt(sp_eff_sky*dit*N_dit+outspec_stray)
        #residual sky subtraction:
        noiseskyres=float(skyres)/100.0*outspec_sky
        #total NOISE:
        noisetot=np.sqrt(noiseback**2+noisedet**2+sp_src*N_dit)
        outnoise=noisetot+noiseskyres
        print("#########################################")
        print(" ")
        print("RESULTS:")
        print(" ")
        print("Spectral resolution: ",int(rs))
            #print("Spectral sampling: ",round(samp,2))
        # estimates the S/N on lines by estimating the noise within twice the width of the line (for the moment given by the sum in quadrature of the spectral resolution and intrinsic line width)
        sn_cont_all = sp_src/outnoise*N_dit#*np.sqrt(npix_y*2) #extracting spectra over 2*FWHM pixels
        if emission_lines:
            tot_fwhm = np.sqrt(rs_fwhm_A*rs_fwhm_A+fwhm_A*fwhm_A)
            sn=np.zeros(np.size(inwl))
            print(" ")
            print("SN from emission line:")
            print(" ")
            #print("Number of lines:",np.size(inwl))
            for iline in range(0,np.size(inwl)):
                i=np.where((outputwl >= inwl[iline]*1.0e4-tot_fwhm[iline]) & (outputwl <= inwl[iline]*1.0e4+tot_fwhm[iline]))[0]
                #tot_N=np.sqrt(np.sum(outnoise[i]**2))
                sn[iline]=np.max(sn_cont_all[i]) #flux_line_dit[iline]/tot_N
                print(round(inwl[iline],4),"um, S/N = ",round(sn[iline],2))
        # estimates median S/N on continuum
            #sn_cont_all = sp_src/outnoise*np.sqrt(N_dit)
        if emission_lines:
            for iline in range(0,np.size(inwl)):
                i=np.where((outputwl <= inwl[iline]*1.0e4-3*tot_fwhm[iline]) | (outputwl >= inwl[iline]*1.0e4+3*tot_fwhm[iline]))[0]
                sn_cont = np.median(sn_cont_all[i])
        else:
            sn_cont = np.median(sn_cont_all)
        # resampling to resolution element
        sn_cont_res = sn_cont*np.sqrt(3.0)
        cent_range=np.where((outputwl >= (cen_wav-wav_range_length*0.02)*10.0) & (outputwl <= (cen_wav+wav_range_length*0.02)*10.0))[0]
        sn_central=sn_cont_all[cent_range].max()
        print(" ")
        #print("Continuum noise contribution from detector: ",round(np.median(noisedet/noisetot)*100.0,2))
        #print("Continuum noise contribution from background: ",round(np.min(noiseback/noisetot)*100.0,2))
        print(" ")
        print('Median S/N on continuum per pixel step in wavelength at R~%s = %s'%(int(rs),round(sn_cont,2)))
        print(" ")
        print('S/N at central wavelength per pixel at R~%s = %s'%(int(rs),round(sn_central,2)))
        print(" ")
        print("Median S/N on continuum per resolution element = ",round(sn_cont_res,2))
        #print('Resolution element of %s pixels'%round(2.7,1)) #samp=2.7
        print(" ")
        print("#########################################")
        # Creates the output observed spectrum
        normnoise=np.random.random((npix))*2.0-1.0 #Pixel Array with values between -1 and 1
        res_noise=normnoise*outnoise*np.sqrt(N_dit)
        #sim_spectrum=(sp_src*N_dit+res_noise)/(dit*N_dit)/atminterp # Perfectly telluric corrected, sky subtracted spectrum
        sim_spectrum=sp_src
        if (simmode=='batch'):
            self.refreshFigure_batch(sn_cont_res,npix,outnoise,outputwl,sp_src,sn_cont_all,sn_central,outspec_sky,eff,atminterp,rs,sim_spectrum)
        else:
            self.refreshFigure(sn_cont_res,npix,outnoise,outputwl,sp_src,sn_cont_all,sn_central,outspec_sky,eff,atminterp,rs,sim_spectrum)
    else:
        print("Error in input parameters, please check your input values.")
if __name__ == "__main__":
    if (simmode=='batch'):
        run=initialize_batch()
