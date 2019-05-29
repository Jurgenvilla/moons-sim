import numpy as np
from astropy.io import fits
skyname=np.loadtxt('list_fits',dtype=str, unpack = True)
contador=0
airmass=['a1.0','a1.2','a1.4','a1.6','a1.8','a2.0']
for name in skyname:
    hdu_sky=fits.open(name+'.fits')
    rwl0=hdu_sky[1].data['lam']*1.e4
    rfn0=hdu_sky[1].data['flux']
    atmwl=hdu_sky[1].data['lam']
    atmtr=hdu_sky[1].data['trans']
    rdwl=0.1
    comments=hdu_sky[0].header['COMMENT'][30:]
    hdu1 = fits.PrimaryHDU()
    hdu1.data=rfn0/((1.e2)**2)/1.e4 # 1e2 cm = 1 m 1.e4 A = 1 um
    hdu1.header['CDELT1']=rdwl
    hdu1.header['CRVAL1']=rwl0[0]
    hdu1.header['TUNIT1']='Angstroms'
    hdu1.header['MTYPE'] = 'Skymodel'
    hdu1.header['MNAME'] = 'ESO_skycalc'
    hdu1.header['TUNIT2']='ph/s/cm2/A/arcsec2'
    hdu1.header['R']='60000-180000'
    hdu1.header['SAMPLING']=1
    for comms in comments:
        hdu1.header.add_comment(comms)
    hdu2 = fits.ImageHDU()
    hdu2.data=atmtr
    hdu2.header['CDELT1']=rdwl
    hdu2.header['CRVAL1']=rwl0[0]
    hdu2.header['TUNIT1']='Angstroms'
    hdu2.header['MTYPE'] = 'Transmission'
    hdu2.header['MNAME'] = 'ESO_skycalc'
    hdu2.header['TUNIT2']='sky transmission fraction'
    hdu2.header['R']='60000-180000'
    hdu2.header['SAMPLING']=1
    for comms in comments:
        hdu2.header.add_comment(comms)
    new_hdu=fits.HDUList([hdu1,hdu2])
    new_hdu.writeto('SkyTemplate_ESO_'+airmass[contador]+'.fits',clobber=True,overwrite=True)
    contador=contador+1
print('done')
