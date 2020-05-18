import sys, os
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from scipy.special import j1
import matplotlib.pyplot as plt
from scipy.constants import c
import pandas as pd

point_centre = ['12h36m55s','+62d14m15s']
pbtype = "Bessel"
D = 25.48 
return_rms = True
return_pb_corr = True
slow_load = True
def Bessel_pb(D,wl,sep):
	insert = ((D*np.pi)/wl)*np.sin(sep)
	I = ((2*j1(insert))/insert)**2
	return I

fitsfiles=sys.argv[1:]
if fitsfiles == []:
	print('No fitsfiles specified')
	sys.exit()
else:
	print('primary beam correcting: %s' % ", ".join(fitsfiles))

rms = []
rms_fits = []
ra = []
dec = []
pb_cor = []
pb_cor_test = []
for i in fitsfiles:
	hdu = fits.open(i)
	header = hdu['PRIMARY'].header
	wcs = WCS(header)
	if slow_load == True:
		x = np.arange(1,int(header['NAXIS1'])+1,1)
		for j in range(int(header['NAXIS2'])):
			if j % 200 == 0:
				print("Correcting row %d/%d"%(j,header['NAXIS2']))
			xv, yv = np.meshgrid(x,j+1,sparse=True, indexing='xy')
			c0 = wcs.celestial.wcs_pix2world(xv,yv,1)
			c1 = SkyCoord(c0[0],c0[1],unit=('deg','deg'))
			pc = SkyCoord(point_centre[0],point_centre[1])
			sep = np.float64(c1.separation(pc).rad)
			wl = c/header['CRVAL3']
			if pbtype == "Bessel":
				I = Bessel_pb(D,wl,sep)
			else:
				print('Please specify pb correction type')
				sys.exit()
			hdu['PRIMARY'].data[:,:,j,:] = hdu['PRIMARY'].data[:,:,j,:]/I
			if j == 0:
				if return_pb_corr == True:
					pb_cor_test.append(I[0,0])
			if j == int(header['NAXIS2']/2):
				if return_pb_corr == True:
					pb_cor.append(I[0,int(I.shape[0]/2)])

	else:
		x = np.arange(1,int(header['NAXIS1'])+1,1)
		y = np.arange(1,int(header['NAXIS2'])+1,1)
		xv, yv = np.meshgrid(x, y, sparse=True, indexing='xy')
		c0 = wcs.celestial.wcs_pix2world(xv,yv,1)
		c1 = SkyCoord(c0[0],c0[1],unit=('deg','deg'))
		pc = SkyCoord(point_centre[0],point_centre[1])
		sep = np.float64(c1.separation(pc).rad)
		wl = c/header['CRVAL3']
		if pbtype == "Bessel":
			I = Bessel_pb(D,wl,sep)
		else:
			print('Please specify pb correction type')
			sys.exit()
		hdu['PRIMARY'].data = hdu['PRIMARY'].data/I
		if return_pb_corr == True:
			pb_cor.append(I[int(I.shape[0]/2),int(I.shape[0]/2)])
			pb_cor_test.append(I[0,0])

	hdu.writeto('%s_PB.fits'%i.split('.fits')[0],overwrite=True)
	if return_rms == True:
		rms_fits.append(i)
		rms.append(np.std(hdu['PRIMARY'].data.squeeze()[0:250,0:250]))
		ra.append(header['CRVAL1']) 
		dec.append(header['CRVAL2'])
	hdu.close()

if return_rms == True and return_pb_corr==False:
	pd.DataFrame({'fitsfile':i,'ra':ra,'dec':dec,'rms':rms}).to_csv('rms_vals.csv')
elif return_rms == True and return_pb_corr==True:
	pd.DataFrame({'fitsfile':i,'ra':ra,'dec':dec,'rms':rms,'pbcor':pb_cor,'pbcorner':pb_cor_test}).to_csv('rms_vals.csv')


