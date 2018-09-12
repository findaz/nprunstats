#!/usr/bin/env python

import os
import numpy as np
import fitsio

from astrometry.util.fits import fits_table
from tractor.psfex import PsfExModel, PixelizedPsfEx

class DepthStats():

    def skynoise(self, image, blanton = True):

        if not blanton:
            '''
            Mask sources and find rms
            '''
            raise NotImplementedError

        '''
        Estimate per-pixel noise via Blanton's 5-pixel MAD
        '''

        slice1 = (slice(0,-5,10),slice(0,-5,10))
        slice2 = (slice(5,None,10),slice(5,None,10))
        mad = np.median(np.abs(image[slice1] - image[slice2]).ravel())
        return 1.4826 * mad / np.sqrt(2.)
    
    def gauss_psf_depth(self, skynoise, fwhmpsf, magzero = 0.0):
    
        '''
        depth from gaussian approximation to the psf 
        '''
            
        skynoise/=10.**((magzero- 22.5) / 2.5)
        psf_sigma = fwhmpsf / 2.35
        gnorm = 1./(2. * np.sqrt(np.pi) * psf_sigma)
        detsig1 = skynoise / gnorm
        depth = 5. * detsig1
        # nanomaggies to magnitudes
        return -2.5 * (np.log10(depth) - 9)

    def psf_depth(self, skynoise, psf, magzero = 0.0):
    
        '''
        psf depth
        '''

        skynoise/=10.**((magzero- 22.5) / 2.5)
        detsig1 = skynoise / np.mean(psf)
        depth = 5. * detsig1
        # nanomaggies to magnitudes
        return -2.5 * (np.log10(depth) - 9)
    
class NinetyPrimeDepthStats(DepthStats):
    def __init__(self, imagename, sxdir = os.getcwd(), cpdir = os.getcwd()):

        self.img_name = os.path.basename(imagename)
        self.wht_name = self.img_name.replace('_ooi_', '_oow_')
        self.msk_name = self.img_name.replace('_ooi_', '_ood_')
        # FIX ME -- compress cat and psf file by default
        self.psf_name = self.img_name.replace('.fz','').replace('.fits', '.psf.fits')
        self.cat_name = self.img_name.replace('.fz','').replace('.fits', '.cat.fits')
        
        self.cpdir = cpdir
        self.sxdir = sxdir
        
        try:
            self.img_hdus = fitsio.FITS( os.path.join(self.cpdir, self.img_name) )
        except Exception:
            self.img_hdus = None
        try:    
            self.wht_hdus = fitsio.FITS( os.path.join(self.cpdir, self.wht_name) )
        except Exception:
            self.wht_hdus = None
        try:
            self.msk_hdus = fitsio.FITS( os.path.join(self.cpdir, self.msk_name) )
        except Exception:
            self.msk_hdus = None
        try:
            self.psf_hdus = fitsio.FITS( os.path.join(self.sxdir, self.psf_name) )
        except:
            self.psf_hdus = None
        try:
            self.cat_hdus = fitsio.FITS( os.path.join(self.sxdir, self.cat_name) )
        except:
            self.cat_hdus = None
        
        self.ampNum = [ [ 4, 2, 3, 1 ],  # CCD1 
                        [ 7, 5, 8, 6 ],  # CCD2 
                        [ 10,12,9,11 ],  # CCD3 
                        [ 13,15,14,16] ] # CCD4 
        
        self.aspp = 0.47 # plate scale
        
        self.coeffs = dict(
            g = [-0.00672, +0.00958, +0.06630, 0.0],
            r = [-0.00563, +0.01100, -0.04836, 0.0])
    
    @property
    def filtname(self):
        try:
            return self.img_hdus[0].read_header()['FILTER'].strip()[-1]
        except Exception:
            return None

    def get_invvar_map(self, ext):
        
        wht = self.wht_hdus[ext].read()
        msk = self.msk_hdus[ext].read()
        
        wht[np.isfinite(wht) == False] = 0.
        wht[np.isnan(wht)] = 0.
        wht[msk != 0] = 0.

        wht[wht > 0.0] = wht[wht > 0.0]
        
        return wht

    def get_sxcat_stds(self, ext, fluxName = 'flux_psf'):
        
        cat = fits_table( self.cat_hdus[ext] )
        return cat[ (cat.flags < 1) & (cat.get(fluxName) > 0.0) ]
        
    def psf_fwhm(self, ext):

        return self.psf_hdus[ext].read_header()['PSF_FWHM']
        
    def psf(self, ext):

        path = os.path.join(self.sxdir, self.psf_name)
        return PixelizedPsfEx(path)

def init_stats_table(imgnames):
        
    try:
        imgnames+[]
    except TypeError:
        imgnames = [imgnames]

    nrows=len(imgnames)
    tab=fits_table()
    tab.filename=np.array([os.path.basename(name) for name in imgnames])
    tab.set('object', np.zeros(nrows, np.int16))
    tab.ra=np.zeros(nrows)
    tab.dec=np.zeros(nrows)
    tab.set('filter',np.repeat(''.ljust(4),nrows))
    tab.exptime=np.zeros(nrows)
    tab.set('utc-obs', np.repeat(''.ljust(12),nrows))
    tab.ha=np.repeat(''.ljust(9),nrows)
    tab.jd=np.zeros(nrows)
    tab.airmass=np.zeros(nrows)
    tab.gain=np.zeros((nrows,4,4))
    tab.readnoise=np.zeros((nrows,4,4))
    tab.camtemp=np.zeros(nrows)
    tab.dewtemp=np.zeros(nrows)
    tab.magzptccd=np.zeros((nrows,4))
    tab.magzpt=np.zeros(nrows)
    tab.magzptrmsccd=np.zeros((nrows,4))
    tab.magzptrms=np.zeros(nrows)
    tab.raoffccd=np.zeros((nrows,4))
    tab.raoff=np.zeros(nrows)
    tab.decoffccd=np.zeros((nrows,4))
    tab.decoff=np.zeros(nrows)
    tab.nstdsccd=np.zeros((nrows,4), np.int16)
    tab.nstds=np.zeros(nrows, np.int16)
    tab.fwhmpsfccd=np.zeros((nrows,4))
    tab.fwhmpsf=np.zeros(nrows)
    tab.fwhmimgccd=np.zeros((nrows,4))
    tab.fwhmimg=np.zeros(nrows)
    tab.sigmaccd=np.zeros((nrows,4))
    tab.sigma=np.zeros(nrows)
    tab.skynoiseccd=np.zeros((nrows,4))
    tab.skynoise=np.zeros(nrows)
    tab.gausspsfdepthccd=np.zeros((nrows,4))	
    tab.gausspsfdepth=np.zeros(nrows)	
    tab.psfdepthccd=np.zeros((nrows,4))	
    tab.psfdepth=np.zeros(nrows)	
    tab.bassframe=np.ones(nrows, np.int16)	
    return tab
