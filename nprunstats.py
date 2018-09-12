#!/usr/bin/env python

import os
import sys
import glob
import numpy as np

import multiprocessing
from functools import partial

import fitsio

from scipy.stats import sigmaclip, tmean
import ps1cal

from astrometry.util.starutil import hms2ra, dms2dec
from astrometry.util.fits import fits_table, merge_tables

from npstats import (NinetyPrimeDepthStats,
                     init_stats_table)

from argparse import RawTextHelpFormatter
import argparse


def zeromean(arr):

    arr = np.array(arr)
    if np.sum(arr) == 0:
        return np.zeros(1)

    arr[arr == 0] = np.nan
    mean = np.nanmean(arr)
    return np.array([mean])

def get_key_val(hdr, key):
    
    '''
    header keys are not always present,
    avoid repeated exception handling
    '''

    try:
        return hdr[key.upper()]
    except NameError:
        pass
    except ValueError:
        pass
    except KeyError:
        pass

sxdir = os.path.join(os.environ['CPCATS'],'bass')
cpdir = os.path.join(os.environ['BASSDATA'], 'BOK_CP')

def _main(imn, opt):
    
    photodir = opt.get('photodir', os.getcwd())
    nightdir = opt.get('nightdir', os.getcwd())

    statstab = init_stats_table(imn)

    # header keys to record
    keys = [
        'filter',
        'exptime',
        'jd',
        'airmass',
        'camtemp',
        'dewtemp']
    # OBJECT, RA, DEC, GAIN, RDNOIS are special cases see below
    
    # loop over images and collect/compute statistcis
    # for n, imn in enumerate(images[0:5]):

    npstats = NinetyPrimeDepthStats(imn, sxdir = photodir, cpdir = nightdir)
    
    # populate table with header stats
    try:
        phdr = npstats.img_hdus[0].read_header()
    except TypeError:
        # failed to open image, nothing more can be done
        return statstab

    # OBJECT and RA/DEC are special cases
    keyval = get_key_val(phdr, 'OBJECT')

    print(npstats.img_name, end=" ")
    if not keyval is None: 
        try:
            # object is tile number + pass
            tile = int(str(keyval).strip()[:-1])
            # is the frame from the bass survey
            if 0 < tile < 8136:
                statstab.set( 'object', np.int32([keyval]) )
                #return statstab
            else:
                # not a bass frame
                statstab.bassframe = np.zeros(1, dtype = np.int16) 
                return statstab

        except ValueError:
            # not a bass frame
            statstab.bassframe = np.zeros(1, dtype = np.int16) 
            return statstab

    keyval = get_key_val(phdr, 'RA')
    if not keyval is None:
        hh,hm,hs = phdr['RA'].split(':')
        ra = hms2ra(int(hh),int(hm),float(hs))
        statstab.ra = np.array([ra])

    keyval = get_key_val(phdr, 'DEC')
    if not keyval is None:
        dd,dm,ds = keyval.split(':')
        sign = -1 if dd[0] == '-' else 1
        dec = dms2dec(sign,int(dd.replace('-','')[0:]),int(dm),float(ds))
        statstab.dec = np.array([dec])

    for key in keys:
        keyval = get_key_val(phdr, key)
        if not keyval is None:
            statstab.set(key, np.array([keyval]))

    for key in ['utc-obs','ha']:
        keyval = get_key_val(phdr, key)
        if not keyval is None:
            statstab.set(key, np.array([keyval]))

    # GAIN and RDNOIS are special cases
    statstab.gain = [[[ get_key_val(phdr, 'GAIN'+str(i))
                        if not(keyval is None) else 0.0 for i in amps ] for amps in npstats.ampNum ]]

    statstab.readnoise = [[[ get_key_val(phdr, 'RDNOIS'+str(i))
                             if not(keyval is None) else 0.0 for i in amps ] for amps in npstats.ampNum ]]

    filt = npstats.filtname
    if filt is None:
        # can't do much without a filter
        return statstab

    # loop over ccds
    for ext in range(1,5):
        
        try:
            n_hdus = len(npstats.wht_hdus)
        except TypeError:
            n_hdus = 0
        
        if n_hdus == 5:
            # if there are missing hdus, don't assume those present are healthy
            invvar = npstats.get_invvar_map(ext).ravel()
            invvar, _, _ = sigmaclip(invvar[invvar > 0.], high = 3, low = 3)
            sigma = 1. / np.sqrt(np.median(invvar))
            statstab.sigmaccd[0][ext-1] = sigma
        else:
            # we assume there is at least an healthy image
            sigma = npstats.skynoise( np.stats.img_hdus[ext].read() )
            statstab.sigmaccd[0][ext-1] = sigma
        
        try:
            n_hdus = len(npstats.cat_hdus)
        except TypeError:
            n_hdus = 0

        if n_hdus == 5:
            try:
                # FIX ME -- some of the catalog extensions seem
                #           to be corrupted -- investigate!
                bokstds = npstats.get_sxcat_stds(ext)
            except OSError:
                continue

            if not len(bokstds):
                # no standards stars, nothing more can be done here
                continue
            fwhm_img, _, _ = sigmaclip(bokstds.fwhm_image, low = 3, high = 3)
            statstab.fwhmimgccd[0][ext-1] = np.median(fwhm_img)
        else:
            # something went wrong during sextracting
            continue
        
        zero = 0.0 # set null zeropiont

        try:
            # cross_match() raises a RuntimeError when there are no ps1 stndards
            bokstds = ps1cal.cross_match(bokstds, filt = filt)
            coeffs = npstats.coeffs[filt]
            trancol = ps1cal.poly_color_transform(bokstds, filt, coeffs)
            bokmags = - 2.5 * np.log10(bokstds.flux_psf)
            dmags, _, _ = sigmaclip(trancol - bokmags, low = 3, high = 3)
            zero = np.median(dmags)
            zerorms = np.std(dmags)
            
            statstab.magzptccd[0][ext-1] = zero
            statstab.magzptrmsccd[0][ext-1] = zerorms
            statstab.nstdsccd[0][ext-1] = np.int16(len(bokstds))
            
            cosdec = np.cos( np.radians(bokstds.dec) )
            dalpha = (bokstds.ra - bokstds.alpha_j2000) * cosdec
            ddelta = bokstds.dec - bokstds.delta_j2000
            dalpha, _, _ = sigmaclip(dalpha, low = 3, high = 3)
            ddelta, _, _ = sigmaclip(ddelta, low = 3, high = 3)
        
            statstab.raoffccd[0][ext-1] = np.median(dalpha) * 3600.
            statstab.decoffccd[0][ext-1] = np.median(ddelta) * 3600.
        except RuntimeError:
            # no standards stars
            pass
        
        if zero > 0.0:
            skyrms = sigma / 10.**((zero - 22.5) / 2.5)
            statstab.skynoiseccd[0][ext-1] = skyrms
            
            try:
                n_hdus = len(npstats.psf_hdus)
            except TypeError:
                n_hdus = 0

            if not(n_hdus == 5):
                # something went wrong during psfex run
                continue

            statstab.fwhmpsfccd[0][ext-1] = npstats.psf_fwhm(ext)

            statstab.gausspsfdepthccd[0][ext-1] = npstats.gauss_psf_depth(sigma, npstats.psf_fwhm(ext), zero)

            psf = npstats.psf(ext)
            h, w = npstats.img_hdus[ext].read().shape
            patch = psf.getPointSourcePatch(w//2, h//2).patch
            psfnorm = np.sqrt(np.sum(patch**2))
            statstab.psfdepthccd[0][ext-1] = npstats.psf_depth(sigma, psfnorm, zero)

    statstab.magzpt = zeromean(statstab.magzptccd)
    statstab.magzptrms = zeromean(statstab.magzptrmsccd)
    statstab.raoff = zeromean(statstab.raoffccd)
    statstab.decoff = zeromean(statstab.decoffccd)
    statstab.nstds = np.array([np.sum(statstab.nstdsccd)], dtype = np.int16)
    statstab.fwhmpsf = zeromean(statstab.fwhmpsfccd)
    statstab.fwhmimg = zeromean(statstab.fwhmimgccd)
    statstab.sigma = zeromean(statstab.sigmaccd)
    statstab.skynoise = zeromean(statstab.skynoiseccd)
    statstab.gausspsfdepth = zeromean(statstab.gausspsfdepthccd)
    statstab.psfdepth = zeromean(statstab.psfdepthccd)

    return statstab

if __name__ == "__main__":
    
    import argparse

    parser = argparse.ArgumentParser(description='Tabulate statistics from BASS imaging nights',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('night', type = str, help = 'Observing night, YYYYMMNN')
    parser.add_argument('ofn', type = str, help = 'Name of output file')
    parser.add_argument('--sxdir', type = str, default = sxdir, help = 'Path to photometry products')
    parser.add_argument('--cpdir', type = str, default = cpdir, help = 'Path to CP reduced imaging')
    parser.add_argument('--nproc', type = int, default = 1, help='muliprocessing, run n processes')

    args = parser.parse_args()
    
    photodir = os.path.join(args.sxdir, 'CP' + args.night)
    nightdir = os.path.join(args.cpdir, 'CP' + args.night)

    opt = {'photodir': photodir, 'nightdir': nightdir}

    images = glob.glob( os.path.join(nightdir, '*_ooi_*' ))
                   
    if not(images):
        print('*Error* No reduced images')
        sys.exit()

    if args.nproc > 1:
        pool = multiprocessing.Pool(processes = args.nproc)
        t = [ pool.apply(_main, args=(imn, opt,)) for imn in images ]
    else:
        t = [ _main(imn, opt) for imn in images ]

    t = merge_tables(t)
    # print(t.about())
    t.writeto(args.ofn)

