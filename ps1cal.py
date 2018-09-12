import sys
import os
import numpy as np
import healpy
from astrometry.util.fits import fits_table, merge_tables
from astrometry.libkd.spherematch import match_radec

def read_chunks(ra, dec, chunkdir = os.environ['PS1CHUNKS']):

	phi = np.radians(ra)
	theta = np.radians(90-dec)
	pixels = healpy.pixelfunc.ang2pix(32,theta,phi)
	pixels = np.unique(pixels)
	cat = []

	for pix in pixels:
		fn = chunkdir+'ps1-%05d.fits'%pix
		try: 
			cat.append( fits_table(fn, ext=1) )
		except ValueError: 
			continue
		except IOError:
			continue
	if not cat: 
		raise RuntimeError('No calibration stars in field')
	return merge_tables(cat)

def read_stars(ra, dec, filt='r', mlo=None, mhi=None):
	
        j = 'gr'.find(filt) # index into PS1 arrays

        cat = read_chunks(ra, dec)
	
        gmi = cat.median[:,0] - cat.median[:,2]

        I = (gmi>0.4) & (gmi<2.7) & \
            (cat.nmag_ok[:,0] >= 1)  & \
            (cat.nmag_ok[:,2] >= 1) 

        if mlo is not None:
                J = cat.median[:,j]>mlo
        else:
                J = np.ones(len(cat),dtype=bool)
        if mhi is not None:
                K = cat.median[:,j]<mhi
        else:
                K = np.ones(len(cat),dtype=bool)

        return cat[I & J & K]

def cross_match(cat, stars=True, rad=2.0, **kwargs):
        
        if stars:
            objs = read_stars(cat.alpha_j2000, cat.delta_j2000, 
                              filt=kwargs.get('filt','r'), 
                              mlo=kwargs.get('mlo'),
                              mhi=kwargs.get('mhi'))
        else:
            objs = read_chunks(cat.alpha_j2000, cat.delta_j2000)

        if len(objs) is 0:
            raise RuntimeError('No calibration stars in field')

        m1, m2, d12 = match_radec(objs.ra, objs.dec, cat.alpha_j2000, cat.delta_j2000, rad/3600.0, nearest=True)

        if not(m1.size):
            raise RuntimeError('No calibration stars in field')

        cat = cat[m2]; objs=objs[m1]
        cat.add_columns_from(objs)
        
        cat.angdist = d12*3600.
        return cat

def poly_color_transform(mtch, filt, coeff, ref = ['g', 'i']):
    
    r1, r2 = ref
    j = 'grizy'.find(filt)
    l = 'grizy'.find(r1)
    m = 'grizy'.find(r2)

    if j < 0: 
        raise ValueError('Unknown filter', filt)
    if l < 0: 
        raise ValueError('Unknown filter', r1)
    if m < 0: 
        raise ValueError('Unknown filter', r2)

    psmag = mtch.median[:,j]
    psref = mtch.median[:,l]-mtch.median[:,m]

    return psmag + np.polyval(coeff,psref)

