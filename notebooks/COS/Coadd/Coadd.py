# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
from IPython import get_ipython

# %%
# Wide set of imports
import numpy as np
import matplotlib.pyplot as plt
import os
# import pandas as pd
# import scipy as sp
# import seaborn as sb
# import scipy.signal as signal
# from astropy.io import fits
# from astropy import units as u
from astroquery.mast import Observations as Obs
from astropy.io import fits
from astropy.table import Table
import glob


# %%
import urllib
import tarfile

from sys import path as system_path
# system_path.append('/Users/nkerman/Software/quickScripts/')
import natpy

system_path.append('/Users/nkerman/Projects/ullyses_dp/high_level_science_products/high_level_science_products')
from coadd import COSSegmentList, STISSegmentList, FUSESegmentList, CCDSegmentList
from coadd import abut

# %% [markdown]
# The source is [VFTS 72 *AKA* BI 253](https://simbad.u-strasbg.fr/simbad/sim-id?Ident=BI+253&submit=submit+id)

# %%
!mkdir -p ./zipdata ./extracted_data ./processed_data
# %%
ull_data_dl_path = 'https://ullyses.stsci.edu/files/vfts72.tar.gz'
zip_data_path, headers = urllib.request.urlretrieve(ull_data_dl_path, './zipdata/vfts72')

print(f"Downloaded zipped ULLYSES data files to {zip_data_path}")
# %%
tar = tarfile.open(zip_data_path, "r:gz")
tar.extractall('./extracted_data/')
tar.close()

datadir = './extracted_data/vfts72/'
# %%
allexposures = glob.glob(datadir+'*x1d*')+glob.glob(datadir+'*vo*')


# %%
for efile in allexposures:
    prihdr = fits.getheader(efile)
    tab = Table.read(efile)
    print(prihdr['TELESCOP'],prihdr['INSTRUME'], "\t: ", len(tab), '\trows of data')


# # %%
# ! /Users/nkerman/miniconda3/envs/astroconda/bin/python /Users/nkerman/Projects/ullyses_dp/high_level_science_products/high_level_science_products/wrapper.py -i ./extracted_data/vfts72/ -o ./processed_data/vfts72/output1


# # %%
# Table.read('./NatTest/vfts72/output2/hlsp_ullyses_hst-fuse_fuse-cos-stis_vfts72_uv_dr2_sed.fits')


# %%
def simple_plot(filename, ext = 0, ls = '-'):
    tab = Table.read(filename)[ext]
    wvln, flux = tab['WAVELENGTH'], tab['FLUX']
    plt.xlabel("Wavelength [$\AA$]")
    plt.ylabel("Flux $[ergs/s/cm^2/\AA]$")
    plt.plot(wvln,flux, label = os.path.basename(filename), linestyle=ls)


# %%
# simple_plot('./NatTest/vfts72/output2/hlsp_ullyses_hst-fuse_fuse-cos-stis_vfts72_uv_dr2_sed.fits')

# %%
extract_dir = './extracted_data/vfts72/'
fplist = glob.glob(extract_dir+'*x1d*')
fflist = [natpy.fitsfile.fitsfile(fpath) for fpath in fplist]
# %%
[print(ff.vitalstats) for ff in fflist]
# %%
g130_ex = 'lda930fhq_x1d.fits'
g160_ex = 'lda930fdq_x1d.fits'
g130_exA = Table.read(extract_dir + g130_ex)[0]
g130_exB = Table.read(extract_dir + g130_ex)[1]
g160_exA = Table.read(extract_dir + g160_ex)[0]
g160_exB = Table.read(extract_dir + g160_ex)[1]
# %%
# Here is the data we want to eventually coadd and abut.
plt.plot(g130_exA['WAVELENGTH'], g130_exA['FLUX'], label = 'G130M/1291,FUVA')
plt.plot(g130_exB['WAVELENGTH'], g130_exB['FLUX'], label = 'G130M/1291,FUVB')
plt.plot(g160_exA['WAVELENGTH'], g160_exA['FLUX'], label = 'G160M/1589,FUVA')
plt.plot(g160_exB['WAVELENGTH'], g160_exB['FLUX'], label = 'G160M/1589,FUVB')

plt.legend()
# %%
# %%
prod_g130 = COSSegmentList('G130M', path = extract_dir)
prod_g130.create_output_wavelength_grid()
prod_g130.coadd()
prod_g130.target = prod_g130.ull_targname(alias_dir='/Users/nkerman/Projects/ullyses_dp/high_level_science_products/high_level_science_products/')
prod_g130.write('/Users/nkerman/Desktop/prod_g130.fits', overwrite=True)
# %%
prod_g160 = COSSegmentList('G160M', path = extract_dir)
prod_g160.create_output_wavelength_grid()
prod_g160.coadd()
prod_g160.target = prod_g160.ull_targname(alias_dir='/Users/nkerman/Projects/ullyses_dp/high_level_science_products/high_level_science_products/')
prod_g160.write('/Users/nkerman/Desktop/prod_g160.fits', overwrite=True)
# %%
abutted_prod = abut(prod_g130,prod_g160)
abutted_prod.write('/Users/nkerman/Desktop/g130g160.fits', overwrite=True)
# %%
for i,f in enumerate(glob.glob('/Users/nkerman/Desktop/*0.fits')):
    simple_plot(f, ls = '-::'[i])
plt.legend()
plt.savefig('delete.png', dpi = 300)

# %%
read_in_abutted = Table.read('/Users/nkerman/Desktop/g130g160.fits')[0]
# %%
wvln, sampling = read_in_abutted['WAVELENGTH'][1:],np.diff(read_in_abutted['WAVELENGTH'])
plt.plot(wvln, sampling)
plt.title("Sampling rate as fn of Å")

print(f"The transition wavelength auto-chosen was around {wvln[np.argmax(sampling)]} with sampling {sampling[np.argmax(sampling)]}")
# %%
abutted_prod = abut(prod_g130,prod_g160, transition_wavelength = 1412)
abutted_prod.write('/Users/nkerman/Desktop/g130g160_manual_tw.fits', overwrite=True)
read_in_abutted = Table.read('/Users/nkerman/Desktop/g130g160_manual_tw.fits')[0]
# %%
wvln, sampling = read_in_abutted['WAVELENGTH'][1:],np.diff(read_in_abutted['WAVELENGTH'])
plt.plot(wvln, sampling)
plt.title("Sampling rate as fn of Å")

print(f"The transition wavelength user-chosen was around {wvln[np.argmax(sampling)]} with sampling {sampling[np.argmax(sampling)]}")
# %%

# %%
