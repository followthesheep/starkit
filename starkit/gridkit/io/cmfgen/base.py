import re
import os

from glob import glob

from astropy.io import fits
import pandas as pd
import numpy as np
from progressbar import ProgressBar



cmfgen_bibtex = """
@ARTICLE{CMFGEN REFERENCE HERE
}

"""

cmfgen_meta = {'bibtex':cmfgen_bibtex,
             'parameters':['teff', 'logg', 'mh', 'micro','heh'],
             'wavelength_unit':'Angstrom',
             'wavelength_type':'vacuum',
             'flux_unit': 'erg/s/cm^2/angstrom'}


def make_raw_index(mh=-0.08,alpha=-0.05,res=300000.0,spectra_dir='spectra',
                   metadata='ogrid_metadata'):
    """
    Read all Phoenix files and generate a raw index with filename association.

    Returns
    -------
        bosz_index : pd.DataFrame
    """
    all_fnames = glob(os.path.join(spectra_dir,'*.tsv'))

    nfiles = len(all_fnames)
    mh_arr = np.zeros(nfiles)
    alpha_arr = np.zeros(nfiles)
    teff_arr = np.zeros(nfiles)
    logg_arr = np.zeros(nfiles)
    micro_arr = np.zeros(nfiles)
    #rot_arr = np.zeros(nfiles)
    res_arr = np.zeros(nfiles)
    he_h_arr = np.zeros(nfiles)
    #pattern = re.compile('a(mp|mm)(\d+)(cp|cm)(\d+)+(op|om)(\d+)t(\d+)g(\d+)v(\d+)modrt(\d+)b(\d+)')
    #pattern_dir = re.compile('metal_(.....)\/carbon_(.....)\/alpha_(.....)')
    tab = pd.read_csv(metadata,delim_whitespace=True,skiprows=23)

    for i in np.arange(nfiles):
        filename = all_fnames[i]
        p = os.path.split(filename)[1]
        base = p.split('_')
        search_str = '_'.join(base[0:3])
        print(search_str)
        indx = tab[tab['MODEL'].str.contains(search_str)].index[0]
        print(indx)
        logg = tab['Logg'].iloc[indx]
        micro = tab['V(km/s)'].iloc[indx]
        he_h = tab['He/H'].iloc[indx]
        mh_arr[i] = float(mh)
        #ch_arr[i] = float(ch)
        alpha_arr[i] = float(alpha)
        teff_arr[i] = tab['Teff'].iloc[indx]
        logg_arr[i] = float(logg)
        micro_arr[i] = float(micro)
        he_h_arr[i] = float(he_h)
        #rot_arr[i] = float(rot)
        res_arr[i] = float(res)

    return pd.DataFrame({'mh':mh_arr,'alpha':alpha_arr,'teff':teff_arr,
                      'logg':logg_arr,'micro':micro_arr,
                      'res':res_arr,'heh':he_h_arr,'filename':all_fnames})

def make_grid_info(fname,spectra_dir='spectra',metadata='ogrid_metadata'):
    """
    Make the HDF5 Grid Info file

    Parameters
    ----------
    fname: str

    """

    raw_index = make_raw_index(spectra_dir=spectra_dir,metadata=metadata)
    wtab = pd.read_csv(raw_index.loc[0, 'filename'], delim_whitespace=True)
    wavelength = wtab['wavelength']

    with pd.HDFStore(fname) as fh:
        fh['index'] = raw_index
        fh['wavelength'] = pd.DataFrame(wavelength)
        fh['meta'] = pd.Series(cmfgen_meta)

def convert_bz2_memmap(fname):
    """
    Convert a bz2 file to memmap
    Parameters
    ----------
    fname : str

    Returns
    -------

    """
    fname_npy = fname.replace('.tsv', '.v1.npy')
    if os.path.exists(fname_npy):
        pass
    else:
        flux = pd.read_csv(fname, usecols=(1,), delim_whitespace=True, dtype=np.float64)
        flux = flux.values[:,0]
        np.save(fname_npy, flux)

def cache_cmfgen_grid(delete=False,spectra_dir='spectra'):
    """
    Extract and cache BOSZ grid
    Parameters
    ----------
    delete: bool
        will delete the existing file

    Returns
    -------

    """
    all_fnames = glob(os.path.join(spectra_dir,'*.tsv'))
    bar = ProgressBar(maxval=len(all_fnames))
    for i, fname in bar(enumerate(all_fnames)):
        convert_bz2_memmap(fname)
