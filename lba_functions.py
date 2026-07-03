import logging
import copy
from astropy.table import Table
from astroquery.utils.tap.core import TapPlus
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
import numpy as np
from pathlib import Path

def get_lba_catalogues(cat_path):
    '''Find catalogues in the given path. Also set the CASDA url for later use.

    Parameters
    -----------
    cat_path: str
        path to directory containing catalogues to be read in.
    '''

    # Find all catalogues in the input directory
    catalogues = [
            str(file) for file in Path(cat_path).rglob("*") if
            file.is_file()]

    lba_catalogues = {}
    emu_cat = []
    vlass_cat = []
    casdatap = []
    # only want to read in the catalogues once if we can
    for catname in catalogues:
        if 'EMU' in catname:
            emu_cat = Table.read(catname)
            lba_catalogues['emu'] = emu_cat
            logging.info(f'Read in EMU catalogue: {catname}')
        elif 'VLASS' in catname:
            vlass_cat = Table.read(catname, format='csv')
            lba_catalogues['vlass'] = vlass_cat
            logging.info(f'Read in VLASS catalogue: {catname}')
        else:
           logging.info(f"{catname} not recognised, skipping")
    # we also want to search RACS Mid
    casda_url = "https://casda.csiro.au/casda_vo_tools/tap"
    lba_catalogues['racs'] = casda_url
    logging.info(f"Setting up CASDA tap (for RACS Mid): {casda_url}")

    return lba_catalogues

def select_lba_catalogue(source_name, lba_catalogues, radius, pointing_centre):
    '''Select the catalogue that gives the most phase centres for the
    pointing_centre of this target source

    Parameters
    -----------
    source_name: str
        the name of the source in the Vex file
    lba_catalogues: list of catalogues
        the catalogues to be checked. The names/formats must be known to this routine.
    radius: float
        test radius within which to search for sources
    pointing_centre: list(float) of length 2
        the coordinates in degrees of the search centre

    Returns
    --------
    catalogue: chosen astro catalogue of sources giving the most matches
    RA_column, Dec_column, flux_colum: str
        the field name for the RA, Dec and flux column respectively in the
        chosen catalogue
    surv: str
        the name of the survey from which the chosen catalogue derives
    '''

    logging.info('Estimating the phase centres for source: %s'%source_name)
    #ras = []
    #decs = []
    #ra_center = Angle(vexfile.source[source_name]['ra']).degree
    #dec_center = Angle(vexfile.source[source_name]['dec']).degree
    ra_center, dec_center = pointing_centre
    #ras.append(ra_center)
    #decs.append(dec_center)
    #pointing_centre = [ra_center, dec_center]
    total_phase_centres = 0
    catalogue = None
    surv = None
    keep_cols = []
    #print(lba_catalogues.keys())
    if 'racs' in lba_catalogues.keys():
        casda_url = lba_catalogues['racs']
        casdatap = TapPlus(url=casda_url)
        job = casdatap.launch_job_async(
                "SELECT * FROM AS110.racs_mid_components_v01 where 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{},{},{}))".format(ra_center,dec_center,radius/60.))
        r = job.get_results()
        logging.info(f'Number of potential phase centers in RACS Mid: {len(r)}')
        if len(r) > total_phase_centres:
            keep_cols = ['ra','dec','total_flux']
            total_phase_centres = len(r)
            #catalogue='RACS_potential_phase_centers_{}.csv'.format(source_name)
            catalogue = r[keep_cols]
            surv = 'RACS'
    if 'emu' in lba_catalogues.keys():
        emu_cat = lba_catalogues['emu']
        center = SkyCoord(ra_center*u.deg, dec_center*u.deg)
        emu_coords = SkyCoord(
                ra=emu_cat['ra_deg_cont'], dec=emu_cat['dec_deg_cont'])
        sep = emu_coords.separation(center)
        emu_potential_phase = emu_cat[np.where((sep<=radius*u.arcmin))]
        logging.info(
                f'Number of potential phase centers in EMU: {len(emu_potential_phase)}')
        # EMU wins in a tie
        if len(emu_potential_phase) >= total_phase_centres:
            keep_cols = ['ra_deg_cont', 'dec_deg_cont', 'flux_int']
            catalogue = emu_potential_phase[keep_cols]
            total_phase_centres = len(emu_potential_phase)
            surv = 'EMU'
    if 'vlass' in lba_catalogues.keys():
        vlass_cat = lba_catalogues['vlass']
        vlass_coords = SkyCoord(
                ra=vlass_cat['RA']*u.deg, dec=vlass_cat['DEC']*u.deg)
        vlass_sep = vlass_coords.separation(center)
        vlass_potential_phase = vlass_cat[np.where((vlass_sep<=radius*u.arcmin))]
        logging.info(
                f'Number of potential phase centers in VLASS: {len(vlass_potential_phase)}')
        if len(vlass_potential_phase) > total_phase_centres:
            keep_cols = ['RA', 'DEC', 'Total_flux']
            catalogue = vlass_potential_phase[keep_cols]
            surv = 'VLASS'

    logging.info(
            f'Chosen catalogue {surv} has {len(catalogue)} potential phase centres'
            )
    #print(keep_cols)
    RA_column = keep_cols[0]
    Dec_column = keep_cols[1]
    flux_column = keep_cols[2]
    #catalogue =  catalogue[keep_cols]
    return(catalogue, RA_column, Dec_column, flux_column, surv)
