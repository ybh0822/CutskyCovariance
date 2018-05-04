
from nbodykit.lab import *
from nbodykit import setup_logging
import os
import argparse

# where the data is
CSCRATCH = '/global/homes/y/ybh0822/BOSSDR12comb'

setup_logging()

# the DR12 redshift bins
ZBINS = [(0.2, 0.5), (0.4, 0.6), (0.5, 0.75)]

# the fiducial DR12 cosmology
cosmo = cosmology.Cosmology(H0=67.6, Om0=0.31, flat=True, Tcmb0=0.)

def read_randoms(sample, zbin):
    """
    Return a Source holding the DR12 combined sample randoms file for
    the specified sample and redshift bin
    """
    # load the catalog
    dirname = '%s/combined_sample/Randoms' %CSCRATCH
    path = os.path.join(dirname, 'random0_DR12v5_CMASSLOWZTOT_%s.fits' %sample)
    s = FITSCatalog(path, use_cache=True)

    # add the Position column
    s['Position'] = transform.SkyToCartesion(s['RA'], s['DEC'], s['Z'], cosmo, degrees=True)

    # randoms get a weight of unity
    s['WEIGHT'] = 1.0

    return s

def get_redshift_bin(source, zbin):

    # do the redshift bin selection
    zmin, zmax = ZBINS[zbin]
    valid = (source['Z'] > zmin)&(source['Z'] < zmax)
    return source[valid]

def window_paircount(ns):

    # the redshift bin
    zbin = ns.zbin

    # load the randoms
    randoms = read_randoms(ns.sample, zbin)
    randoms = randoms[::10]

    # trim to this redshift bin
    this_randoms = get_redshift_bin(randoms, zbin)

    # do the pair count
    redges = numpy.logspace(0, 4, 500)
    r = SurveyDataPairCount('2d', this_randoms, redges, cosmo, Nmu=100, ra='RA', dec='DEC', redshift='Z')

    # and save!
    output_dir =  "/global/homes/y/ybh0822/BOSSDR12comb/output"
    args = (ns.sample, ZBINS[zbin][0], ZBINS[zbin][1])
    output = os.path.join(output_dir, "more10_randoms_DDsmu_boss_dr12_combined_sample_%s_%.1f-%.1f" % args)
    output += '.json'
    r.save(output)

if __name__ == '__main__':

    desc = 'compute the pair count of the eBOSS randoms file'
    parser = argparse.ArgumentParser(description=desc)

    h = 'the redshift bin'
    parser.add_argument('zbin', type=int, choices=list(range(len(ZBINS))), help=h)

    h = 'the sample, either North or South'
    parser.add_argument('sample', type=str, choices=['North', 'South'], help=h)

    ns = parser.parse_args()
    window_paircount(ns)
