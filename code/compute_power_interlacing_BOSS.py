
from nbodykit.lab import *
from nbodykit import setup_logging

import os
import argparse

# where the data is
CSCRATCH = '/global/homes/y/ybh0822/BOSSDR12comb'

setup_logging()

# the DR12 redshift bins
ZBINS = [(0.2, 0.5), (0.4, 0.6), (0.5, 0.75)]

# fiducial boss cosmology
h = 0.676
Ob0 = 0.022/h**2
Ocdm0 = 0.31 - Ob0
cosmo = cosmology.Cosmology(h=h, Omega_b=Ob0, Omega_cdm=Ocdm0, m_ncdm=None, T_cmb=2.7255)

def read_randoms(sample, zbin):
	"""
	Return a Source holding the DR12 combined sample randoms file for
	the specified sample and redshift bin
	"""
	# load the catalog
	dirname = '%s/combined_sample/Randoms' %CSCRATCH
	path = os.path.join(dirname, 'random0_DR12v5_CMASSLOWZTOT_%s.fits' %sample)
	s = FITSCatalog(path)
	
	# do the redshift bin selection
	zmin, zmax = ZBINS[zbin-1]
	valid = (s['Z'] > zmin)&(s['Z'] < zmax)
	s = s[valid]

	# add the Position column
	s['Position'] = transform.SkyToCartesian(s['RA'], s['DEC'], s['Z'], cosmo, degrees=True)

	# randoms get a weight of unity
	s['WEIGHT'] = 1.0

	return s

def read_data(sample, zbin, no_cp_weights=False):
	"""
	Return a Source holding the DR12 combined sample data file for
	the specified sample and redshift bin
	"""
	# load the catalog
	dirname = '%s/combined_sample/Data' %CSCRATCH
	path = os.path.join(dirname, 'galaxy_DR12v5_CMASSLOWZTOT_%s.fits' %sample)
	s = FITSCatalog(path)

	# do the redshift bin selection
	zmin, zmax = ZBINS[zbin-1]
	valid = (s['Z'] > zmin)&(s['Z'] < zmax)
	s = s[valid]

	# add the Position column
	s['Position'] = transform.SkyToCartesian(s['RA'], s['DEC'], s['Z'], cosmo, degrees=True)

	# add the systematic + cp weights (see eq 48 of Reid et al 2016)
	if not no_cp_weights :
		s['WEIGHT'] = s['WEIGHT_SYSTOT'] * (s['WEIGHT_NOZ'] + s['WEIGHT_CP'] - 1.0)
	else:
		s['WEIGHT'] = s['WEIGHT_SYSTOT'] * (s['WEIGHT_NOZ'] + 1.0 - 1.0)

	return s


def compute_power(ns):

	# load the data first
	data = read_data(ns.sample, ns.zbin, no_cp_weights=ns.no_cp_weights)

	# load the randoms
	randoms = read_randoms(ns.sample, ns.zbin)

	# combine data and randoms into the FKP source
	fkp = FKPCatalog(data=data, randoms=randoms, BoxPad=0.1)

	# specify painting params and relevant columns
	fkp = fkp.to_mesh(Nmesh=512, window='tsc', interlaced=True, dtype='f8', nbar='NZ', fkp_weight='WEIGHT_FKP', comp_weight='WEIGHT')

	# compute the power
	LMAX = 16
	poles = list(range(0, LMAX+1, 2))
	r = ConvolvedFFTPower(fkp, poles=poles, dk=0.005, kmin=0.)

	# and save!
	output_dir =  "/global/homes/y/ybh0822/BOSSDR12comb/Combined_output"
	output = os.path.join(output_dir, "CMASSLOWZTOT_poles_%d_%s_dk005_kmin0" %(ns.zbin, ns.sample))
	if ns.no_cp_weights:
		output += '_no_cp_weights'
	output += '.json'
	r.save(output)

if __name__ == '__main__':

	desc = 'compute the power spectrum of the BOSS DR12 combined sample'
	parser = argparse.ArgumentParser(description=desc)

	h = 'the sample, either North or South'
	parser.add_argument('sample', type=str, choices=['North', 'South'], help=h)

	h = 'the redshift bin, one of 1,2,3'
	parser.add_argument('zbin', type=int, choices=[1,2,3], help=h)

	h = 'whether to include close pair weights for fiber collisions'
	parser.add_argument('--no_cp_weights', action='store_true', help=h)

	ns = parser.parse_args()
	compute_power(ns)
