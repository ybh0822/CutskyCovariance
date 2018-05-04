from nbodykit.lab import *
from nbodykit import setup_logging

import os
import argparse

setup_logging()

# the DR12 redshift bins
ZBINS = [(0.2, 0.5), (0.4, 0.6), (0.5, 0.75)]

# Planck 2013 cosmology
h = 0.6777
Ob0 = 0.048
Ocdm0 = 0.307115 - Ob0
cosmo = cosmology.Cosmology(h=h, Omega_b=Ob0, Omega_cdm=Ocdm0, m_ncdm=None, T_cmb=2.7255)

def read_data(sample, mock_num, zbin, rm_fb_obj=False, no_cp_weights=False):
	
	# the directory 
	if sample == 'N':
		dr = 'CMASSLOWZTOT-NGC'
	else:
		dr = 'SGC'

	# the data file path
	filename = "Patchy-Mocks-DR12CMASSLOWZTOT-%s-V6S-Portsmouth-mass_%04d.dat" %(sample, mock_num)
	data_path = os.path.join("/global/project/projectdirs/boss/galaxy/PATCHY/PATCHY-V6S/CMASSLOWZ", dr)
	data_path = os.path.join(data_path, filename)

	# load the CSV file
	names = ['ra', 'dec', 'z', 'mstar', 'NZ', 'bias', 'veto_flag', 'cp_weight']
	d = CSVCatalog(data_path, names=names)
	
	# do the redshift bin selection
	zmin, zmax = ZBINS[zbin-1]
	valid = (d['z'] > zmin)&(d['z'] < zmax)
	d = d[valid]

	# Remove fiber-collided object
	if rm_fb_obj :
		valid = (d['cp_weight'] > 0)
		d = d[valid] # slice the source

	# add the Position column
	d['Position'] = transform.SkyToCartesian(d['ra'], d['dec'], d['z'], cosmo, degrees=True)

	# FKP weight
	d['WEIGHT_FKP'] = 1./(1. + 10000.*d['NZ'])

	# complete weight
	if not no_cp_weights :
		d['WEIGHT'] = d['veto_flag']*d['cp_weight']
	else:
		d['WEIGHT'] = d['veto_flag']

	return d

def read_randoms(sample, zbin, rm_fb_obj=False, no_cp_weights=False):

	# load the catalog
	randoms_path = "/global/project/projectdirs/boss/galaxy/PATCHY/PATCHY-V6S/CMASSLOWZ/RANDOMS/CMASSLOWZTOT/Random-DR12CMASSLOWZTOT-%s-V6S-x50.dat" %(sample)
	names = ['ra', 'dec', 'z', 'NZ', 'bias', 'veto_flag', 'cp_weight']
	r = CSVCatalog(randoms_path, names = names)

	# do the redshift bin selection
	zmin, zmax = ZBINS[zbin-1]
	valid = (r['z'] > zmin)&(r['z'] < zmax)
	r = r[valid]

	# Remove fiber-collided object
	if rm_fb_obj :
		valid = (r['cp_weight'] > 0)
		r = r[valid] # slice the source

	# add the Position column
	r['Position'] = transform.SkyToCartesian(r['ra'], r['dec'], r['z'], cosmo, degrees=True)

	# FKP weight
	r['WEIGHT_FKP'] = 1./(1. + 10000.*r['NZ'])

	# complete weight
	if not no_cp_weights :
		r['WEIGHT'] = r['veto_flag']*r['cp_weight']
	else:
		r['WEIGHT'] = r['veto_flag']

	return r

def compute_power(ns):

	# load the data first
	data = read_data(ns.sample, ns.mock_num, ns.zbin, rm_fb_obj=ns.rm_fb_obj, no_cp_weights=ns.no_cp_weights)

	# load the randoms
	randoms = read_randoms(ns.sample, ns.zbin, rm_fb_obj=ns.rm_fb_obj, no_cp_weights=ns.no_cp_weights)

	# combine data and randoms into the FKP source
	fkp = FKPCatalog(data=data, randoms=randoms, BoxPad=0.1)

	# specify painting params and relevant columns
	fkp = fkp.to_mesh(Nmesh=512, window='tsc', interlaced=True, dtype='f8', nbar='NZ', fkp_weight='WEIGHT_FKP', comp_weight='WEIGHT')

	# compute the power
	LMAX = 16
	poles = list(range(0, LMAX+1, 2))
	s = ConvolvedFFTPower(fkp, poles=poles, dk=0.005, kmin=0.)

	# and save!
	output_dir =  "/global/homes/y/ybh0822/BOSSDR12comb/PATCHYmock/output2_bin1"
	output = os.path.join(output_dir, "Poles_%s_mock%d_bin%d_dk005_kmin0" %(ns.sample, ns.mock_num, ns.zbin))
	if ns.no_cp_weights:
		output += '_no_cp_weights'
	if ns.rm_fb_obj:
		output += '_rm_fb_obj'
	output += '.json'
	s.save(output)

if __name__ == '__main__':

	desc = 'compute the power spectrum of the BOSS DR12 combined sample'
	parser = argparse.ArgumentParser(description=desc)

	h = 'the sample, either North or South'
	parser.add_argument('sample', type=str, choices=['N', 'S'], help=h)

	h = 'the mock sample number'
	parser.add_argument('mock_num', type=int, choices=list(range(1, 1001)), help=h)

	h = 'the redshift bin, one of 1,2,3'
	parser.add_argument('zbin', type=int, choices=[1,2,3], help=h)

	h = 'whether to include fiber-collided objects'
	parser.add_argument('--rm_fb_obj', action='store_true', help=h)

	h = 'whether to include close pair weights for fiber collisions'
	parser.add_argument('--no_cp_weights', action='store_true', help=h)

	ns = parser.parse_args()
	compute_power(ns)























