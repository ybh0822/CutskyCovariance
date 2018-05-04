
from nbodykit.lab import *
from nbodykit import setup_logging

import os
import argparse

# where the data is
CSCRATCH = '/global/homes/y/ybh0822/BOSSDR12comb'

setup_logging()


# load the full catalog
fullcat = FITSCatalog('/global/homes/y/ybh0822/BOSSDR12comb/LOWZ/galaxy_DR12v5_LOWZ_South.fits')

# add the completeness weights
fullcat['WEIGHT'] = fullcat['WEIGHT_SYSTOT'] * (fullcat['WEIGHT_NOZ'] + fullcat['WEIGHT_CP'] - 1.0)

# select the right redshift range 
ZMIN = 0.16
ZMAX = 0.36
valid = (fullcat['Z'] > ZMIN)&(fullcat['Z'] < ZMAX) 

# select the valid objects 
subcat = fullcat[valid]

# total weight is FKP weight * completeness weight
total_weight =  subcat['WEIGHT']*subcat['WEIGHT_FKP']

# effective redshift
zeff = (total_weight*subcat['Z']).sum() / total_weight.sum()

# effective nbar
nbar =  (total_weight*subcat['NZ']).sum() / total_weight.sum()

# finally, compute
zeff = zeff.compute()
nbar = nbar.compute()

print(zeff)
print(nbar)
