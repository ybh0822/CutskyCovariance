from nbodykit.lab import *
import numpy as np
from pyRSD.rsdfit.data import PowerMeasurements

CSCRATCH = '/Users/Byeonghee/Desktop/PowerSpectrum/BOSS_DR12_combined/Cov_run_PATHCY/N_bin3_1000_newrun'
RESULTD = '/Users/Byeonghee/Desktop/PowerSpectrum/BOSS_DR12_combined/Cov_run_PATHCY/N_bin3_1000_newrun/2018-05-02_nlopt_30i_chain0__1.npz'

sample = 'N'
#zbin = 1, 2 or 3
zbin = 3
mock_num = 1000
rm_fb = 0
no_cp = 1
Run = 2


from pyRSD.rsd import GalaxySpectrum
from pyRSD.rsdfit import FittingDriver
d = FittingDriver.from_directory('/Users/Byeonghee/Desktop/PowerSpectrum/BOSS_DR12_combined/Cov_run_PATHCY/N_bin3_1000_newrun', model_file='Run1_Model_bin3.npy', results_file=RESULTD)
d.set_fit_results()
model = d.theory.model
# save the initialized model to disk
model.to_npy('Run2_Model_bin3.npy')




import numpy
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from pyRSD.rsd import GalaxySpectrum
from pyRSD.rsdfit.data import PoleCovarianceMatrix

if sample == 'N':
	region = 'North'
elif sample == 'S':
	region = 'South'

# load n(z) from file and interpolate it
filename = 'nbar_DR12v5_CMASSLOWZ_%s_om0p31_Pfkp10000.dat' %(region)
nbar = numpy.loadtxt(filename, skiprows=3)
nbar = spline(nbar[:,0], nbar[:,3])

# the sky fraction the survey covers
if sample == 'N':
	fsky = 0.143599315554
elif sample == 'S':
	fsky = 0.0610296309343


Y = numpy.loadtxt('krange_%s_%d.dat' %(sample, zbin))
k = Y[:,0]

# # the k values to compute covariance at
# k = numpy.arange(0., 0.4, 0.005) + 0.005/2

# the multipoles to compute covariance of
ells = [0,2,4]

# the redshift range of the survey
ZBINS = [(0.2, 0.5), (0.4, 0.6), (0.5, 0.75)]
zmin, zmax = ZBINS[zbin-1]

# the FKP weight P0
P0_FKP = 1e4

# the PoleCovarianceMatrix holding the Gaussian covariance
C = PoleCovarianceMatrix.cutsky_gaussian_covariance(model, k, ells, nbar, fsky, zmin, zmax, FKP_P0=P0_FKP)

C.to_plaintext('Cov_Run%d_%s_zbin%d_mock%d_rmfb%d_nocp%d.dat' %(Run, sample, zbin, mock_num, rm_fb, no_cp))


if zbin == 1:
	zeff = 0.383929
elif zbin == 2:
	zeff = 0.510029
elif zbin == 3:
	zeff = 0.60598

if zbin == 1:
	nbar = 0.000311824
elif zbin == 2:
	nbar = 0.000330795
elif zbin == 3:
	nbar = 0.000222788


DataFileIn = open("%s/params.dat" %CSCRATCH, "r")

DataList = DataFileIn.readlines()

Testoutput = open("%s/params.dat" %CSCRATCH, "w")

DataList[5] = "driver.init_from = 'result'\n"
DataList[11] = "driver.start_from = '%s'\n" %(RESULTD)
DataList[18] = "data.covariance = '%s/Cov_Run%d_%s_zbin%d_mock%d_rmfb%d_nocp%d.dat'\n" %(CSCRATCH, Run, sample, zbin, mock_num, rm_fb, no_cp)
DataList[21] = "data.data_file = '%s/pole_%s_zbin%d_mock%d_rmfb%d_nocp%d.dat'\n" %(CSCRATCH ,sample, zbin, mock_num, rm_fb, no_cp)
DataList[29] = "data.window_file = '%s/Window_%s_%d.dat'\n" %(CSCRATCH, sample, zbin)
DataList[64] = "theory.nbar = {'name': 'nbar', 'vary': False, 'value': %.10f, 'analytic': False, 'fiducial': %.10f}\n" %(nbar, nbar)
DataList[106] = 'model.z = %.6f\n' %(zeff)

Testoutput.writelines(DataList)
Testoutput.close()

print('rsdfit nlopt -m Run2_Model_bin%d.npy -p %s/params.dat -i 400 -o %s' %(zbin, CSCRATCH, CSCRATCH))

