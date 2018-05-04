from nbodykit.lab import *
import numpy as np
from pyRSD.rsdfit.data import PowerMeasurements

CSCRATCH = '/Users/Byeonghee/Desktop/PowerSpectrum/BOSS_DR12_combined/Cov_run_PATHCY/N_bin3_1000_newrun'

sample = 'N'
#zbin = 1, 2 or 3
zbin = 3
mock_num = 1000
num = 1000
rm_fb = 0
no_cp = 1
Run = 1

Nell = 3
r_N1 = ConvolvedFFTPower.load('Poles_N_mock1_bin3_dk005_kmin0_no_cp_weights.json')
k_N1 = r_N1.poles['k']
test_pole = [np.zeros(k_N1.shape) for i in range(Nell)]

# num = 100
for i in range(1, num+1):
    r_N1 = ConvolvedFFTPower.load('Poles_N_mock%d_bin3_dk005_kmin0_no_cp_weights.json' %i)

    k_N1 = r_N1.poles['k']
    p0_N1 = r_N1.poles['power_0'].real - r_N1.attrs['shotnoise']
    p2_N1 = r_N1.poles['power_2'].real
    p4_N1 = r_N1.poles['power_4'].real
    # p6_N1 = r_N1.poles['power_6'].real
    # p8_N1 = r_N1.poles['power_8'].real
    # p10_N1 = r_N1.poles['power_10'].real
    # p12_N1 = r_N1.poles['power_12'].real
    # p14_N1 = r_N1.poles['power_14'].real
    # p16_N1 = r_N1.poles['power_16'].real

    Nk = len(k_N1)
    ells = [0, 2, 4]

    test_pole[0] += p0_N1
    test_pole[1] += p2_N1
    test_pole[2] += p4_N1
    # test_pole[3] += p6_N1

test_pole = np.array(test_pole)/num

data = np.empty((Nk,Nell), dtype=[('k', 'f8'), ('power', 'f8')])

for i, ell in enumerate(ells):
	data['k'][:,i] = k_N1[:]
	data['power'][:,i] = test_pole[i]

names = ['pole_0', 'pole_2', 'pole_4']

measurements = PowerMeasurements.from_array(names, data)

measurements.to_plaintext("pole_%s_zbin%d_mock%d_rmfb%d_nocp%d.dat" %(sample, zbin, mock_num, rm_fb, no_cp))

DataOut = np.column_stack((k_N1, k_N1))
np.savetxt('krange_%s_%d.dat' %(sample, zbin), DataOut)

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

# the model instance to compute P(k,mu)
model = GalaxySpectrum.from_npy('Run1_Model_bin3.npy')

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

DataList[5] = "driver.init_from = 'fiducial'\n"
DataList[11] = 'driver.start_from = None\n'
DataList[18] = "data.covariance = '%s/Cov_Run%d_%s_zbin%d_mock%d_rmfb%d_nocp%d.dat'\n" %(CSCRATCH, Run, sample, zbin, mock_num, rm_fb, no_cp)
DataList[21] = "data.data_file = '%s/pole_%s_zbin%d_mock%d_rmfb%d_nocp%d.dat'\n" %(CSCRATCH ,sample, zbin, mock_num, rm_fb, no_cp)
DataList[29] = "data.window_file = '%s/Window_%s_%d.dat'\n" %(CSCRATCH, sample, zbin)
DataList[64] = "theory.nbar = {'name': 'nbar', 'vary': False, 'value': %.10f, 'analytic': False, 'fiducial': %.10f}\n" %(nbar, nbar)
DataList[106] = 'model.z = %.6f\n' %(zeff)

Testoutput.writelines(DataList)
Testoutput.close()

print('rsdfit nlopt -m Run1_Model_bin%d.npy -p %s/params.dat -i 400 -o %s' %(zbin, CSCRATCH, CSCRATCH))

