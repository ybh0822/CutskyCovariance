from pyRSD.rsd import GalaxySpectrum

#zbin = 1, 2 or 3
zbin = 1


if zbin == 1:
	zeff = 0.383929
elif zbin == 2:
	zeff = 0.510029
elif zbin == 3:
	zeff = 0.60598

model = GalaxySpectrum(z = zeff, params='boss_dr12_fidcosmo.ini')

model.initialize()

model.kmax = 1.

print(model.kmax)

# save the initialized model to disk
model.to_npy('Run1_Model_bin%d.npy' %(zbin))