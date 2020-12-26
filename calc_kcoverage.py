#! /usr/bin/env python
'''
Calculates thermal noise levels of a 21cm experiment given uv coverage
data calculation with mk_array_file.py.
'''
import aipy as a, numpy as n, optparse, sys
from scipy import interpolate
import numpy as np, h5py
from math import sqrt
from matplotlib import pyplot as plt
from copy import deepcopy

o = optparse.OptionParser()
o.set_usage('calc_kcoverage.py [options] *.npz')
o.set_description(__doc__)
o.add_option('-m', '--model', dest='model', default='mod',
    help="The model of the foreground wedge to use.  Three options are 'pess' (all k modes inside horizon + buffer are excluded, and all baselines are added incoherently), 'mod' (all k modes inside horizon + buffer are excluded, but all baselines within a uv pixel are added coherently), and 'opt' (all modes k modes inside the primary field of view are excluded).  See Pober et al. 2014 for more details.")
o.add_option('-b', '--buff', dest='buff', default=0.1, type=float,
    help="The size of the additive buffer outside the horizon to exclude in the pessimistic and moderate models.")
o.add_option('--ndays', dest='ndays', default=180., type=float,
    help="The total number of days observed.  The default is 180, which is the maximum a particular R.A. can be observed in one year if one only observes at night.  The total observing time is ndays*n_per_day.")
o.add_option('--bwidth', dest='bwidth', default=0.008, type=float,
    help="Cosmological bandwidth in GHz.  Note this is not the total instrument bandwidth, but the redshift range that can be considered co-eval.  Default is 0.008 (8 MHz).")
opts, args = o.parse_args(sys.argv[1:])

#=========================COSMOLOGY/BINNING FUNCTIONS=========================

#Convert frequency (GHz) to redshift for 21cm line.
def f2z(fq):
    F21 = 1.42040575177
    return (F21 / fq - 1)

#Multiply by this to convert an angle on the sky to a transverse distance in Mpc/h at redshift z
def dL_dth(z):
    '''[h^-1 Mpc]/radian, from Furlanetto et al. (2006)'''
    return 1.9 * (1./a.const.arcmin) * ((1+z) / 10.)**.2

#Multiply by this to convert a bandwidth in GHz to a line of sight distance in Mpc/h at redshift z
def dL_df(z, omega_m=0.31):
    '''[h^-1 Mpc]/GHz, from Furlanetto et al. (2006)'''
    return (1.7 / 0.1) * ((1+z) / 10.)**.5 * (omega_m/0.15)**-0.5 * 1e3

#Multiply by this to convert a baseline length in wavelengths (at the frequency corresponding to redshift z) into a tranverse k mode in h/Mpc at redshift z
def dk_du(z):
    '''2pi * [h Mpc^-1] / [wavelengths], valid for u >> 1.'''
    return 2*n.pi / dL_dth(z) # from du = 1/dth, which derives from du = d(sin(th)) using the small-angle approx

#Multiply by this to convert eta (FT of freq.; in 1/GHz) to line of sight k mode in h/Mpc at redshift z
def dk_deta(z):
    '''2pi * [h Mpc^-1] / [GHz^-1]'''
    return 2*n.pi / dL_df(z)

#scalar conversion between observing and cosmological coordinates
def X2Y(z):
    '''[h^-3 Mpc^3] / [str * GHz]'''
    return dL_dth(z)**2 * dL_df(z)

#A function used for binning
def find_nearest(array,value):
    idx = (n.abs(array-value)).argmin()
    return idx


#====================OBSERVATION/COSMOLOGY PARAMETER VALUES====================

#Load in data from array file; see mk_array_file.py for definitions of the parameters
array = n.load(args[0])
name = array['name']
obs_duration = array['obs_duration']
dish_size_in_lambda = array['dish_size_in_lambda']
Trx = array['Trx']
t_int = array['t_int']
if opts.model == 'pess':
    uv_coverage = array['uv_coverage_pess']
else:
    uv_coverage = array['uv_coverage']

h = 0.68
B = opts.bwidth
z = f2z(array['freq'])

dish_size_in_lambda = dish_size_in_lambda*(array['freq']/.150) # linear frequency evolution, relative to 150 MHz
first_null = 1.22/dish_size_in_lambda #for an airy disk, even though beam model is Gaussian
bm = 1.13*(2.35*(0.45/dish_size_in_lambda))**2

Tsky = 60e3 * (3e8/(array['freq']*1e9))**2.55  # sky temperature in mK

#=================================MAIN CODE===================================

#set up blank arrays/dictionaries
kprs = []
#sense will include sample variance, Tsense will be Thermal only
sense, Tsense = {}, {}

uv_coverage *= t_int
SIZE = uv_coverage.shape[0]

#Before cutting out the half plane, make a deep copy of the uv coverage array
uv_covfull = deepcopy(uv_coverage)

#Compute coverage in k-space (k in 1/Mpc)
N = 128 #grid size
Nhalf = int(N/2)
L = 300 #Mpc
dk = 2*np.pi / L #spacing in frequency space (1/Mpc)

#Vectors of wavenumbers (in 1/Mpc)
modes = np.arange(N)
modes[modes > N/2] -= N
kx_vec = modes * dk
ky_vec = modes * dk
kz_vec = modes[:Nhalf+1] * dk

#The noise cubes in mK^2, with and without accounting for the horizon limit
noise_cube = np.zeros((N, N, Nhalf+1))
noise_cube_horizon = np.zeros((N, N, Nhalf+1)) #also accounting for the horizon
noise_cube_horizon_opt = np.zeros((N, N, Nhalf+1)) #with optimal removal

#The amount of coverage in k space, with and without horizon accounting
k_coverage = np.zeros((N, N, Nhalf+1))
k_coverage_horizon = np.zeros((N, N, Nhalf+1))
k_coverage_horizon_opt = np.zeros((N, N, Nhalf+1))

#Calculate beam term
bm2 = bm/2. #beam^2 term calculated for Gaussian; see Parsons et al. 2014
bm_eff = bm**2 / bm2

#System temperature
Tsys = Tsky + Trx

#Temperature squared, only needs to be divided by the uv_coverage
T2 = (X2Y(z) / h**3) * bm_eff * Tsys**2 / (2*(1e9) * opts.ndays)

#Run over all u,v modes
where = np.where(uv_covfull > 0)
for iu,iv in zip(where[1], where[0]):
    #Convert (u,v) to (kx,ky) in 1/Mpc
    u, v = (iu - SIZE/2) * dish_size_in_lambda, (iv - SIZE/2) * dish_size_in_lambda
    umag = n.sqrt(u**2 + v**2)
    kx = u * dk_du(z) * h
    ky = v * dk_du(z) * h

    #Calculate horizon limit for baseline of length umag
    if opts.model in ['mod','pess']: hor = dk_deta(z) * umag/array['freq'] + opts.buff
    elif opts.model in ['opt']: hor = dk_deta(z) * (umag/array['freq'])*n.sin(first_null/2)
    else: print '%s is not a valid foreground model; Aborting...' % opts.model; sys.exit()

    #The horizon limit assuming the optimal scenario
    hor_opt = dk_deta(z) * (umag/array['freq'])*n.sin(first_null/2)

    #Convert the horizons to 1/Mpc
    hor = hor * h
    hor_opt = hor_opt * h

    #Run over all frequency channels
    for kz in kz_vec:

        #Find the closest index in the k_cube to deposit the coverage
        kx_id = find_nearest(kx_vec, kx)
        ky_id = find_nearest(ky_vec, ky)
        kz_id = find_nearest(kz_vec, kz)

        #The amount of uv coverage (baseline intensity)
        cov = uv_covfull[iv,iu]

        #Deposit the coverage
        k_coverage[kx_id, ky_id, kz_id] += cov

        #Deposit the coverage if we are above the horizon limit
        if kz >= hor:
            k_coverage_horizon[kx_id, ky_id, kz_id] += cov

        #Deposit the coverage, assuming optimal foreground removal
        if kz >= hor_opt:
            k_coverage_horizon_opt[kx_id, ky_id, kz_id] += cov


#Coverage per kpl
total_cov_uv = uv_covfull.sum()
total_cov_k = k_coverage[:,:,0].sum()

print 'Total coverage (uv): %e' % total_cov_uv
print 'Total coverage (k): %e' % total_cov_k

#Where do we have nonzero coverage?
nonz = np.where(k_coverage > 0)
nonz_hor = np.where(k_coverage_horizon > 0)
nonz_hor_opt = np.where(k_coverage_horizon_opt > 0)

#Compute the noise by dividing T^2 by the coverage
noise_cube[nonz] = T2 / k_coverage[nonz]
noise_cube_horizon[nonz_hor] = T2 / k_coverage_horizon[nonz_hor]
noise_cube_horizon_opt[nonz_hor_opt] = T2 / k_coverage_horizon_opt[nonz_hor_opt]

#The redshift as a string
stylized_z = int(round(z))
stylized_z = str(stylized_z)

#Store the k-coverage array in HDF5 format
cube_fname = "hera_350_k_coverage_z" + stylized_z + ".h5"
f = h5py.File(cube_fname, mode="w")
f.create_group("Header")
f["Header"].attrs["Redshift"] = z
f["Header"].attrs["BoxLength"] = L
f["Header"].attrs["Temperature (mK)"] = sqrt(T2)
f["Noise"] = noise_cube
f["Noise_Horizon"] = noise_cube_horizon
f["Noise_Horizon_Optimal"] = noise_cube_horizon_opt
f.close()

print 'Exported noise cube to %s' %  cube_fname
