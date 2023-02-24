import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.time import Time

from orbitize import system, read_input, priors, sampler
from orbitize import results

"""
Fits to run:
1. lit astrometry only (True, False, False, False, False)
2. lit astrometry + first GRAVITY epoch (True, True, False, False, False)
2.5. lit astrometry + second GRAVITY epoch (True, False, True, False, False)
3. lit astrometry +2 GRAVITY epochs (True, True, True, False, False)
4. GRAVITY astrometry only (False, True, True, False, False)
5. all astrometry with e fixed to 0 (True, True, True, True, False)
6. accepted fit (#3) with linearly decreasing prior on e (True, True, True, False, True)

Fits converged:
2
5

Fits need to be run:
1 (running)
2.5 (running)
3 (running)
4 (running)
6 (running)
"""


"""
Begin keywords <<
"""
run_fit = False

lit_astrom = True
first_grav = True
second_grav = True
fix_ecc = True
lin_ecc_prior = False

savedir = 'results/'

if not os.path.exists(savedir):
    os.mkdir(savedir)

if lit_astrom:
    savedir += 'with_literature_astrom'
if first_grav:
    savedir += 'with_first_vlti_point'
if second_grav:
    savedir += 'with_second_vlti_point'
if fix_ecc:
    savedir += '_fixed_ecc'
if lin_ecc_prior:
    savedir += '_linear_ecc'
"""
>> End keywords
"""

if not os.path.exists(savedir):
    os.mkdir(savedir)

input_file = 'HIP65426.csv'
plx = 9.303088053872793 # [mas] (Gaia eDR3)
plx_err = 0.034559656
m_st = 1.96 # [M_sol] (Bowler+ 2020)
mass_err = 0.04

num_secondary_bodies = 1
data_table = read_input.read_file(input_file)

if not lit_astrom:
    data_table = data_table[-3:]
if not first_grav:
    data_table = data_table[0:-1 & -1]
if not second_grav:
    data_table = data_table[:-1]

print(data_table)

HIP654_system = system.System(
    num_secondary_bodies, data_table, m_st, plx, fit_secondary_mass=False, mass_err=mass_err, 
    plx_err=plx_err
)

# fix eccentricity to 0
if fix_ecc:
    HIP654_system.sys_priors[HIP654_system.param_idx['ecc1']] = 0

# set a linearly decreasing prior on ecc
if lin_ecc_prior:
    HIP654_system.sys_priors[HIP654_system.param_idx['ecc1']] = priors.LinearPrior(-2.18, 2.01)

# Check that orbitizie! initialized everything correctly.
# (I wrote the code, therefore I do not trust the code.)
assert not HIP654_system.fit_secondary_mass
assert not HIP654_system.track_planet_perturbs


if run_fit:

    # run MCMC
    num_threads = 15
    num_temps = 20
    num_walkers = 1000
    num_steps = 50_000_000 # n_walkers x n_steps_per_walker
    burn_steps = 10_000
    thin = 100

    HIP654_sampler = sampler.MCMC(
        HIP654_system, num_threads=num_threads, num_temps=num_temps, 
        num_walkers=num_walkers
    )
    HIP654_sampler.run_sampler(num_steps, burn_steps=burn_steps, thin=thin)

    # save chains
    HIP654_sampler.results.save_results(
        '{}/chains.hdf5'.format(savedir)
    )

# make corner plot
HIP654_results = results.Results() # create blank results object for loading
HIP654_results.load_results('{}/chains.hdf5'.format(savedir))

if fix_ecc:
    param_list = ['sma1','inc1','aop1','pan1','tau1','mtot','plx']
else:
    param_list = None
# fig = HIP654_results.plot_corner(param_list=param_list)
# plt.savefig('{}/corner.png'.format(savedir), dpi=250)

# make orbit plot
fig = HIP654_results.plot_orbits(num_epochs_to_plot=500, plot_astrometry=False)
radec_ax, sep_ax, pa_ax, cbar_ax = fig.axes

sep, serr, pa, paerr = data_table['quant1'],  data_table['quant1_err'], data_table['quant2'], data_table['quant2_err']
epoch = Time(data_table['epoch'], format='mjd').decimalyear

sphere_mask = [0,1,2,3,6,7]
naco_mask = [4,5]
grav_mask = [9,10]
sep_ax.errorbar(epoch[sphere_mask], sep[sphere_mask], serr[sphere_mask], marker='^', ec='purple', facecolor='white',ls='', label='SPHERE')
sep_ax.errorbar(epoch[naco_mask], sep[naco_mask], serr[naco_mask], marker='s', ls='', color='purple', label='NACO')
sep_ax.legend()
pa_ax.errorbar(epoch[sphere_mask], pa[sphere_mask], paerr[sphere_mask], marker='^', ec='purple', facecolor='white',ls='')
pa_ax.errorbar(epoch[naco_mask], pa[naco_mask], paerr[naco_mask], marker='s', ls='', color='purple')
plt.savefig('{}/orbit.png'.format(savedir), dpi=250)