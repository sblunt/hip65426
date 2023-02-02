import numpy as np
import os
import matplotlib.pyplot as plt

from orbitize import system, read_input, priors, sampler
from orbitize import results

"""
Fits to run:
1. lit astrometry only (True, False, False)
2. lit astrometry +1 GRAVITY epoch (True, True, False)
3. lit astrometry +2 GRAVITY epochs (True, True, True)
4. GRAVITY astrometry only (False, True, True)
5. accepted fit with Kelly's prior
"""


"""
Begin keywords <<
"""
run_fit = True

lit_astrom = False
first_grav = True
second_grav = True

savedir = 'results/'

if not os.path.exists(savedir):
    os.mkdir(savedir)

if lit_astrom:
    savedir += 'with_literature_astrom'
if first_grav:
    savedir += 'with_first_vlti_point'
if second_grav:
    savedir += 'with_second_vlti_point'
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

HIP654_system = system.System(
    num_secondary_bodies, data_table, m_st, plx, fit_secondary_mass=False, mass_err=mass_err, 
    plx_err=plx_err
)

# Check that orbitizie! initialized everything correctly.
# (I wrote the code, therefore I do not trust the code.)
assert not HIP654_system.fit_secondary_mass
assert not HIP654_system.track_planet_perturbs


if run_fit:

    # run MCMC
    num_threads = 20
    num_temps = 20
    num_walkers = 1000
    num_steps = 1000000 # n_walkers x n_steps_per_walker
    burn_steps = 10000
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
fig = HIP654_results.plot_corner()
plt.savefig('{}/corner.png'.format(savedir), dpi=250)

# make orbit plot
fig = HIP654_results.plot_orbits()
plt.savefig('{}/orbit_IAD{}.png'.format(savedir), dpi=250)