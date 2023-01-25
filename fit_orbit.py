import numpy as np
import os
import matplotlib.pyplot as plt

from orbitize import system, read_input, priors, sampler
from orbitize.hipparcos import HipparcosLogProb
from orbitize.gaia import GaiaLogProb
from orbitize import results

"""
Fits to run:
- lit astrometry only
- lit astrometry +1 GRAVITY epoch
- lit astrometry +2 GRAVITY epochs
- GRAVITY astrometry only
- accepted fit with Kelly's prior
"""


"""
Begin keywords <<
"""
fit_IAD = True 
gaia_dr = 'edr3'

if fit_IAD:
    savedir = 'plots/HG'
else:
    savedir = 'plots/noHG'
"""
>> End keywords
"""

# if not os.path.exists(savedir):
#     os.mkdir(savedir)

# input_file = 'HIP65426.csv'
# plx = 9.303088053872793 # [mas] (Gaia eDR3)
# plx_err = 0.034559656
# m_st = 1.96 # [M_sol] (Bowler+ 2020)

# num_secondary_bodies = 1
# data_table = read_input.read_file(input_file)

# if fit_IAD:
#     hipparcos_number='065426'
#     gaia_nums = {'edr3': 6070080754075553792, 'dr2':6070080754075553792}
#     fit_secondary_mass=True
#     hipparcos_filename='/data/user/sblunt/HipIAD-2021/ResRec_JavaTool_2014/H{}/H{}.d'.format(
#         hipparcos_number[0:3], hipparcos_number
#     )
#     HIP654_Hip = HipparcosLogProb(
#         hipparcos_filename, hipparcos_number, num_secondary_bodies
#     )
#     HIP654_gaia = GaiaLogProb(
#         gaia_nums[gaia_dr], HIP654_Hip, dr=gaia_dr
#     )

# else:
#     fit_secondary_mass=False
#     HIP654_Hip = None
#     HIP654_gaia = None

# HIP654_system = system.System(
#     num_secondary_bodies, data_table, m_st, plx, hipparcos_IAD=HIP654_Hip, 
#     gaia=HIP654_gaia, fit_secondary_mass=fit_secondary_mass, mass_err=0.04, 
#     plx_err=plx_err
# )

# if fit_IAD:
#     assert HIP654_system.fit_secondary_mass
#     assert HIP654_system.track_planet_perturbs

#     HIP654_system.sys_priors[HIP654_system.param_idx['plx']] = priors.UniformPrior(plx - 3 * plx_err, plx + 3 * plx_err)


# else:
#     assert not HIP654_system.fit_secondary_mass
#     assert not HIP654_system.track_planet_perturbs

# # run MCMC
# num_threads = 100
# num_temps = 20
# num_walkers = 1000
# num_steps = 1000000 # n_walkers x n_steps_per_walker
# burn_steps = 10000
# thin = 100

# HIP654_sampler = sampler.MCMC(
#     HIP654_system, num_threads=num_threads, num_temps=num_temps, 
#     num_walkers=num_walkers
# )
# HIP654_sampler.run_sampler(num_steps, burn_steps=burn_steps, thin=thin)

# # save chains
# HIP654_sampler.results.save_results(
#     '{}/HIP654_IAD{}.hdf5'.format(savedir, fit_IAD)
# )

# make corner plot
HIP654_results = results.Results() # create blank results object for loading
HIP654_results.load_results('{}/HIP654_IAD{}.hdf5'.format(savedir, fit_IAD))
# fig = HIP654_results.plot_corner()
# plt.savefig('{}/corner_IAD{}.png'.format(savedir, fit_IAD), dpi=250)

# make orbit plot
fig = HIP654_results.plot_orbits()
plt.savefig('{}/orbit_IAD{}.png'.format(savedir, fit_IAD), dpi=250)