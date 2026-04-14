import os
import sys
import emcee

import numpy as np

from scipy import stats
from pyazr import azure2
from brick.azr import AZR
from multiprocess import Pool, current_process

# Restrict processes to one thread only
os.environ['OMP_NUM_THREADS'] = '1'

# Define the parameters prior distributions
priors = [ ]

# Minimization variables
nsteps   = 10000     # How many steps should each walker take?
nprocs   = 20        # How many Python processes do you want to allocate?

# We read the .azr file and set the external capture file to speed up the calculation
azr = azure2('15O.azr', nprocs=nprocs, port=60000)
ndatasets = azr.nsegments

# Get the parameters from AZURE2
params = azr.params
params_rwa = azr.params_rwa

# params_rwa contaisn the normalization factors for the datasets, so we need to concatenate it with the R-matrix parameters
params = np.concatenate( (params, np.ones(ndatasets)) )

# We'll read the data from the output file since it's already in the center-of-mass frame
y = azr.cross
yerr = azr.cross_err

# Calculate the lnL offset basing on the errors only
offset = 0 
for i in range(len(yerr)):
    offset += np.sum( np.log( 2 * np.pi * pow(yerr[i], 2) ) )

# We define variables
ndim     = len(params)  # How many parameters are you fitting?
nwalkers = 2 * ndim     # How many walkers do you want to use?

# Count how many R-matrix parameters
ntheta = len(params) - ndatasets

from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d
from scipy.stats import rv_continuous

lifetimbes = np.loadtxt('lifetime_samples.txt', usecols=0)
widths = 0.6582119514 / lifetimbes

# Clip widths to [0, 10]
widths_clipped = widths[(widths >= 0) & (widths <= 10)]

# Build KDE from the clipped data
kde = gaussian_kde(widths_clipped)

# Create a fine grid and normalize the PDF over [0, 20]
x_grid = np.linspace(0, 50, 10000)
pdf_vals = kde(x_grid)
pdf_vals = pdf_vals / np.trapz(pdf_vals, x_grid)  # ensure it integrates to 1

# Build CDF by cumulative integration
cdf_vals = np.cumsum(pdf_vals) * (x_grid[1] - x_grid[0])
cdf_vals = cdf_vals / cdf_vals[-1]  # normalize to exactly 1

# Create interpolators
pdf_interp = interp1d(x_grid, pdf_vals, bounds_error=False, fill_value=0.0)
cdf_interp = interp1d(x_grid, cdf_vals, bounds_error=False, fill_value=(0.0, 1.0))
ppf_interp = interp1d(cdf_vals, x_grid, bounds_error=False, fill_value=(0.0, 20.0))

# Define a custom scipy distribution
class WidthPrior(rv_continuous):
    def _pdf(self, x):
        return pdf_interp(x)
    def _cdf(self, x):
        return cdf_interp(x)
    def _ppf(self, q):
        return ppf_interp(q)

prior = WidthPrior(a=0, b=20, name='width_prior')

# Select the priors
priors = [ 

        # 0 MeV
        stats.norm(loc=7, scale=0.39), # ANC
        #stats.uniform(loc=-1e8, scale=2e8),

        # 7.5 MeV state
        stats.uniform(loc=-1e8, scale=2e8), # Proton Width
        stats.uniform(loc=-1e8, scale=2e8), # Gamma Width gs
        stats.uniform(loc=-1e8, scale=2e8), # Gamma Width 6.79 MeV

        # 8.7 MeV state
        stats.uniform(loc=-1e8, scale=2e8), # Proton Width

        # 9.6 MeV state
        stats.uniform(loc=-1e8, scale=2e8), # Proton Width
        stats.uniform(loc=-1e8, scale=2e8), # Gamma Width GS
        stats.uniform(loc=-1e8, scale=2e8), # Gamma Width 6.79 MeV

        # 6.79 MeV 
        stats.norm(loc=4.9, scale=0.51), # ANC
        #stats.uniform(loc=-1e8, scale=2e8),
        prior, # Width
        #stats.uniform(loc=-1e8, scale=2e8),

        # 8.28 MeV state
        stats.uniform(loc=-1e8, scale=2e8), # Proton Width 1
        stats.uniform(loc=-1e8, scale=2e8), # Proton Width 2
        stats.uniform(loc=-1e8, scale=2e8), # Proton Width 3
        stats.uniform(loc=-1e8, scale=2e8), # Gamma Width

        # 9.49 MeV state
        stats.uniform(loc=-1e8, scale=2e8), # Proton Width 1
        stats.uniform(loc=-1e8, scale=2e8), # Proton Width 2
        stats.uniform(loc=-1e8, scale=2e8), # Proton Width 3
        stats.uniform(loc=-1e8, scale=2e8), # Gamma Width

        # 3/2+ Background Pole
        stats.uniform(loc=-1e8, scale=2e8), # Gamma Width

        # 8.98 MeV state
        stats.uniform(loc=-1e8, scale=2e8), # Proton Width
        stats.uniform(loc=-1e8, scale=2e8), # Gamma Width

        # 9.48 MeV state
        stats.uniform(loc=-1e8, scale=2e8), # Proton Width 1
        stats.uniform(loc=-1e8, scale=2e8), # Proton Width 2
        stats.uniform(loc=-1e8, scale=2e8), # Proton Width 3
        stats.uniform(loc=-1e8, scale=2e8), # Gamma Width

        # 8.92 MeV state
        stats.uniform(loc=-1e8, scale=2e8), # Proton Width

        # Frentz uncertainty
        stats.halfcauchy(scale=0.07),

        # Marta uncertainty
        stats.halfcauchy(scale=0.08),

        # Chen uncertainty
        stats.halfcauchy(scale=0.08),

        # Li uncertainty
        stats.halfcauchy(scale=0.07),
        stats.halfcauchy(scale=0.07),
        stats.halfcauchy(scale=0.07),
        stats.halfcauchy(scale=0.07),
        stats.halfcauchy(scale=0.07),
        stats.halfcauchy(scale=0.07),

        # Runkle uncertainty
        stats.halfcauchy(scale=0.10),

        # Wagner uncertainty
        stats.halfcauchy(scale=0.08),

        # deBoer uncertainty
        stats.halfcauchy(scale=0.05),
        stats.halfcauchy(scale=0.05),
        stats.halfcauchy(scale=0.05),
        stats.halfcauchy(scale=0.05),
        stats.halfcauchy(scale=0.05),

        # Imbriani uncertainty
        stats.halfcauchy(scale=0.08),

        # Chen 6.79 MeV state
        stats.halfcauchy(scale=0.08),

        # Frentz 6.79 MeV state
        stats.halfcauchy(scale=0.07),

        # Imbriani 6.79 MeV state
        stats.halfcauchy(scale=0.08),

        # Wagner 6.79 MeV uncertainty
        stats.halfcauchy(scale=0.08),

        # Li 6.79 MeV state
        stats.halfcauchy(scale=0.07),
        stats.halfcauchy(scale=0.07),
        stats.halfcauchy(scale=0.07),
        stats.halfcauchy(scale=0.07),
        stats.halfcauchy(scale=0.07),
        stats.halfcauchy(scale=0.07),

        # Marta 6.79 MeV uncertainty
        stats.halfcauchy(scale=0.08),

        # Runkle 6.79 MeV uncertainty
        stats.halfcauchy(scale=0.10)

           ]


print(f"Fitting {ndim} parameters with {len(priors)} priors, {nwalkers} walkers and {nsteps} steps each on {nprocs} processes, with {ndatasets} datasets.")

# Prior log probability
def lnPi( theta ):
    classic_prior = np.sum( [pi.logpdf(t) for (pi, t) in zip(priors, theta)] )
    omega = 1 / 3
    omega_gamma = np.abs( theta[1] * theta[2] ) / (np.abs(theta[1]) + np.abs(theta[2])) * omega
    omega_gamma_prior = stats.norm(loc=12.9 * 0.0149 * 1e-3, scale=0.5 * 0.0149 * 1e-3).logpdf(omega_gamma)
    return classic_prior + omega_gamma_prior

# Log likelihood
def lnL( theta, proc=0 ):
    norms = theta[ntheta:]
    model = azr.calculate( theta, proc=proc )
    lnl = 0
    for i in range(len(y)):
        err_tot = np.sqrt(yerr[i]**2 + (model[i] * norms[i])**2) # Add in quadrature the systematic uncertainty
        lnl += -0.5 * np.sum( (y[i] - model[i])**2 / (2 * pow(err_tot, 2)) )
        lnl += -0.5 * np.sum(np.log(2 * np.pi * err_tot**2))
    return lnl

# Posterior log probability
def lnP( theta ):
    try: proc = int(current_process().name.split('-')[1]) - 1 # We want to get the numbe r of the process to call the right AZURE2 port
    except: proc = 0 
    lnpi = lnPi( theta )
    if np.any(np.isnan(lnpi)) or np.any(np.isinf(lnpi)): return -np.inf
    lnl = lnL( theta, proc=proc )
    if np.any(np.isnan(lnl)) or np.any(np.isinf(lnl)): return -np.inf
    return lnl + lnpi

# Prepare initial walker positions
p0 = np.zeros( (nwalkers, ndim) )
for i in range(nwalkers):
    for j in range(ndim):
        if( j >= ntheta ): # For the normalization factors, we start from 1
            p0[i, j] = stats.halfcauchy(scale=0.05).rvs() # Start from a random value drawn from the prior distribution
        else:
            p0[i, j] = np.random.normal(params[j], 0.05 * np.abs(params[j]))

# Prepare the file to write the chains
backend = emcee.backends.HDFBackend('results/samples-complete.h5') 

try:
    n_iterations_done = backend.iteration
    print(f"Restarting from iteration {n_iterations_done}")
    with Pool(processes=nprocs) as pool:
        sampler = emcee.EnsembleSampler( nwalkers, ndim, lnP, pool=pool, backend=backend )
        state = sampler.run_mcmc( None, nsteps, progress=True, tune=True )

except:
    print("Starting from scratch")
    with Pool(processes=nprocs) as pool:
        sampler = emcee.EnsembleSampler( nwalkers, ndim, lnP, pool=pool, backend=backend ) 
        state = sampler.run_mcmc( p0, nsteps, progress=True, tune=True )

# Kill all AZURE2 servers
os.system("killall AZURE2")

# Exit with code 0
sys.exit(0)