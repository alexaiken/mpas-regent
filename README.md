# mpas-regent





Plotting:
We turn the mesh into a 'patch' using mpas_patches.py. 

The file mpas_patches.py is a helper script that is used to create a MatPlotLib ‘patch collection’ of each of the individual grid cells.

We then use matplotlib to plot the 'patch' object. Edits to 

The script can be run by doing: python mpas-plotting.py <netcdf output file> -v <variable>: e.g. python mpas_plot_pressure.py timestep_output.nc -v pressure_p.
  
The plot headers can be edited in mpas_plot_pressure.py - in general there is no need to touch mpas_patches.py.

More information on plotting with MPAS can be found at https://github.com/MiCurry/MPAS-Plotting (although this is not very comprehensive either).



#Instructions to run regent on Sherlock:

**load modules**
module load python
module load openmpi/2.0.2
module load netcdf

**clone legion repo**
git clone -b master https://github.com/StanfordLegion/legion.git
cd legion/language

https://github.com/StanfordLegion/legion.git

**launch SLURM job**
salloc --partition=aaiken --tasks 1 --nodes=1 --cpus-per-task=20 --time=05:00:00
srun --pty bash

**untar terra and llvm builds**
wget sapling.stanford.edu/~eslaught/terra.build.tar.gz
tar -xzf terra.build.tar.gz
wget sapling.stanford.edu/~eslaught/llvm.tar.gz
tar -xzf llvm.tar.gz

**setup**
CC=gcc CXX=g++ CONDUIT=ibv ./scripts/setup_env.py


It should be good to run now: Regent is not added to the path by default, so when running regent scripts, you have to invoke the regent.py file directly:

LAUNCHER="srun" ~/legion/language/regent.py <file_name>.rg


In the future, when you login to Sherlock, you have to do the following:

module load python
module load openmpi/2.0.2
module load netcdf

salloc --partition=aaiken --tasks 1 --nodes=1 --cpus-per-task=20 --time=02:00:00

LAUNCHER="srun" ~/legion/language/regent.py <file_name>.rg
