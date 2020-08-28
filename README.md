# mpas-regent


## Instructions to run regent on Sherlock:

**load modules** <br />
module load python
module load openmpi/2.0.2
module load netcdf

**clone legion repo** <br />
git clone -b master https://github.com/StanfordLegion/legion.git
cd legion/language

https://github.com/StanfordLegion/legion.git

**launch SLURM job** <br />
salloc --partition=aaiken --tasks 1 --nodes=1 --cpus-per-task=20 --time=05:00:00
srun --pty bash

**untar terra and llvm builds** <br />
wget sapling.stanford.edu/~eslaught/terra.build.tar.gz
tar -xzf terra.build.tar.gz
wget sapling.stanford.edu/~eslaught/llvm.tar.gz
tar -xzf llvm.tar.gz

**setup** <br />
CC=gcc CXX=g++ CONDUIT=ibv ./scripts/setup_env.py


It should be good to run now: Regent is not added to the path by default, so when running regent scripts, you have to invoke the regent.py file directly:

LAUNCHER="srun" ~/legion/language/regent.py <file_name>.rg


In the future, when you login to Sherlock, you have to do the following:

module load python
module load openmpi/2.0.2
module load netcdf

salloc --partition=aaiken --tasks 1 --nodes=1 --cpus-per-task=20 --time=02:00:00

LAUNCHER="srun" ~/legion/language/regent.py <file_name>.rg

## Installing MPAS

Step 1: **Environment Variables**  <br />

You’ll need to add the openmpi installation directory to your path, as I have done. In general, you can also modify the PNETCDF and PIO environment variables to be wherever you want to install them.

My ~/.bash_profile file looked like this:
PATH=$PATH:$HOME/.local/bin:$HOME/bin:$HOME/openmpi_install/bin  <br />
export PATH  <br />

export CC=gcc  <br />
export FC=gfortran  <br />
export F77=gfortran  <br />
export MPICC=mpicc  <br />
export MPIF90=mpif90  <br />
export MPIF77=mpif90  <br />

export OPENMPI=$HOME/openmpi_install  <br />
export PNETCDF=$HOME/pnetcdf_install  <br />
export MPIFC=mpif90  <br />
export PNETCDF_PATH=$PNETCDF  <br />
export PIO=$HOME/pio_install  <br />


￼

Step 2: **Other modules**  <br />
ml purge
ml libfabric/1.10.1 gcc/9.1.0



Step 3: **Install OpenMPI**  <br />
wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.4.tar.bz2
tar -xvf openmpi-4.0.4.tar.bz2
cd openmpi-4.0.4
./configure --prefix=$OPENMPI --with-pmi-libdir=/usr/lib64 --with-pmix=internal --with-libevent=internal --with-slurm --without-verbs
make 
make install

(You can change —prefix=$OPENMPI to be wherever you want openMPI to be installed to)
Be sure to add the openmpi install directory to your path, as I have done above.

Step 4: **Install PNETCDF**  <br />
Wget https://parallel-netcdf.github.io/Release/parallel-netcdf-1.8.1.tar.gz
tar -xvf parallel-netcdf-1.8.1.tar
cd parallel-netcdf-1.8.1
./configure --prefix=$PNETCDF --disable-cxx 
make 
make install

Wget https://parallel-netcdf.github.io/Release/pnetcdf-1.12.1.tar.gz
tar -xvf pnetcdf-1.12.1.tar.gz
cd pnetcdf-1.12.1 


Step 5: **Install PIO**  <br />
wget https://github.com/NCAR/ParallelIO/archive/pio1_7_1.tar.gz
tar -xvf ParallelIO-pio1_7_1.tar
cd ParallelIO-pio1_7_1/pio
./configure --prefix=$PIO --disable-netcdf --disable-mpiio 
make 
make install

Step 6: **Install MPAS**   <br />
git clone https://github.com/MPAS-Dev/MPAS-Model.git 
cd MPAS-Model
make gfortran CORE=init_atmosphere 



## Installing MPAS

## Plotting
We turn the mesh into a 'patch' using mpas_patches.py. 

The file mpas_patches.py is a helper script that is used to create a MatPlotLib ‘patch collection’ of each of the individual grid cells.

We then use matplotlib to plot the 'patch' object. Edits to 

The script can be run by doing: python mpas-plotting.py <netcdf output file> -v <variable>: e.g. python mpas_plot_pressure.py timestep_output.nc -v pressure_p.
  
The plot headers can be edited in mpas_plot_pressure.py - in general there is no need to touch mpas_patches.py.

More information on plotting with MPAS can be found at https://github.com/MiCurry/MPAS-Plotting (although this is not very comprehensive either).


