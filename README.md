# mpas-regent

## Resources to learn about MPAS

I would start by going to https://mpas-dev.github.io/ and poking around, especially the links about the atmospheric model.

Then, you can skim through the tutorial here: http://www2.mmm.ucar.edu/projects/mpas/tutorial/Boulder2019/index.html <br /> 
File 4 in particular is very helpful for understanding the mesh structure.

I have uploaded a Google Drive folder with a bunch of PDFs I found helpful to understand things. The tutorial PDFs are also located there.
The link is https://drive.google.com/drive/folders/1d3mViA53ELeKhiph5kzJndGQwXw7zL_W?usp=sharing. <br />

The user guide is also very helpful (direct link: http://www2.mmm.ucar.edu/projects/mpas/mpas_atmosphere_users_guide_7.0.pdf). It is also in the google drive. Pages 65-69 are particularly useful for understanding the Voronoi mesh, as well as making sense of the variables in the mesh files. It might also provide some motivation for why we designed the data structures the way we did.

When you get to the stage where you want to understand the MPAS codebase, we have a google doc here: https://docs.google.com/document/d/1yF4sEZyL1xUkHHfx-uYcRANpZyboRIZo3iiaydksnlw/edit



## Instructions to run regent on Sherlock:


First, get Prof. Aiken to invite you to sherlock.

You can then log onto sherlock by doing ssh <sunetID>@login.sherlock.stanford.edu <br />
, and then typing in your Stanford password and 2FA. 

**load modules** <br />
module load python <br />
module load openmpi/2.0.2 <br />
module load netcdf <br />

**clone legion repo** <br />
git clone -b master https://github.com/StanfordLegion/legion.git <br />
cd legion/language <br />



**launch SLURM job** <br />
salloc --partition=aaiken --tasks 1 --nodes=1 --cpus-per-task=20 --time=05:00:00 <br />
srun --pty bash <br />

**untar terra and llvm builds** <br />
wget sapling.stanford.edu/~eslaught/terra.build.tar.gz <br />
tar -xzf terra.build.tar.gz <br />
wget sapling.stanford.edu/~eslaught/llvm.tar.gz <br />
tar -xzf llvm.tar.gz <br />

**setup** <br />
CC=gcc CXX=g++ CONDUIT=ibv ./scripts/setup_env.py <br />


It should be good to run now: Regent is not added to the path by default, so when running regent scripts, you have to invoke the regent.py file directly:

LAUNCHER="srun" ~/legion/language/regent.py <file_name>.rg <br />


In the future, when you login to Sherlock, you have to do the following:

module load python <br />
module load openmpi/2.0.2 <br />
module load netcdf <br />

salloc --partition=aaiken --tasks 1 --nodes=1 --cpus-per-task=20 --time=02:00:00 <br />

LAUNCHER="srun" ~/legion/language/regent.py <file_name>.rg <br />

## Running GPU Nodes
Elliott has put up a Regent build with CUDA here on Sherlock: /home/groups/aaiken/eslaught/regent_build_cuda_2020-09-03/language <br />
Thus, all we need to do is:

module load cuda <br />
cd /home/groups/aaiken/eslaught/regent_build_cuda_2020-09-03/language <br />
source env.sh <br />

<navigate to the folder with the main.rg regent file i want to run>
  
salloc --partition=aaiken --tasks 1 --nodes=1 --cpus-per-task=10 --gres=gpu:4 --time=02:00:00 <br />
LAUNCHER="srun" /home/groups/aaiken/eslaught/regent_build_cuda_2020-09-03/language/regent.py main.rg <br />


## Helpful tricks for working in Sherlock
You can avoid having to 2FA multiple times when logging into Sherlock by following the instructions here: 
https://www.sherlock.stanford.edu/docs/advanced-topics/connection/#avoiding-multiple-duo-prompts <br />


You can also edit files in Sherlock using VSCODE by doing something similar to this: 
https://www.youtube.com/watch?v=vpK4rXLc0WY&feature=youtu.be&ab_channel=RyanEberhardt <br />


## Installing MPAS (This still does not work :( )

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

Unfortunately however, I still have not been able to get MPAS to work.  I have an ongoing thread on the MPAS forum and they are trying to help me troubleshoot there. The thread is https://forum.mmm.ucar.edu/phpBB3/viewtopic.php?f=12&t=9462.

## Running regent-mpas
In the top level directory, run LAUNCHER="srun" ~/legion/language/regent.py main.rg. <br />
 
Please also add the following to your ~/.bash_profile so that terra knows where to look for the files we "require".<br />

export TERRA_PATH="$HOME/regent_project_2020/mpas-regent/mesh_loading/?.rg;$HOME/regent_project_2020/mpas-regent/dynamics/?.rg;$HOME/regent_project_2020/mpas-regent/?.rg;$HOME/regent_project_2020/mpas-regent/vertical_init/?.rg" <br />

## Overview of project:

### Mesh loading
If you navigate to MPAS-Atmosphere -> MPAS-Atmosphere meshes, you will find some MPAS meshes. If you download them, you will see that they have a grid.nc file. This contains the mesh data in netcdf format. If you have netcdf installed (you can load it easily on Sherlock), you can manipulate these files and see their contents easy using ncdump. Syntax for ncdump can be found here: http://www.bic.mni.mcgill.ca/users/sean/Docs/netcdf/guide.txn_79.html#:~:text=The%20ncdump%20tool%20generates%20an,variable%20data%20in%20the%20file.&text=Thus%20ncdump%20and%20ncgen%20can,between%20binary%20and%20ASCII%20representations.

For example, you can do ncdump x1.2562.grid.nc >> output.txt to dump the contents of the grid file into a txt file called output.txt, and ncdump -v "latCell" to dump the contents of variable latCell.

We have converted these netcdf files into regent data structures. This was possible because there is a C library to manipulate netcdf files (I think this is the one: https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html). Regent has support for calling C functions, so that is basically what we do to read the grid files into the data structures in mesh_loading.rg. Many of the C functions are wrapped in terra functions in netcdf_tasks.rg


### Partitioning
In the grid folder that you download, there is a .graph.info file. The graph.info files can be partitioned using a software called METIS.

If you'd like to partition the cells into N partitions, you run gpmetis graph.info N, which creates a file graph.info.part.N. 

This file has nCells rows, and each row has a number from 0-(N-1), which I assume to mean the partition that each cell is split into. 

I have a task called read_file in mesh_loading.rg that parses this graph.info file and returns an array where each element is the partition number of that cell index. We then assign this partition number to the cell and partition in regent based on this.

We also create halo regions around each partition.
partition_s_1 is the immediate ring around each partition, i.e. 10 neighbour cells
partition_s_2 is two rings around each partition, so it has the neighbours of all neighbours. I.e. 100 cells.

partition_halo_1 is just the inner halo, so partition_s_1 - cell_partition_initial
and partition_halo_2 is the outer halo, so partition_s2  - partition_halo_1 - cell_partition_initial

To understand the halo code, I would recommend reading the sections about images and preimages at http://regent-lang.org/reference/.  (And ask any questions that may come up, it's a little confusing).

You can also read about dependent partitioning here: https://drive.google.com/drive/u/1/folders/1d3mViA53ELeKhiph5kzJndGQwXw7zL_W

A significant TODO would be to find a better way to do dependent partitioning/haloes - because we require each of the 10^2 neighbours to be values in the cell region for dependent partitioning, our cell data structure has 100 fields which is very messy...
### 


## File by file overview

### main.rg
This is the overview files that calls all of our sub-tasks.

### data_structures.rg
This is the file in which we define our regions. We currently have 4 main regions - a cell region, an edge region, a vertex region, and a vertical region. We find that most variables either are parameterized by a cell, edge, or vertex, while a small number are only a property of the vertical level we are at (hence the need for the fourth region).

### constants.rg
We declare constants here for use in other file.

### mesh_loading/mesh_loading.rg
This file has the task to load the data from the netcdf grid file. 

task load_mesh: This loads the data from the netcdf grid file. You can change the grid name in constants.rg.

task partition_regions: This partitions the regions. It works correctly, however I still have not been able to return the partition objects.

task write_output: This writes an output file to test that we have read the file correctly. To be called just after load_mesh if verification is required.

task write_output_plotting: This writes an output file to test that we have read the file correctly. To be called after running timestep, etc.

### vertical_init.rg/init_atm_cases.rg
This initializes many variables needed before we run the RK timestep (including much of the vertical grid). We currently use the Jablonowski and Williamson baroclinic wave test case.

### mesh_loading/netcdf_tasks.rg
This is a helper file employed by mesh_loading.rg that has the terra wrapper tasks around the C/netcdf functions we use to read the grid file.

### dynamics/dynamics_tasks.rg
This is where the meat of our kernels are. 

### dynamics/rk_timestep.rg
This is the file that contains the logic for taking a time step.


## Plotting
We turn the mesh into a 'patch' using mpas_patches.py. 

The file mpas_patches.py is a helper script that is used to create a MatPlotLib ‘patch collection’ of each of the individual grid cells.

We then use matplotlib to plot the 'patch' object. Edits to 

The script can be run by doing: python mpas-plotting.py <netcdf output file> -v <variable>: e.g. python mpas_plot_pressure.py timestep_output.nc -v pressure_p.
  
The plot headers can be edited in mpas_plot_pressure.py - in general there is no need to touch mpas_patches.py.

More information on plotting with MPAS can be found at https://github.com/MiCurry/MPAS-Plotting (although this is not very comprehensive either).


