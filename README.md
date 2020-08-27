# mpas-regent





Plotting:
We turn the mesh into a 'patch' using mpas_patches.py. 

The file mpas_patches.py is a helper script that is used to create a MatPlotLib ‘patch collection’ of each of the individual grid cells.

We then use matplotlib to plot the 'patch' object. Edits to 

The script can be run by doing: python mpas-plotting.py <netcdf output file> -v <variable>: e.g. python mpas_plot_pressure.py timestep_output.nc -v pressure_p.
  
The plot headers can be edited in mpas_plot_pressure.py - in general there is no need to touch mpas_patches.py.

More information on plotting with MPAS can be found at https://github.com/MiCurry/MPAS-Plotting (although this is not very comprehensive either).
