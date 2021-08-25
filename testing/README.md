# Testing
Instruction for using the testing infrastructure. 
## Dependencies
Before running the testing script, if you are on sapling, please make sure you have netCDF4 for python installed. 
If not, you can run the following commands to first install pip on sapling: 
```bash
curl [https://bootstrap.pypa.io/get-pip.py -o get-pip.py](https://bootstrap.pypa.io/get-pip.py%20-o%20get-pip.py)
python get-pip.py
export PATH=$HOME/.local/bin:$PATH
```
Then, run the following command to install netCDF4 for python:
```bash
pip install netCDF4
```
## Run Test
The default test case is JW baroclinic wave on the 480km mesh. The Fortran outputs are already stored in the testing/fortran folder and need not be modified.
Before you run the testing script, please modify the following paths in testing/testing.py, changing the paths to your home directory on sapling:
```python
FORTRAN_FILE_PATH = '/home/luyuan/mpas-regent/testing/fortran/'
REGENT_FILE_PATH = '/home/luyuan/mpas-regent/testing/regent/'
```
Then, to run the test, simply run the following command in the mpas-regent folder: 
```bash
python testing/testing.py
```
There are two flags: 
**\-\-cell** indicates that the test is only going to be run on the cell region variables (pressure, pressure_p, surface_pressure, theta, and rho). 
**\-\-edge** indicates that the test is only going to be run on the cell region variable (v). 