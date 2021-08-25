import netCDF4 as nc
from collections import defaultdict
import numpy as np
import sys
import subprocess
import os


# setting for default dict
def default():
    return 'No such variable.'


# define global variables
# is_same records if all variables are the same for the Fortran and the regent versions
# INPUT_VARS are variables read from command line arguments, specifying which variables are to be compared,
# if empty then compare all variables.
MPAS_DIR = '.'
FORTRAN_FILE_PATH = '/home/luyuan/mpas-regent/testing/fortran/'
REGENT_FILE_PATH = '/home/luyuan/mpas-regent/testing/regent/'
FORTRAN_VARS = defaultdict(default)
REGENT_VARS = defaultdict(default)
NOT_IN_FORTRAN = []
DIFF_VARS = [] 
SAME_VARS = []
INPUT_VARS = [] 
is_same = 1
tol = 0.005
final = True
total_loc = 37
num_timesteps = 15
ALIGNMENT = [1] * total_loc

LOCATION_DICT = {0: "mesh_loading.rg: load_mesh",
                 1: "init_atm_cases.rg: init_atm_case_jw",
                 2: "atm_core.rg: atm_core_init",
                 3: "atm_core.rg: atm_do_timestep: 0",
                 4: "atm_core.rg: atm_do_timestep: 1",
                 5: "atm_core.rg: atm_do_timestep: 2",
                 6: "atm_core.rg: atm_do_timestep: 3",
                 7: "atm_core.rg: atm_do_timestep: 4",
                 8: "atm_core.rg: atm_do_timestep: 5",
                 9: "atm_core.rg: atm_do_timestep: 6",
                 10: "atm_core.rg: atm_do_timestep: 7",
                 11: "atm_core.rg: atm_do_timestep: 8",
                 12: "atm_core.rg: atm_do_timestep: 9",
                 13: "atm_core.rg: atm_do_timestep: 10",
                 14: "atm_core.rg: atm_do_timestep: 11",
                 15: "atm_core.rg: atm_do_timestep: 12",
                 16: "atm_core.rg: atm_do_timestep: 13",
                 17: "atm_core.rg: atm_do_timestep: 14",
                 18: "dynamics_tasks.rg: atm_compute_output_diagnostics",
                 19: "mesh_loading.rg: write_output_plotting",
                 20: "dynamics_tasks.rg: atm_init_coupled_diagnostics",
                 21: "dynamics_tasks.rg: atm_compute_solve_diagnostics",
                 22: "dynamics_tasks.rg: mpas_reconstruct_2d",
                 23: "atmphys_init.rg: physics_init",
                 24: "dynamics_tasks.rg: atm_compute_mesh_scaling",
                 25: "dynamics_tasks.rg: atm_compute_damping_coefs",
                 26: "dynamics_tasks.rg: atm_rk_integration_setup",
                 27: "dynamics_tasks.rg: atm_compute_moist_coefficients",
                 28: "dynamics_tasks.rg: atm_compute_vert_imp_coefs",
                 29: "dynamics_tasks.rg: atm_compute_dyn_tend",
                 30: "dynamics_tasks.rg: atm_set_smlstep_pert_variables",
                 31: "dynamics_tasks.rg: atm_advance_acoustic_step",
                 32: "dynamics_tasks.rg: atm_divergence_damping_3d",
                 33: "dynamics_tasks.rg: atm_recover_large_step_variables",
                 34: "dynamics_tasks.rg: atm_compute_solve_diagnostics",
                 35: "dynamics_tasks.rg: atm_rk_dynamics_substep_finish",
                 36: "rk_timestep.rg: summarize_timestep"
                 }


num_timesteps = 15


def replace(original_code, new_code, place): 
    subprocess.call(["sed", "-i", "s*" + original_code + "*" + new_code + "*g", place])

def change_to_debug_mode():
    for i in range(20):
        if i == 3:
            replace("-- This is location " + str(i) + " in debug mode.",
                    "print_var_to_file(cell_region, edge_region, \"time_steps\", j + 1)", "main.rg")
            #replace("-- This is location " + str(i) + " in debug mode.",
            #        "print_var_to_file(cell_region, \"time_steps\", j / 8)", "main.rg")
            statements = ["if (not (j == 0)) then"]
            for i in range(num_timesteps): 
                statements.append("if (j == " + str(i+1) + ") then\\\nfile_name = \"testing/regent/regent_output_" + str(i + 3) + ".nc\"\\\nend")
            replace("-- temp location for modifying file names", "\\\n".join(statements)+"\\\nend", "mesh_loading/netcdf_tasks.rg")
        else:
            replace("-- This is location " + str(i) + " in debug mode.",
                    "print_var_to_file(cell_region, edge_region, \"testing/regent/regent_output_" + str(i) + ".nc\", 0)", "main.rg")
    for i in range (20, 26): 
        replace("-- This is location " + str(i) + " in debug mode.",
            "print_var_to_file(cr, er, \"testing/regent/regent_output_" + str(i) + ".nc\", 0)", "atm_core.rg")
    for i in range(26, total_loc): 
        if i == 26: 
            replace("-- This is location " + str(i) + " in debug mode.",
                    "var file = constants.c.fopen(\"testing/regent/regent_output_" + str(i) + ".nc\", \"r\")\\\nif (not isnull(file)) then\\\nconstants.c.fclose(file)\\\nelse\\\nprint_var_to_file(cr, er, \"testing/regent/regent_output_" + str(i) + ".nc\", 0)\\\nend", "dynamics/rk_timestep.rg")
        elif i == 33: 
            replace("-- This is location " + str(i) + " in debug mode.", 
                    "if (rk_step == 2) then\\\nfile = constants.c.fopen(\"testing/regent/regent_output_" + str(i) + ".nc\", \"r\")\\\nif (not isnull(file)) then\\\nconstants.c.fclose(file)\\\nelse\\\nprint_var_to_file(cr, er, \"testing/regent/regent_output_" + str(i) + ".nc\", 0)\\\nend\\\nend", "dynamics/rk_timestep.rg")
        else: 
            replace("-- This is location " + str(i) + " in debug mode.",
                    "file = constants.c.fopen(\"testing/regent/regent_output_" + str(i) + ".nc\", \"r\")\\\nif (not isnull(file)) then\\\nconstants.c.fclose(file)\\\nelse\\\nprint_var_to_file(cr, er, \"testing/regent/regent_output_" + str(i) + ".nc\", 0)\\\nend", "dynamics/rk_timestep.rg")


def change_to_normal_mode():
    for i in range(20):
        if i == 3:
            replace("print_var_to_file(cell_region, edge_region, \"time_steps\", j + 1)", 
                    "-- This is location " + str(i) + " in debug mode.", "main.rg")
            with open("mesh_loading/netcdf_tasks.rg", 'r') as f:
                linum = 0
                for line in f:
                    if "if (not (j == 0)) then" in line:
                        break
                    linum += 1
            replace("if (not (j == 0)) then", "-- temp location for modifying file names", "mesh_loading/netcdf_tasks.rg")
            os.system("sed -i \'" + str(linum + 2) + "," + str(linum + 2 + 3 * num_timesteps) + "d\' mesh_loading/netcdf_tasks.rg")
        else:
            replace("print_var_to_file(cell_region, edge_region, \"testing/regent/regent_output_" + str(i) + ".nc\", 0)",
                    "-- This is location " + str(i) + " in debug mode.", "main.rg")
    for i in range(20, 26): 
        replace("print_var_to_file(cr, er, \"testing/regent/regent_output_" + str(i) + ".nc\", 0)",
                "-- This is location " + str(i) + " in debug mode.", "atm_core.rg")
    for i in range(26, total_loc): 
        if i == 26: 
            with open("dynamics/rk_timestep.rg", 'r') as f:
                linum = 0
                for line in f:
                    linum += 1
                    if "var file = constants.c.fopen(\"testing/regent/regent_output_26.nc\", \"r\")" in line:
                        replace("var file = constants.c.fopen(\"testing/regent/regent_output_26.nc\", \"r\")", "-- This is location 26 in debug mode.", "dynamics/rk_timestep.rg")
                        os.system("sed -i \'" + str(linum + 1) + "," + str(linum + 5) + "d\' dynamics/rk_timestep.rg")
                        break

        else:
            with open("dynamics/rk_timestep.rg", 'r') as f:
                linum = 0
                for line in f:
                    linum += 1
                    if "file = constants.c.fopen(\"testing/regent/regent_output_" + str(i) + ".nc\", \"r\")" in line:
                        if i != 33: 
                            replace("file = constants.c.fopen(\"testing/regent/regent_output_" + str(i) + ".nc\", \"r\")", "-- This is location " + str(i) + " in debug mode.", "dynamics/rk_timestep.rg")
                            os.system("sed -i \'" + str(linum + 1) + "," + str(linum + 5) + "d\' dynamics/rk_timestep.rg")
                        else: 
                            replace("if (rk_step == 2) then", "-- This is location 33 in debug mode.", "dynamics/rk_timestep.rg")
                            os.system("sed -i \'" + str(linum) + "," + str(linum + 6) + "d\' dynamics/rk_timestep.rg")
                        break       

# function to parse input arguments
# if there are any entries following the file name, each should be a variable name
# add these variables to INPUT_VARS
def parse_input():
    global final
    global INPUT_VARS
    args = sys.argv
    if len(args) > 1:
        if "--cell" in sys.argv:
            INPUT_VARS = ["pressure_p", "pressure", "theta", "rho", "surface_pressure"]
        elif "--edge" in sys.argv:
            INPUT_VARS = ["v"]
        for arg in args[1:]:
            if arg[0] != "-":
                INPUT_VARS.append(arg)


# function to read in a nc file
# this function reads in variable names and values and put the values in different dicts
# according to specified version
def read_nc_file(file_path, version):
    dataset = nc.Dataset(file_path)
    for var in dataset.variables.values():
        if version == 'fortran':
            if final:
                if len(dataset[var.name].shape) < 3:
                    FORTRAN_VARS[var.name] = dataset[var.name][:]
                else:
                    FORTRAN_VARS[var.name] = dataset[var.name][-1, :, :]
            else:
                FORTRAN_VARS[var.name] = dataset[var.name][:].T
        else:
            REGENT_VARS[var.name] = dataset[var.name][:]
    # if version == 'regent':
    #     print(REGENT_VARS["pressure"][0][0])
    #     print(FORTRAN_VARS["pressure"][0][0])


# function to compare two values
# compare two values up to a specific tolerance threshold
# if not a float, the values should be exactly equal
def comp(v1, v2):
    if isinstance(v1, float):
        if v1 != 0.0:
            return abs((v1 - v2) / v1) <= tol
        elif v2 != 0.0:
            return abs((v1 - v2) / v2) <= tol
        else:
            return 1
    else:
        return v1 == v2


# function to compare two variables
# consider different situations in which the two variables may be different types
def compare_var(var_name, regent_var, fortran_var):
    global is_same

    # if neither is an array, compare the values directly
    if not np.ma.isarray(regent_var) and not np.ma.isarray(fortran_var):
        if not comp(regent_var, fortran_var):
            is_same = 0
            DIFF_VARS.append(var_name)
            print('[' + var_name + ']: different values for regent and fortran, regent: '
                  + regent_var + ' \nfortran: ' + fortran_var + '.')

    # if only one is an array, naturally they cannot be equal
    elif (np.ma.isarray(regent_var) and not np.ma.isarray(fortran_var)) or (not np.ma.isarray(regent_var)
                                                                            and np.ma.isarray(fortran_var)):
        is_same = 0
        DIFF_VARS.append(var_name)
        print('[' + var_name + ']: different types for regent and fortran, regent: '
              , regent_var, ' fortran: ', fortran_var)

    # if both are arrays, consider different sub-possibilities
    elif np.ma.isarray(regent_var) and np.ma.isarray(fortran_var):

        # if the two arrays have different shapes
        if regent_var.shape != fortran_var.shape:
            is_same = 0
            print('[' + var_name + ']: different shapes for regent and fortran, regent: '
                  + str(regent_var.shape) + ' fortran: ' + str(fortran_var.shape) + '.')
            if len(fortran_var.shape) == 3:
                compare_var(var_name, regent_var[:], fortran_var[-1, :, -1])
            if len(fortran_var.shape) == 2:
                if regent_var.shape[0] == fortran_var.shape[0]:
                    compare_var(var_name, regent_var[:], fortran_var[:, 0])
                elif regent_var.shape[0] == fortran_var.shape[1]:
                    compare_var(var_name, regent_var[:], fortran_var[0, :])
        else:
            # if both have empty shape, compare their values directly
            if regent_var.shape == ():
                if not comp(regent_var, fortran_var):
                    is_same = 0
                    DIFF_VARS.append(var_name)
                    print('[' + var_name + ']: different values for regent and fortran, regent: '
                          + regent_var + ' fortran: ' + fortran_var + '.')
            else:

                # iterate through the arrays to compare their values
                # change value of is_same to 0 if one of the values does not match up
                # to save time, we break the loop after finding one difference
                for index in range(regent_var.shape[0]):
                    same = 1
                    if np.ma.isarray(regent_var[index]):
                        if isinstance(regent_var[index], float):
                            if np.linalg.norm(regent_var[index] - fortran_var[index]) / \
                                    np.linalg.norm(regent_var[index]) > tol:
                                same = 0
                        elif not(regent_var[index] is np.ma.masked and fortran_var[index] is np.ma.masked):
                            for i in range(len(regent_var[index])):
                                if not comp(regent_var[index][i], fortran_var[index][i]) \
                                        and not (regent_var[index][i] is np.ma.masked and
                                                 fortran_var[index][i] is np.ma.masked):
                                    same = 0
                    else:
                        same = comp(regent_var[index], fortran_var[index])
                    if not same:
                        is_same = 0
                        DIFF_VARS.append(var_name)
                        print('[' + var_name + ']: different values for regent and fortran, regent: '
                              + str(regent_var[index][0]) + '  fortran: ' + str(fortran_var[index][0]) + ', at index ' + str(index) + '.')
                        break

    # record the variable name as different or same
    if var_name not in DIFF_VARS and var_name not in SAME_VARS:
        print('[' + var_name + ']: passed test.')
        SAME_VARS.append(var_name)


def comp_files(location: int):

    print("----------------------------------------------\nAt location", location, "in", LOCATION_DICT[location])
    # read both files and record variable values
    if location < 10:
        read_nc_file(FORTRAN_FILE_PATH + "fortran_output_0" + str(location) + ".nc", 'fortran')
    else:
        read_nc_file(FORTRAN_FILE_PATH + "fortran_output_" + str(location) + ".nc", 'fortran')
    read_nc_file(REGENT_FILE_PATH + "regent_output_" + str(location) + ".nc", 'regent')
    # print("Fortran pressure_p =", FORTRAN_VARS["pressure_p"][0][0])
    # print("Regent pressure_p =", REGENT_VARS["pressure_p"][0][0])
    # print("Fortran pressure =", FORTRAN_VARS["pressure"][0][0])
    # print("Regent pressure =", REGENT_VARS["pressure"][0][0])
    # print("Fortran surface_pressure =", FORTRAN_VARS["surface_pressure"][0])
    # print("Regent surface_pressure =", REGENT_VARS["surface_pressure"][0])

    # compare variables from regent and fortran
    if not INPUT_VARS:
        for var in REGENT_VARS.keys():
            compare_var(var, REGENT_VARS[var], FORTRAN_VARS[var])
    else:
        for var in INPUT_VARS:
            compare_var(var, REGENT_VARS[var], FORTRAN_VARS[var])

    if DIFF_VARS:
        print('The following variables from regent differ from fortran:', DIFF_VARS)
    if SAME_VARS:
        print('The following variables are the same for both versions:', SAME_VARS)

    if DIFF_VARS:
        REGENT_VARS.clear()
        FORTRAN_VARS.clear()
        DIFF_VARS.clear()
        SAME_VARS.clear()
        return location
    else:
        REGENT_VARS.clear()
        FORTRAN_VARS.clear()
        DIFF_VARS.clear()
        SAME_VARS.clear()
        return -1


def comp_final_output():

    # read both files and record variable values
    print('Reading both output files...')
    read_nc_file(FORTRAN_FILE_PATH + "480km_fortran_output.nc", 'fortran')
    read_nc_file(REGENT_FILE_PATH + "timestep_output.nc", 'regent')

    # print("Fortran pressure_p =", FORTRAN_VARS["pressure_p"][0][0])
    # print("Regent pressure_p =", REGENT_VARS["pressure_p"][0][0])
    # print("Fortran pressure =", FORTRAN_VARS["pressure"][0][0])
    # print("Regent pressure =", REGENT_VARS["pressure"][0][0])
    # print("Fortran surface_pressure =", FORTRAN_VARS["surface_pressure"][0])
    # print("Regent surface_pressure =", REGENT_VARS["surface_pressure"][0])

    # compare variables from regent and fortran
    print('Comparing values of each variable...')
    if not INPUT_VARS:
        for var in REGENT_VARS.keys():
            if var in FORTRAN_VARS.keys():
                compare_var(var, REGENT_VARS[var], FORTRAN_VARS[var])
            else:
                NOT_IN_FORTRAN.append(var)
    else:
        for var in INPUT_VARS:
            if var in FORTRAN_VARS.keys() and var in REGENT_VARS.keys():
                compare_var(var, REGENT_VARS[var], FORTRAN_VARS[var])
            else:
                NOT_IN_FORTRAN.append(var)

    # print summary of the overall result
    print('-------------------------------\nCompleted comparing variables.\nSummary:')
    if NOT_IN_FORTRAN:
        print('The following variables are not in the Fortran version:', NOT_IN_FORTRAN)

    if DIFF_VARS:
        print('The following variables from regent differ from fortran:', DIFF_VARS)

    if SAME_VARS:
        print('The following variables are the same for both versions:', SAME_VARS)

    if is_same:
        print('The regent output matches the fortran output.')


if __name__ == '__main__':

    change_to_debug_mode()
    r = os.system("/home/zengcs/regent main.rg -ll:csize 1024")
    change_to_normal_mode()

    if not r: 
        os.system("cp timestep_output.nc testing/regent")
        parse_input()
        final = False
        for i in range(20):
            loc_diff = comp_files(i)
            if loc_diff >= 0:
                print("Regent vars differ from fortran at location " + str(loc_diff))
                ALIGNMENT[loc_diff] = 0
                break
        if not ALIGNMENT[2]:
            for i in range(20, 26):
                loc_diff = comp_files(i)
                if loc_diff >= 0:
                    print("Regent vars differ from fortran at location " + str(loc_diff))
                    ALIGNMENT[loc_diff] = 0
                    break
        if not ALIGNMENT[3]:
            for i in range(26, total_loc):
                loc_diff = comp_files(i)
                if loc_diff >= 0:
                    print("Regent vars differ from fortran at location " + str(loc_diff))
                    ALIGNMENT[loc_diff] = 0
                    break
        print("***************************************************************")
        if not (0 in ALIGNMENT):
            final = True
            comp_final_output()
        else:
            print("Error occured in location", ALIGNMENT.index(0), "in", LOCATION_DICT[ALIGNMENT.index(0)])
            if ALIGNMENT.index(0) == 2:
                print("More specifically, error is in location", ALIGNMENT[20:].index(0) + 20, "in",
                    LOCATION_DICT[ALIGNMENT[20:].index(0) + 20])
            if ALIGNMENT.index(0) == 3:
                print("More specifically, error is in location", ALIGNMENT[20:].index(0) + 20, "in",
                    LOCATION_DICT[ALIGNMENT[20:].index(0) + 20])
        for i in range(26, total_loc): 
            os.system("rm -f testing/regent/regent_output_" + str(i) + ".nc")
    else:
        print("***************************************************************\nThe model did not run successfully, testing process not initiated.")
