from graphviz import Digraph
from var_reader import populate_read_vars
from var_reader import populate_write_vars
from var_reader import populate_readwrite_vars
import re
import sys

#Initialize our graph
dot = Digraph(comment='MPAS', engine ='fdp')

#Open our file to parse
textfile = open('input.rg', 'r')
filetext = textfile.read()
textfile.close()

#Get all the tasks
task_arr = re.findall('(?s)(?<=task)(.*?)(?=do\n)', filetext)

#All vars are the variables we've added to the graph
all_vars = []

#Check whether we want to build the graph of all the variables or just a subset
all_vars_bool = True if sys.argv[1] == "all" else False

#If just a subset, extract the variables we want and append to a list. Also create those nodes in our graph
if (not all_vars_bool): 
    desired_vars = [str(arg) for arg in sys.argv[1:]]

    for var in desired_vars:
        dot.node(var, var)
        all_vars.append(var)

#Loop over every task
for task in task_arr:
    
    #Extract the task name
    task_name = re.search('^[^\(]+', task).group(0).strip()

    #Data structures to store variables that are written to and read by task
    reads_vars_output = []
    writes_vars_output = []
 
    #Call our helper functions to populate read and write vars
    populate_read_vars(task, reads_vars_output)

    populate_write_vars(task, writes_vars_output)

    populate_readwrite_vars(task, reads_vars_output, writes_vars_output)

    print("reads vars:", reads_vars_output)
    print("writes vars:", writes_vars_output)

    #Loop over variables
    for reads_var in reads_vars_output:
        for write_var in writes_vars_output:
            #If the variable we are writing to is a variable we want (or if we want all variables)
            if all_vars_bool or (write_var in desired_vars): 
                #add tha variable node to the graph if it doesn't already exist
                if write_var not in all_vars:
                    dot.node(write_var, write_var)
                    all_vars.append(write_var)
                #same for the read variable
                if reads_var not in all_vars:
                    dot.node(reads_var, reads_var)
                    all_vars.append(reads_var)
                #and also add an edge between the read and write var
                dot.edge(reads_var, write_var, URL=task_name, Tooltip = task_name)

u = dot.unflatten(stagger=5)
u.render('output/dynamics_tasks.gv', view=True)

