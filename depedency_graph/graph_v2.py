from graphviz import Digraph
from var_reader import populate_read_vars
from var_reader import populate_write_vars
from var_reader import populate_readwrite_vars
import re
import sys
import mesh_vars



#desired var is the variable we want inbound arrows to
def add_layer(desired_var, all_vars, dot, task_arr, seen_vars, end_vars):
    
    seen_vars.append(desired_var)

    #if we've not seen the variable before, add a node to the graph
    if desired_var not in all_vars:
        all_vars.append(desired_var)
        dot.node(desired_var, desired_var)

    end_var_bool = True

    #iterate over tasks
    for task in task_arr:
        #extract the task name
        task_name = re.search('^[^\(]+', task).group(0).strip()
        reads_vars_output = []
        writes_vars_output = []

        #call our helper functions to extract the variables that are read and written
        populate_write_vars(task, writes_vars_output)
        populate_readwrite_vars(task, reads_vars_output, writes_vars_output)

        #If this task writes to the variable we want, we want to add inbound edges. It also means the variable can't be an 'end variable'
        if (desired_var in writes_vars_output): 
            populate_read_vars(task, reads_vars_output)
            end_var_bool = False

            #Add edges between all of the variables that are read in the task and our desired variable
            for read_var in reads_vars_output:
                #Check that the variables aren't from the mesh, though
                if (read_var not in mesh_vars.mesh_vars_set):
                    #If we haven't seen the variable before, add a node
                    if read_var not in all_vars:
                        dot.node(read_var, read_var)
                        all_vars.append(read_var)
                    dot.edge(read_var, desired_var, URL=task_name, Tooltip = task_name)

                    #If we 
                    if (read_var not in seen_vars):
                        add_layer(read_var, all_vars, dot, task_arr, seen_vars, end_vars)

    if (end_var_bool): 
        end_vars.append(desired_var)
        print("Edge var is", desired_var)

#initialize our graph                
dot = Digraph(comment='MPAS', engine ='fdp')

#open our file to parse
textfile = open('input.rg', 'r')
filetext = textfile.read()
textfile.close()

#This creates an array of all the tasks
task_arr = re.findall('(?s)(?<=task)(.*?)(?=do\n)', filetext)

#All_vars are the variables that have been added to the graph
all_vars = []

#Seen_vars are the variables that have been recursed on
seen_vars = []

#End_Vars are the variables that are not written to
end_vars = []

#our variable in question is the argument passed to the script
start_var = sys.argv[1] 

#Call our first call of the function
add_layer(start_var, all_vars, dot, task_arr, seen_vars, end_vars)

u = dot.unflatten(stagger=5)
u.render('output/dynamics_tasks.gv', view=True)

