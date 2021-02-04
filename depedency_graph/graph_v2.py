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

    #print("Desired var is", desired_var)

    if desired_var not in all_vars:
        all_vars.append(desired_var)
        dot.node(desired_var, desired_var)


    end_var_bool = True
    for task in task_arr:
        task_name = re.search('^[^\(]+', task).group(0).strip()
        reads_vars_output = []
        writes_vars_output = []

        populate_write_vars(task, writes_vars_output)
        populate_readwrite_vars(task, reads_vars_output, writes_vars_output)


        if (desired_var in writes_vars_output): 
            populate_read_vars(task, reads_vars_output)
            end_var_bool = False

            
            for read_var in reads_vars_output:
                if (read_var not in mesh_vars.mesh_vars_set):
                    if read_var not in all_vars:
                        dot.node(read_var, read_var)
                        all_vars.append(read_var)


                    dot.edge(read_var, desired_var, URL=task_name, Tooltip = task_name)
                    if (read_var not in seen_vars):
                        add_layer(read_var, all_vars, dot, task_arr, seen_vars, end_vars)

    if (end_var_bool): 
        end_vars.append(desired_var)
        print("Edge var is", desired_var)

                
dot = Digraph(comment='MPAS', engine ='fdp')


textfile = open('input.rg', 'r')
filetext = textfile.read()
textfile.close()

task_arr = re.findall('(?s)(?<=task)(.*?)(?=do\n)', filetext)
all_vars = []
seen_vars = []
end_vars = []


start_var = sys.argv[1] 

add_layer(start_var, all_vars, dot, task_arr, seen_vars, end_vars)

u = dot.unflatten(stagger=5)
u.render('output/dynamics_tasks.gv', view=True)

