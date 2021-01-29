from graphviz import Digraph
from var_reader import populate_read_vars
from var_reader import populate_write_vars
from var_reader import populate_readwrite_vars
import re
import sys

dot = Digraph(comment='MPAS', engine ='fdp')

textfile = open('input.rg', 'r')
filetext = textfile.read()
textfile.close()

task_arr = re.findall('(?s)(?<=task)(.*?)(?=do\n)', filetext)
all_vars = []

all_vars_bool = True if sys.argv[1] == "all" else False

if (not all_vars_bool ): 
    desired_vars = [str(arg) for arg in sys.argv[1:]]

    for var in desired_vars:
        dot.node(var, var)
        all_vars.append(var)

for task in task_arr:
    
    task_name = re.search('^[^\(]+', task).group(0).strip()
    reads_vars_output = []
    writes_vars_output = []
 
    populate_read_vars(task, reads_vars_output)

    populate_write_vars(task, writes_vars_output)

    populate_readwrite_vars(task, reads_vars_output, writes_vars_output)

    print("reads vars:", reads_vars_output)
    print("writes vars:", writes_vars_output)

    for reads_var in reads_vars_output:
        for write_var in writes_vars_output:
            if all_vars_bool or (write_var in desired_vars): 
                if write_var not in all_vars:
                    dot.node(write_var, write_var)
                    all_vars.append(write_var)
                if reads_var not in all_vars:
                    dot.node(reads_var, reads_var)
                    all_vars.append(reads_var)
                dot.edge(reads_var, write_var, URL=task_name, Tooltip = task_name)

    

u = dot.unflatten(stagger=5)
u.render('output/dynamics_tasks.gv', view=True)

