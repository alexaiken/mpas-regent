from graphviz import Digraph
import re
import sys


def add_layer(desired_var, all_vars, dot, task_arr, seen_vars):
    
    seen_vars.append(desired_var)

    print("Desired var is", desired_var)

    if desired_var not in all_vars:
        all_vars.append(desired_var)
        dot.node(desired_var, desired_var)

    for task in task_arr:
        
        task_name = re.search('^[^\(]+', task).group(0).strip()
        reads_vars_output = []
        writes_vars_output = []

        writes_vars_chunk = re.search('(?s)writes\s+\(([^)]+)', task)
        if writes_vars_chunk is not None:
            writes_var_chunk_condensed = re.search('(?s)(?<=\().*', writes_vars_chunk.group(0)).group(0)

            writes_regions_arr = re.findall('(\w+)(?=\.\{)', writes_var_chunk_condensed)
            writes_vars_by_region= re.findall('\{([^}]+)', writes_var_chunk_condensed)

            for i in range(len(writes_vars_by_region)):
                curr_region = writes_regions_arr[i]
                curr_vars = writes_vars_by_region[i].split(',')

                for var in curr_vars:
                    new_var = curr_region.strip() + "." + var.strip()
                    writes_vars_output.append(new_var.strip())

        readswrites_var = re.search('(?s)reads writes\s+\(([^)]+)', task)
        if readswrites_var is not None:
            readswrites_var = re.search('(?s)(?<=\().*', readswrites_var.group(0)).group(0)

            readswrites_regions = re.findall('(\w+)(?=\.\{)', readswrites_var)
            readswrites_vars = re.findall('\{([^}]+)', readswrites_var)

            for i in range(len(readswrites_vars)):
                curr_region = readswrites_regions[i]
                curr_vars = readswrites_vars[i].split(',')

                for var in curr_vars:
                    new_var = curr_region.strip() + "." + var.strip()
                    writes_vars_output.append(new_var.strip())
                    reads_vars_output.append(new_var.strip())


        if (desired_var in writes_vars_output): 

            reads_var = re.search('(?s)reads\s+\(([^)]+)', task)
            if reads_var is not None:
                reads_var = re.search('(?s)(?<=\().*', reads_var.group(0)).group(0)
                reads_regions = re.findall('(\w+)(?=\.\{)', reads_var)
                reads_vars = re.findall('\{([^}]+)', reads_var)

                for i in range(len(reads_vars)):
                    curr_region = reads_regions[i]
                    curr_vars = reads_vars[i].split(',')

                    for var in curr_vars:
                        new_var = (curr_region.strip() + "." + var.strip()).strip()
                        reads_vars_output.append(new_var)
                        
                        
            
            
            for read_var in reads_vars_output:
                if new_var not in all_vars:
                    dot.node(new_var, new_var)
                    all_vars.append(new_var)
                dot.edge(read_var, desired_var, URL=task_name, Tooltip = task_name)




            for reads_var in reads_vars_output:
                if (reads_var not in seen_vars):
                    add_layer(reads_var, all_vars, dot, task_arr, seen_vars)



dot = Digraph(comment='MPAS', engine ='fdp')

textfile = open('input.rg', 'r')
filetext = textfile.read()
textfile.close()

task_arr = re.findall('(?s)(?<=task)(.*?)(?=do\n)', filetext)
all_vars = []
seen_vars = []

all_vars_bool = True if sys.argv[1] == "all" else False

if (not all_vars_bool): 
    desired_var = sys.argv[1] 
    print("Desired var is ", desired_var)

add_layer(desired_var, all_vars, dot, task_arr, seen_vars)

u = dot.unflatten(stagger=5)
u.render('output/dynamics_tasks.gv', view=True)

