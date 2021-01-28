from graphviz import Digraph
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
    desired_vars = [str(arg) for arg in sys.argv[1:]] #['cr.pressure_p', 'cr.rho']

for task in task_arr:
    
    task_name = re.search('^[^\(]+', task).group(0).strip()
    print(task_name)
    reads_vars_output = []
    writes_vars_output = []
 
    reads_var = re.search('(?s)reads\s+\(([^)]+)', task)
    print(reads_var)
    if reads_var is not None:
        reads_var = re.search('(?s)(?<=\().*', reads_var.group(0)).group(0)
        reads_regions = re.findall('(\w+)(?=\.\{)', reads_var)
        reads_vars = re.findall('\{([^}]+)', reads_var)

        

        for i in range(len(reads_vars)):
            curr_region = reads_regions[i]
            curr_vars = reads_vars[i].split(',')

            for var in curr_vars:
                new_var = curr_region.strip() + "." + var.strip()
                reads_vars_output.append(new_var.strip())
                #if new_var not in all_vars:
                #    dot.node(new_var, new_var)
                print(new_var)


    writes_var = re.search('(?s)writes\s+\(([^)]+)', task)
    if writes_var is not None:
        writes_var = re.search('(?s)(?<=\().*', writes_var.group(0)).group(0)

        writes_regions = re.findall('(\w+)(?=\.\{)', writes_var)
        writes_vars = re.findall('\{([^}]+)', writes_var)

        

        for i in range(len(writes_vars)):
            curr_region = writes_regions[i]
            curr_vars = writes_vars[i].split(',')

            for var in curr_vars:
                new_var = curr_region.strip() + "." + var.strip()
                writes_vars_output.append(new_var.strip())
                #if new_var not in all_vars:
                #    dot.node(new_var, new_var)
                print(new_var)

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
                #if new_var not in all_vars:
                #    dot.node(new_var, new_var)
                print(new_var)

    print("reads vars:", reads_vars_output)
    print("writes vars:", writes_vars_output)

    for reads_var in reads_vars_output:
        for write_var in writes_vars_output:
            if all_vars_bool: 
                if write_var not in all_vars:
                    dot.node(write_var, write_var)
                if reads_var not in all_vars:
                    dot.node(reads_var, reads_var)
                dot.edge(reads_var, write_var, URL=task_name, Tooltip = task_name)

            elif write_var in desired_vars:
                if write_var not in all_vars:
                    dot.node(write_var, write_var)
                if reads_var not in all_vars:
                    dot.node(reads_var, reads_var)
                dot.edge(reads_var, write_var, URL=task_name, Tooltip = task_name)

    


#reads_arr = ['cr.edgesOnCell', 'cr.nEdgesOnCell', 'cr.verticesOnCell', 'er.cellsOnEdge', 'er.cellsOnEdge', 'er.verticesOnEdge', 'er.zb', 'er.zb3', 'vr.cellsOnVertex', 'vr.edgesOnVertex']
#writes_arr = ['cr.edgesOnCellSign', 'cr.kiteForCell', 'cr.zb_cell', 'cr.zb3_cell',  'vr.edgesOnVertexSign']

#task_name = 'atm_compute_signs'

#for var in reads_arr:
    #dot.node(var, var)

#for var in writes_arr:
    #dot.node(var, var)

#for reads_var in reads_arr:
#    for write_var in writes_arr:
#        dot.edge(reads_var, write_var)
        #dot.edge(reads_var, write_var, constraint='false', label=task_name)

#print(dot.source)  



u = dot.unflatten(stagger=5)


u.render('dynamics_tasks.gv', view=True)

