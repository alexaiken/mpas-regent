
import re


def read_var(task, reads_vars_output):
    reads_var = re.search('(?s)reads\s+\(([^)]+)', task)
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