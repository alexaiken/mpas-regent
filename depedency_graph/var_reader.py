
import re


def populate_read_vars(task, reads_vars_output):
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

def populate_write_vars(task, writes_vars_output):

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

def populate_readwrite_vars(task, reads_vars_output, writes_vars_output):
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
