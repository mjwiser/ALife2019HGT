#!/mnt/home/caninoko/anaconda3/bin/python

# Analyze HGT Fragment Information across different backgrounds
# RCK 10/20/18

## read the population file to get a reference of genomes
## for each line in the fitdist file
##  1. calculate fragment information in donor background
##     a. collect the fragment
##     b. findn where the fragment goes in the donor organism
##     c. calculate the information in the whole donor genome using KO
##     d. count the information in just the fragment [i_f]
##  2. calculate the informationn of the recipient
##     a. calculate the info in the whole recipient genome using KO
##     b. count the info in the region replaced by the fragment (as 
##        detailed by columns 8, 9)
##     c. having applied the fragment in the recipient, calculate the 
##        info in the whole recipient using KO
##     d. count the info in the fragment (against its new background)

## in order to do the above, generate a very long analyze file. :/ 


from optparse import OptionParser
import gzip
import pandas as pd
import subprocess
import sys

options = None ## global

def tidy_split(df, column, sep='|', keep=False):
    """
    Split the values of a column and expand so the new DataFrame has one split
    value per row. Filters rows where the column is missing.

    Params
    ------
    df : pandas.DataFrame
        dataframe with the column to split and expand
    column : str
        the column to split and expand
    sep : str
        the string used to split the column's values
    keep : bool
        whether to retain the presplit value as it's own row

    Returns
    -------
    pandas.DataFrame
        Returns a dataframe with the same columns as `df`.
    """
    indexes = list()
    new_values = list()
    df = df.dropna(subset=[column])
    for i, presplit in enumerate(df[column].astype(str)):
        values = presplit.split(sep)
        if keep and len(values) > 1:
            indexes.append(i)
            new_values.append(presplit)
        for value in values:
            indexes.append(i)
            new_values.append(value)
    new_df = df.iloc[indexes, :].copy()
    new_df[column] = new_values
    return new_df

def runProcess(cmd):

    p = subprocess.Popen(cmd.split(' '), stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE)
    return iter(p.stdout.readline, b'')  

def collect_parameters():
    # Set up options
    usage = """
     %prog [options] population.spop[.gz] fitdist.dat[.gz] environment.cfg outputfile.dat
    """
    parser = OptionParser(usage)
    parser.add_option("--debug_messages", action = "store_true", dest = "debug_messages",
                      default = False, help = "print debug messages to stdout")
    parser.add_option("--verbose", action = "store_true", dest = "verbose",
                      default = False, help = "print verbose messages")
    parser.add_option("--dryrun", action = "store_true", dest = "dryrun",
                      default = False, help = "don't run the big commands (avida)") 
    parser.add_option("--seed", dest = "seed",
                      default = 1000, help = "random seed for Avida analyze mode")
    ## fetch the args
    (opt, args) = parser.parse_args()
    if len(args) < 3:
        parser.error("incorrect number of arguments")

    popfile = args[0]
    fitdistfile = args[1]
    envfile = args[2]
    outputfilename = args[3]

    global options # use the global options variable
    options = opt

    return popfile,fitdistfile,envfile,outputfilename

def read_population_data(popfile):

    ## read the population file to get a reference of genomes
    fd = None
    if popfile[-3:] == ".gz":
        fd = gzip.open(popfile)
    else:
        fd = open(popfile)

    #for line in fd:
    #   print(line)

    population_data = pd.read_csv(fd, sep=" ", header=None, comment="#", skip_blank_lines=True,
        usecols=[4,11,16,17], names=['num_living', 'recipient_born', 'genome','occupied_cell_ids']
        )

    population_data['cell_id'] = population_data['occupied_cell_ids']

    ## explode the population into individual rows for each cell 


    population_data = tidy_split(population_data, 'cell_id', sep=',')
    population_data['cell_id'] = population_data['cell_id'].apply(lambda x: int(x))
    population_data['cell_id_num'] = population_data['cell_id']
    population_data = population_data.set_index('cell_id')

    return population_data

def donor_loc(row):

    fraglen = len(row['fragment'])
    #donorlen = len(row['donor_genome'])

    searchspace = row['donor_genome'] + row['donor_genome'][0:fraglen]

    loc = searchspace.find(row['fragment'])

    return loc

def recipient_prep(row):
    fraglen = len(row['fragment'])
    origgenomelen = len(row['genome'])
    insert_loc = row['insert_loc']

    mod_genome = row['genome'][insert_loc:] + row['genome'][:insert_loc] # rotate
    mod_genome = row['fragment'] + mod_genome[row['num_inst_replaced']:] # insert and remove
    mod_genome = mod_genome[-insert_loc:] + mod_genome[:-insert_loc] # rotate back

    return mod_genome

def prepare_recipient_analyze_commands(pop):

    analyzefilename = "TMP_analyze_recipients.cfg"

    command_section = """
LOAD_SEQUENCE $1 $2
"""
    recipient_command = ""
    for index, row in pop.iterrows():
        bit = command_section
        bit = bit.replace("$1", row['genome'])
        bit = bit.replace("$2", "recipient_cellID_" + str(index))
        recipient_command += bit

    recipient_command += """
RECALCULATE
MAP_TASKS phenotype_recipients/ fitness merit gest_time length viable task.0 task.1 task.2 task.3 task.4 task.5 task.6 task.7 task.8 total_task_count sequence
"""

    text_file = open(analyzefilename, "w")
    text_file.write(recipient_command)
    text_file.close()

    return analyzefilename

def prepare_donor_mutant_analyze_commands(df):

    analyzefilename = "TMP_analyze_donors_mutants.cfg"

    command_section = """
###### Analyze $2 and $4
LOAD_SEQUENCE $1 $2
LOAD_SEQUENCE $3 $4
"""
    recipient_command = ""
    for index, row in df.iterrows():
#       print(row)
        bit = command_section
        bit = bit.replace("$1", row['donor_genome'])
        bit = bit.replace("$2", "donor_fragmentnum_" + str(index) + "_cellID_" + str(row['cell_id']))
        bit = bit.replace("$3", row['genome_with_fragment'])
        bit = bit.replace("$4", "mutated_fragmentnum_" + str(index) + "_cellID_" + str(row['cell_id']))
                                                
        recipient_command += bit

    recipient_command += """
RECALCULATE
MAP_TASKS phenotype/ fitness merit gest_time length viable task.0 task.1 task.2 task.3 task.4 task.5 task.6 task.7 task.8 total_task_count sequence
"""

    text_file = open(analyzefilename, "w")
    text_file.write(recipient_command)
    text_file.close()

    return analyzefilename

def run(command):

    if options.verbose:
       print(command)
       print("---------------------------------------------------------")

    if options.dryrun:
        print( "DRY RUN")
    else:
        ct = 0
        sys.stdout.write(">-")
        for line in runProcess(command):
            ct += 1
            if options.verbose:
                print( "    ", line)
            elif ct % 50 == 0:
                sys.stdout.write("\b->")
                sys.stdout.flush()
        sys.stdout.write("\b-")
        print(">| done")

    if options.verbose:
       print( "#########################################################")      

def run_avida(analyzefile, envfile, datadir=None):

    ## run the Avida command. This may take a while.
    #command = "./avida -v4 -a -set DATA_DIR " + datadir + " -set ENVIRONMENT_FILE "+options.envfile
    command = ("time ./avida -a -v4 -s " + options.seed + " -set ANALYZE_FILE " + analyzefile +
               " -set HGT_UPTAKE_RECOMBINATION_P 0.1 -set ENVIRONMENT_FILE " + envfile +
               " -set EVENT_FILE events_savefrag__static.cfg")
    if datadir:
        command += " -set DATA_DIR " + datadir

    print( "#########################################################")
    print( "AVIDA RUN")
    print(command)

    run(command)

def get_recipient_stat(cellid):

    filename = "data_tmp/phenotype_recipients/tasksites.recipient_cellID_" + str(cellid) + ".dat"
    df = pd.read_csv(filename, sep=" ", header=None, 
        comment="#", skip_blank_lines=True,
        usecols=[0,3],
        names=['pos', 'fitness'])

    basefit = df.iloc[0]['fitness']

    tmp = df['fitness'].apply(lambda x: x < basefit) # where this site contributed to an increase in fitness.
    skel = list(map(int, tmp.tolist()[1:]))


    tmp = df['fitness'].apply(lambda x: x != basefit) # where this site contributed to fitness, period.
    skel_sensitive = list(map(int, tmp.tolist()[1:]))

    return (skel, skel_sensitive, basefit)

def collect_recipient_stats(pop):
    ## for each recipient (cell in the final population)
    ## 1. build an information skeleton for:
    ##    a. the original unmodified recipient (A/B)

    skel = pop["cell_id_num"].apply(lambda cellid: get_recipient_stat(cellid))
    pop['skeleton'] = skel.apply(lambda row: row[0])
    pop['skeleton_sensitive'] = skel.apply(lambda row: row[1])    
    pop['fitness'] = skel.apply(lambda row: row[2])

    return pop

def get_donor_mutant_stat(row):

    cell = int(row['cell_id_num'])
    fragnum = int(row.name)

    ####### DONOR SKELETON
    filename = "data_tmp/phenotype/tasksites.donor_fragmentnum_"+str(fragnum)+"_cellID_" + str(cell) + ".dat"
    df = pd.read_csv(filename, sep=" ", header=None, 
        comment="#", skip_blank_lines=True,
        usecols=[0,3],
        names=['pos', 'fitness'])

    basefit_donor = df.iloc[0]['fitness']
    tmp = df['fitness'].apply(lambda x: x < basefit_donor) # where this site contributed to an increase in fitness.
    skel_donor = list(map(int, tmp.tolist()[1:]))

    tmp = df['fitness'].apply(lambda x: x != basefit_donor) # where this site contributed to fitness, period.
    skel_donor_sensitive = list(map(int, tmp.tolist()[1:]))


    ####### MUTANT SKELETON
    filename = "data_tmp/phenotype/tasksites.mutated_fragmentnum_"+str(fragnum)+"_cellID_" + str(cell) + ".dat"
    df = pd.read_csv(filename, sep=" ", header=None, 
        comment="#", skip_blank_lines=True,
        usecols=[0,3],
        names=['pos', 'fitness'])

    basefit_mutant = df.iloc[0]['fitness']
    tmp = df['fitness'].apply(lambda x: x < basefit_mutant) # where this site contributed to an increase in fitness.
    skel_mutant = list(map(int, tmp.tolist()[1:]))

    tmp = df['fitness'].apply(lambda x: x != basefit_mutant) # where this site contributed to fitness, period.
    skel_mutant_sensitive = list(map(int, tmp.tolist()[1:]))


    return (skel_donor, skel_donor_sensitive,basefit_donor,skel_mutant, skel_mutant_sensitive, basefit_mutant)

def calc_donor_fragment_info(row):
    start = row['donor_fragment_loc']
    end = start + len(row['fragment'])
    fragment_info = sum(row['skeleton_donor'][start:end])

    return fragment_info

def calc_recipient_fragment_info(row):
    start = row['insert_loc']
    end = start + len(row['fragment'])
    fragment_info = sum(row['skeleton_mutant'][start:end])

    return fragment_info

def calc_recipient_removed_info(row):
    start = row['insert_loc']
    end = start + row['num_inst_replaced']
    removed_info = sum(row['skeleton'][start:end])

    return removed_info

def collect_donor_mutant_stats(chunk):
    ## for each line in this motherfucker
    ## 1. build an information skeleton for:
    ##    a. the donor
    ##    b. the fragment (based on the donor) 
    ##    c. the fragment (based on the recipient)
    ##    d. the count of removed information that was replaced by the fragment in the recipient

    skel = chunk.apply(lambda row: get_donor_mutant_stat(row), axis=1)
    chunk['skeleton_donor'] = skel.apply(lambda row: row[0])
    chunk['skeleton_donor_sensitive'] = skel.apply(lambda row: row[1])
    chunk['fitness_donor'] = skel.apply(lambda row: row[2])

    chunk['skeleton_mutant'] = skel.apply(lambda row: row[3])
    chunk['skeleton_mutant_sensitive'] = skel.apply(lambda row: row[4])
    chunk['fitness_mutant'] = skel.apply(lambda row: row[5])

    chunk['fragment_info_in_donor'] = chunk.apply(lambda row: calc_donor_fragment_info(row), axis=1)
    chunk['fragment_info_in_recipient'] = chunk.apply(lambda row: calc_recipient_fragment_info(row), axis=1)
    chunk['removed_info_in_recipient'] = chunk.apply(lambda row: calc_recipient_removed_info(row), axis=1)

    return chunk

def compose_output_header():
    line = "#"
    line += ' cell_id' # the cell of the recipient
    line += ' recipient_born' # the update the recipient was born
    line += ' fragment_id' # the unique fragment identifier
    line += ' donor_born' # the update that the donor was born
    line += ' len_skeleton' # the length of the recipient
    line += ' len_skeleton_donor' # the length of the donor
    line += ' len_skeleton_mutant' # the length of the mutant
    line += ' len_fragment' # the length of the fragment
    line += ' recipient_fitness' # recipient original fitness
    line += ' removed_info_in_recipient' # info removed from recipient to apply fragment
    line += ' donor_fitness' # the fitness of the fragment donor
    line += ' fragment_info_in_donor' # info in the fragment in donor background
    line += ' mutant_fitness' # the fitness of the recipient + fragment
    line += ' fragment_info_in_recipient' # info in the fragment in recipient background
    line += ' skeleton_recipient' # the skeleton of the recipient
    line += ' skeleton_donor' # the skeleton of the donor
    line += ' skeleton_mutant' # the skeleton of the mutant
    line += ' skeleton_recipient_sensitive' # the sensitive (any fitness change) skeleton of the recipient
    line += ' skeleton_donor_sensitive' # the sensitive (any fitness change) skeleton of the donor
    line += ' skeleton_mutant_sensitive' # the sensitive (any fitness change) skeleton of the mutant
    return line

def compose_output(row):
    line = []

    line.append( row['cell_id_num'] ) # the cell of the recipient organism
    line.append( row['recipient_born'] ) # the cell of the recipient organism
    line.append( row.name ) # the unique fragment identifier
    line.append( row['donor_born'] ) # the cell of the recipient organism
    line.append( len(row['skeleton']) )  # the length of the recipient
    line.append( len(row['skeleton_donor']) ) # the length of the donor
    line.append( len(row['skeleton_mutant']) ) # the length of the mutant
    line.append( len(row['fragment']) ) # the length of the fragment
    line.append( row['fitness'] ) # recipient original fitness
    line.append( row['removed_info_in_recipient'] ) # info removed from recipient to apply fragment
    line.append( row['fitness_donor'] ) # recipient original fitness
    line.append( row['fragment_info_in_donor'] ) # info in the fragment in donor background
    line.append( row['fitness_mutant'] ) # recipient original fitness
    line.append( row['fragment_info_in_recipient'] ) # info in the fragment in recipient background
    line.append( "".join([str(x) for x in row['skeleton']] ) ) # the skeleton of the recipient
    line.append( "".join([str(x) for x in row['skeleton_donor']] ) ) # the skeleton of the donor
    line.append( "".join([str(x) for x in row['skeleton_mutant']] ) ) # the skeleton of the mutant
    line.append( "".join([str(x) for x in row['skeleton_sensitive']] ) ) # the sensitive (any fitness change) skeleton of the recipient
    line.append( "".join([str(x) for x in row['skeleton_donor_sensitive']] ) ) # the sensitive (any fitness change) skeleton of the donor
    line.append( "".join([str(x) for x in row['skeleton_mutant_sensitive']] ) ) # the sensitive (any fitness change) skeleton of the mutant

    return " ".join([str(x) for x in line])

def output_stats(chunk, filename, header=False):
    text_file = open(filename, "a+")

    output = chunk.apply(lambda row: compose_output(row), axis=1)

    if header:
        text_file.write(compose_output_header())
        text_file.write("\n")

    output = "\n".join(output.tolist())

    text_file.write(output)
    text_file.write("\n")

    text_file.close()

def main():
    ## setup
    popfile, fitdistfile, envfile, outputfilename = collect_parameters()

    print("Collecting Recipient Data")
    print("=========================")

    print(" * Reading avida population")
    ## read the population data .spop file into an indexable dataframe
    population_data = read_population_data(popfile)
    print(" * Generating recipient population analyze mode commands")
    afile = prepare_recipient_analyze_commands(population_data)

    ## calculate the recipient data
    print(" * Generating recipient skeletons in Avida analyze mode")
    run_avida(afile, envfile, datadir="data_tmp/")

    ## collect the stats
    print(" * Collecting analyze mode output")
    population_data = collect_recipient_stats(population_data)

    # delete the temporary directories for this chunk
    run("rm -rf data_tmp/")

    print()
    print("Collecting Donor and Mutant Data")
    print("================================")

    ## prep for output
    run("rm " + outputfilename)

    fd = None
    if fitdistfile[-3:] == ".gz":
        fd = gzip.open(fitdistfile)
    else:
        fd = open(fitdistfile)

    ## split the work up into chunks
    chunkct = 0
    chunksize = 1000
    for chunk in pd.read_csv(fd, chunksize=chunksize, sep=" ", header=None, comment="#", skip_blank_lines=True,
        usecols=[0,1,7,8,9,10,11],
        names=['cell_id', 'donor_born', 'insert_loc', 'num_inst_replaced', 
               'num_inst_inserted', 'fragment', 'donor_genome']
        ):

        print()
        print("Chunk", chunkct)
        print('~~~~~~~~~~~~~~~~~~~~~~~~~')

        # connect the relevant population data ('num_living', 'genome','occupied_cell_ids')
        chunk = chunk.join(population_data, on='cell_id')

        # collect the various bits of information
        chunk['donor_fragment_loc'] = chunk.apply(lambda row: donor_loc(row),axis=1)
        chunk['genome_with_fragment'] = chunk.apply(lambda row: recipient_prep(row),axis=1)
        
        # generate the analyze commands   
        print(" * Generating donor and mutant analyze mode commands [",chunkct,"]")  
        afile = prepare_donor_mutant_analyze_commands(chunk)

        # run avida for each chunk against the two environments
        print(" * Generating donor and mutant skeletons in Avida analyze mode [",chunkct,"]")
        run_avida(afile, envfile, datadir="data_tmp/")

        # collect the stats
        print(" * Collecting analyze mode output [",chunkct,"]")
        chunk = collect_donor_mutant_stats(chunk)

        # output the stats somewhere
        print(" * Exporting the stats to ", outputfilename, " [",chunkct,"]")
        output_stats(chunk, outputfilename, (chunkct == 0))

        # delete the temporary directories for this chunk
        run("rm -rf data_tmp/")

        chunkct += 1

if __name__ == "__main__":
    main()


