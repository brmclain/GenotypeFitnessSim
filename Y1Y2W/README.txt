MultiProcessingCode.py:

    MultiProcessingCode is a simulation to predict the frequencies of sex determining alleles of flys.
        It tracks the frequencies of each genotype, their fitness values, sex's, allele, and the 
        generation in which these frequencioes appreared to stop changing. The program iterates through
        the generations using the genotype frequencies of the previos to determine the current freuquencies. 
        Then it finds the generations where change stops, calculates all otehr returns from that generations,
        and places the results in folder "OutputArchive" with the name "<ID>_data.txt" as a row. Each row 
        is the end of a simulation from the same program execution. The parameters for each successful execution 
        of the program are written as a row in the file "ParameterLog.txt" in the folder "OutputArchive."

        To run this program by itself, several arguments are required to be passed through the command prompt.
        The required aurguments:
            1) simulations (whole number)
            2) Generations (whole number)
            3) CPU cores available (whole number)
            4) The dominance type (must be "dominant", "recessive", "random", or "additive")
        Optional aurguments:
             ) --control_file or -o (text file with a float value for every combination of alleles in the genotype that add up to 1. If none is given, will divide frequency evenly.)
                DEFAULT for Y1Y2W: [0.05556, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556]
             ) --allele_range or -a (value between 0 and 1, exluding 0)
                DEFAULT: 1
             ) --switches or -s (ones or zeros, 7 characters long to determine what data is placed in the final table.)
                [0] allele_fitness
                [1] fitness_array
                [2] dom_result
                [3] final_sim_result
                [4] allele_freq
                [5] sex_ratios
                [6] genotype_equil
                DEFAULT: "1111111"

    For the output, this program will create a folder if it does not exist called "OutputArchive" with a text file labels "data_<jobID>.txt" That contains a 2D array 
        with all the requested outputs. Along with this will be a file called "ParameterLog.txt" that is created if it does not exist that will place the jobID, the time
        it took to complete, and any input aurguments into a row within the file.









