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
             ) --control_file or -o (text file with 18 initial frequencies seperated by a tab and add up to 1)
                DEFAULT: [0.05556, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556]
             ) --allele_range or -a (value between 0 and 1, exluding 0)
                DEFAULT: 1
             ) --switches or -s (ones or zeros, 7 characters long to determine what is returned)
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

    ORDEROFPROGRESS:
        Detailed steps and the purpose of each code block in the order the program will execute. 
        Code blocks are labeled or inside function names

        Section A:
            This section handles the augments passed in the command line execution of this program and places them into variables for later use.
            Order of augments are:
        
        Section B:
            The starting "if" statement in this section is responsible for either retrieving the starting populations from the control file or if one is not present,
                creates default values to use instead.
            The data provided by the simulation, generation, and dominance are checked for uasability. Either raises assert error or for dominance will set value ro random.
            Takes user path for the output file and added root path to current directory to place output file back in desired location from node.

        Section C:
            Creation of dictionaries in order to make dataframes for needed data sets.

        Main:
            Calls the function "simu" to be done in parallel and returns the results to the array "data"
            Skip to function "simu" before continueing. 
            Comments used to select 1 of 4 options. Desired dataframes are pulled out from array "data" and concatinated into a new dataframe "result_data".
            Option 1:
                allele_fitness, fitness_array, dom_results, simulation_result, allele_freq, sex_ratio
            Option 2:
            Option 3:
            Option 4:
                dataout, dataout2, and genotype_equil
            Saves "result_data" to a .txt file as a csv.
            PROGRAM END

        Simu:
            Initializes results array
            Creates an array of random numbers with uniform distribution between the specified min/max for allele_fitness
            fitness_array created using the allele fitness and dominace in the function "genotype_fitness"
            Skip to function "Genotype_fitness" before continueing.
            !!COMMENTED OUT CODE!!
            Append first/initial population frequencies to simulation results.
            For each generation, run function "recursion_seln" with simulation_result with index current generation-1 and the fitness_array.
            Skip to function "recursion_seln" before continueing. 
            Returned array appended to varaible "simulation_result"
            convert variable "simulation_result" to numpy array "geno_sim_results".
            Run function "binary_search" individually for each collumn of geno_sim_results.
            Skip to function call "binary_search" before continueing.
            Place point of equalibrium(beginning of asymptot) in array "geno_asymp".
            !!COMMENTED OUT CODE!!
            Place max value from geno_asymp in varaible "max_geno_asymp".
            !!COMMENTED OUT CODE!!
            Place last value in array "simulation_result" in varaible "final_sim_result"
            Run fucntion "allele_freqs" with f"inal_sim_result" and place in variable "allele_freq"
            Skip to function "allele_freqs"
            Run fucntion "sex_ratio" with "final_sim_result" and place result in varaible "sex_ratios"
            Skip to function "sex_ratio".
            !!COMMENTED OUT CODE!!
            Return "allele_fitness", "fitness_array", "dom_results"(housekeeping done earlier), "final_sim_result",
                    "allele_freq", "sex_ratios", and "geno_asymp"
            Return to function "main"

        Genotype_fitness:
            Magic?
            return to function "simu"

        recursion_seln:
            Pass simulation_result to raw_count_ref and fitness_array to fitness_ref.
            Multiply each index from raw_count_ref and fitness_ref together into new list and normalize results.
            Place each index into seperate variable for each 18 variation.
            Calculate new genotype frequencies. place in new array "after_mating".
            Return normalized "after_mating" array.
            Go back to function "simu".

        binary_search:
            Finds and returns the value where the population frequency begins to reach an asymptots.
            Return to function "simu".

        allele_freqs:
            Calculates and returns the frequencies of each of the 5 alleles are present from the 18 possible combinations.
            Return to fucntion "simu".

        sex_ratio:
            Calculates and returns the ratio of males and females. 
            Return to fucntion "simu"







