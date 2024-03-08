#!/usr/bin/perl
# Run recursions to explore sample space of selection coefficients on individual alleles
#  with arbitrary dominance of fitness effects 
#  and arbitrary starting allele/chromosome frequencies
# System has: 
#	1. X1Y1
#	2. ZW
#	3. X2Y2
# Y1 < W < Y2

# This approach allows for dominant, recessive, or additive fitness effects


# cd ~/Projects/Calyptrates/Musca/PopGenModel2021
# perl seln_allele_dominance_Y1WY2.pl generations simulations output
# perl seln_allele_dominance_Y1WY2.pl 100 10 seln_allele_dominance_Y1WY2.out


# In each parameterization: 
#	1. Draw random frequencies of Y2, Y1, and W
#	2. Perform 10 generations of random mating to determine starting frequencies of each genotype
#	3. Draw random dominance and fitness values for each chromosome
#	4. Calculate multi-chromosome genotype fitness value for all 18 genotypes
#	5. Simulate for 1,000 generations to determine freqs after selection

use strict;
use warnings;
use List::Util qw(sum min max);
use Math::Random qw(random_uniform); # for random number generation -- random_uniform($n, $low, $high)
use POSIX; # for rounding

my $simulations=10; # number of simulations to perform
my $generations=100; # number of generations to run each simulation (can be reset by user)
my $outfile;


# 1	female	f01	X1X1	ZZ	X2X2	0.5*f01*m01 + 0.25*f02*m01 + 0.125*f04*m01 + 0.25*f01*m02 + 0.125*f02*m02 + 0.0625*f04*m02 + 0.5*f01*m04 + 0.25*f02*m04 + 0.125*f04*m04 + 0.25*f01*m05 + 0.125*f02*m05 + 0.0625*f04*m05 + 0.125*f01*m06 + 0.0625*f02*m06 + 0.03125*f04*m06
# 3	female	f02	X1X1	ZW	X2X2	0.25*f02*m01 + 0.5*f03*m01 + 0.125*f04*m01 + 0.25*f05*m01 + 0.25*f01*m02 + 0.25*f02*m02 + 0.25*f03*m02 + 0.125*f04*m02 + 0.125*f05*m02 + 0.5*f01*m03 + 0.25*f02*m03 + 0.125*f04*m03 + 0.25*f02*m04 + 0.5*f03*m04 + 0.125*f04*m04 + 0.25*f05*m04 + 0.125*f02*m05 + 0.25*f03*m05 + 0.0625*f04*m05 + 0.125*f05*m05 + 0.125*f01*m06 + 0.125*f02*m06 + 0.125*f03*m06 + 0.0625*f04*m06 + 0.0625*f05*m06 + 0.25*f01*m07 + 0.125*f02*m07 + 0.0625*f04*m07
# 5	female	f03	X1X1	WW	X2X2	0.125*f02*m02 + 0.25*f03*m02 + 0.0625*f04*m02 + 0.125*f05*m02 + 0.25*f02*m03 + 0.5*f03*m03 + 0.125*f04*m03 + 0.25*f05*m03 + 0.0625*f02*m06 + 0.125*f03*m06 + 0.03125*f04*m06 + 0.0625*f05*m06 + 0.125*f02*m07 + 0.25*f03*m07 + 0.0625*f04*m07 + 0.125*f05*m07
# 9	female	f04	X1Y1	ZW	X2X2	0.125*f04*m01 + 0.25*f05*m01 + 0.25*f06*m01 + 0.5*f07*m01 + 0.125*f04*m02 + 0.125*f05*m02 + 0.25*f06*m02 + 0.25*f07*m02 + 0.125*f04*m03 + 0.25*f06*m03 + 0.25*f02*m04 + 0.5*f03*m04 + 0.25*f04*m04 + 0.5*f05*m04 + 0.25*f06*m04 + 0.5*f07*m04 + 0.125*f02*m05 + 0.25*f03*m05 + 0.125*f04*m05 + 0.25*f05*m05 + 0.125*f06*m05 + 0.25*f07*m05 + 0.125*f01*m06 + 0.125*f02*m06 + 0.125*f03*m06 + 0.125*f04*m06 + 0.125*f05*m06 + 0.125*f06*m06 + 0.125*f07*m06 + 0.25*f01*m07 + 0.125*f02*m07 + 0.125*f04*m07 + 0.125*f06*m07 + 0.5*f02*m08 + 1*f03*m08 + 0.25*f04*m08 + 0.5*f05*m08 + 0.25*f02*m09 + 0.5*f03*m09 + 0.125*f04*m09 + 0.25*f05*m09 + 0.25*f01*m10 + 0.25*f02*m10 + 0.25*f03*m10 + 0.125*f04*m10 + 0.125*f05*m10 + 0.5*f01*m11 + 0.25*f02*m11 + 0.125*f04*m11
#11	female	f05	X1Y1	WW	X2X2	0.0625*f04*m02 + 0.125*f05*m02 + 0.125*f06*m02 + 0.25*f07*m02 + 0.125*f04*m03 + 0.25*f05*m03 + 0.25*f06*m03 + 0.5*f07*m03 + 0.0625*f02*m06 + 0.125*f03*m06 + 0.0625*f04*m06 + 0.125*f05*m06 + 0.0625*f06*m06 + 0.125*f07*m06 + 0.125*f02*m07 + 0.25*f03*m07 + 0.125*f04*m07 + 0.25*f05*m07 + 0.125*f06*m07 + 0.25*f07*m07 + 0.125*f02*m10 + 0.25*f03*m10 + 0.0625*f04*m10 + 0.125*f05*m10 + 0.25*f02*m11 + 0.5*f03*m11 + 0.125*f04*m11 + 0.25*f05*m11
#15	female	f06	Y1Y1	ZW	X2X2	0.125*f04*m04 + 0.25*f05*m04 + 0.25*f06*m04 + 0.5*f07*m04 + 0.0625*f04*m05 + 0.125*f05*m05 + 0.125*f06*m05 + 0.25*f07*m05 + 0.0625*f04*m06 + 0.0625*f05*m06 + 0.125*f06*m06 + 0.125*f07*m06 + 0.0625*f04*m07 + 0.125*f06*m07 + 0.25*f04*m08 + 0.5*f05*m08 + 0.5*f06*m08 + 1*f07*m08 + 0.125*f04*m09 + 0.25*f05*m09 + 0.25*f06*m09 + 0.5*f07*m09 + 0.125*f04*m10 + 0.125*f05*m10 + 0.25*f06*m10 + 0.25*f07*m10 + 0.125*f04*m11 + 0.25*f06*m11
#17	female	f07	Y1Y1	WW	X2X2	0.03125*f04*m06 + 0.0625*f05*m06 + 0.0625*f06*m06 + 0.125*f07*m06 + 0.0625*f04*m07 + 0.125*f05*m07 + 0.125*f06*m07 + 0.25*f07*m07 + 0.0625*f04*m10 + 0.125*f05*m10 + 0.125*f06*m10 + 0.25*f07*m10 + 0.125*f04*m11 + 0.25*f05*m11 + 0.25*f06*m11 + 0.5*f07*m11
# 2	male	m01	X1X1	ZZ	X2Y2	0.5*f01*m01 + 0.25*f02*m01 + 0.125*f04*m01 + 0.25*f01*m02 + 0.125*f02*m02 + 0.0625*f04*m02 + 0.25*f01*m05 + 0.125*f02*m05 + 0.0625*f04*m05 + 0.125*f01*m06 + 0.0625*f02*m06 + 0.03125*f04*m06
# 4	male	m02	X1X1	ZW	X2Y2	0.25*f02*m01 + 0.5*f03*m01 + 0.125*f04*m01 + 0.25*f05*m01 + 0.25*f01*m02 + 0.25*f02*m02 + 0.25*f03*m02 + 0.125*f04*m02 + 0.125*f05*m02 + 0.5*f01*m03 + 0.25*f02*m03 + 0.125*f04*m03 + 0.125*f02*m05 + 0.25*f03*m05 + 0.0625*f04*m05 + 0.125*f05*m05 + 0.125*f01*m06 + 0.125*f02*m06 + 0.125*f03*m06 + 0.0625*f04*m06 + 0.0625*f05*m06 + 0.25*f01*m07 + 0.125*f02*m07 + 0.0625*f04*m07
# 6	male	m03	X1X1	WW	X2Y2	0.125*f02*m02 + 0.25*f03*m02 + 0.0625*f04*m02 + 0.125*f05*m02 + 0.25*f02*m03 + 0.5*f03*m03 + 0.125*f04*m03 + 0.25*f05*m03 + 0.0625*f02*m06 + 0.125*f03*m06 + 0.03125*f04*m06 + 0.0625*f05*m06 + 0.125*f02*m07 + 0.25*f03*m07 + 0.0625*f04*m07 + 0.125*f05*m07
# 7	male	m04	X1Y1	ZZ	X2X2	0.125*f04*m01 + 0.25*f06*m01 + 0.0625*f04*m02 + 0.125*f06*m02 + 0.5*f01*m04 + 0.25*f02*m04 + 0.25*f04*m04 + 0.25*f06*m04 + 0.25*f01*m05 + 0.125*f02*m05 + 0.125*f04*m05 + 0.125*f06*m05 + 0.125*f01*m06 + 0.0625*f02*m06 + 0.0625*f04*m06 + 0.0625*f06*m06 + 1*f01*m08 + 0.5*f02*m08 + 0.25*f04*m08 + 0.5*f01*m09 + 0.25*f02*m09 + 0.125*f04*m09 + 0.25*f01*m10 + 0.125*f02*m10 + 0.0625*f04*m10
# 8	male	m05	X1Y1	ZZ	X2Y2	0.125*f04*m01 + 0.25*f06*m01 + 0.0625*f04*m02 + 0.125*f06*m02 + 0.25*f01*m05 + 0.125*f02*m05 + 0.125*f04*m05 + 0.125*f06*m05 + 0.125*f01*m06 + 0.0625*f02*m06 + 0.0625*f04*m06 + 0.0625*f06*m06 + 0.5*f01*m09 + 0.25*f02*m09 + 0.125*f04*m09 + 0.25*f01*m10 + 0.125*f02*m10 + 0.0625*f04*m10
#10	male	m06	X1Y1	ZW	X2Y2	0.125*f04*m01 + 0.25*f05*m01 + 0.25*f06*m01 + 0.5*f07*m01 + 0.125*f04*m02 + 0.125*f05*m02 + 0.25*f06*m02 + 0.25*f07*m02 + 0.125*f04*m03 + 0.25*f06*m03 + 0.125*f02*m05 + 0.25*f03*m05 + 0.125*f04*m05 + 0.25*f05*m05 + 0.125*f06*m05 + 0.25*f07*m05 + 0.125*f01*m06 + 0.125*f02*m06 + 0.125*f03*m06 + 0.125*f04*m06 + 0.125*f05*m06 + 0.125*f06*m06 + 0.125*f07*m06 + 0.25*f01*m07 + 0.125*f02*m07 + 0.125*f04*m07 + 0.125*f06*m07 + 0.25*f02*m09 + 0.5*f03*m09 + 0.125*f04*m09 + 0.25*f05*m09 + 0.25*f01*m10 + 0.25*f02*m10 + 0.25*f03*m10 + 0.125*f04*m10 + 0.125*f05*m10 + 0.5*f01*m11 + 0.25*f02*m11 + 0.125*f04*m11
#12	male	m07	X1Y1	WW	X2Y2	0.0625*f04*m02 + 0.125*f05*m02 + 0.125*f06*m02 + 0.25*f07*m02 + 0.125*f04*m03 + 0.25*f05*m03 + 0.25*f06*m03 + 0.5*f07*m03 + 0.0625*f02*m06 + 0.125*f03*m06 + 0.0625*f04*m06 + 0.125*f05*m06 + 0.0625*f06*m06 + 0.125*f07*m06 + 0.125*f02*m07 + 0.25*f03*m07 + 0.125*f04*m07 + 0.25*f05*m07 + 0.125*f06*m07 + 0.25*f07*m07 + 0.125*f02*m10 + 0.25*f03*m10 + 0.0625*f04*m10 + 0.125*f05*m10 + 0.25*f02*m11 + 0.5*f03*m11 + 0.125*f04*m11 + 0.25*f05*m11
#13	male	m08	Y1Y1	ZZ	X2X2	0.125*f04*m04 + 0.25*f06*m04 + 0.0625*f04*m05 + 0.125*f06*m05 + 0.03125*f04*m06 + 0.0625*f06*m06 + 0.25*f04*m08 + 0.5*f06*m08 + 0.125*f04*m09 + 0.25*f06*m09 + 0.0625*f04*m10 + 0.125*f06*m10
#14	male	m09	Y1Y1	ZZ	X2Y2	0.0625*f04*m05 + 0.125*f06*m05 + 0.03125*f04*m06 + 0.0625*f06*m06 + 0.125*f04*m09 + 0.25*f06*m09 + 0.0625*f04*m10 + 0.125*f06*m10
#16	male	m10	Y1Y1	ZW	X2Y2	0.0625*f04*m05 + 0.125*f05*m05 + 0.125*f06*m05 + 0.25*f07*m05 + 0.0625*f04*m06 + 0.0625*f05*m06 + 0.125*f06*m06 + 0.125*f07*m06 + 0.0625*f04*m07 + 0.125*f06*m07 + 0.125*f04*m09 + 0.25*f05*m09 + 0.25*f06*m09 + 0.5*f07*m09 + 0.125*f04*m10 + 0.125*f05*m10 + 0.25*f06*m10 + 0.25*f07*m10 + 0.125*f04*m11 + 0.25*f06*m11
#18	male	m11	Y1Y1	WW	X2Y2	0.03125*f04*m06 + 0.0625*f05*m06 + 0.0625*f06*m06 + 0.125*f07*m06 + 0.0625*f04*m07 + 0.125*f05*m07 + 0.125*f06*m07 + 0.25*f07*m07 + 0.0625*f04*m10 + 0.125*f05*m10 + 0.125*f06*m10 + 0.25*f07*m10 + 0.125*f04*m11 + 0.25*f05*m11 + 0.25*f06*m11 + 0.5*f07*m11


($generations, $simulations, $outfile) = @ARGV;

# Make sure number of simulations are > 0
if($simulations < 1){
	die "\nSimulations must be >0\n";
}

# Make sure generations are >0
if($generations<1){
	die "\nGenerations must be >0\n";
}


# Print header
open OUTPUT, ">$outfile" or die "\nCould not open output file, $outfile\n";
print OUTPUT "Y1_fem\tW_fem\tY1_male\tY2_male\tW_male\t"; # fitness values
print OUTPUT "female_Y_dom\tmale_Y_dom\tfemale_W_dom\tmale_W_dom\t"; # dominance effects of proto-sex chromosomes
print OUTPUT "X1_i\tY1_i\tX2_i\tY2_i\tZ_i\tW_i\t"; # initial chromosome freqs
print OUTPUT "f01_i\tf02_i\tf03_i\tf04_i\tf05_i\tf06_i\tf07_i\t"; # initial genotype freqs (female)
print OUTPUT "m01_i\tm02_i\tm03_i\tm04_i\tm05_i\tm06_i\tm07_i\tm08_i\tm09_i\tm10_i\tm11_i\t"; # initial genotype freqs (male)
print OUTPUT "f01_fit\tf02_fit\tf03_fit\tf04_fit\tf05_fit\tf06_fit\tf07_fit\t"; # genotype fitness values (female)
print OUTPUT "m01_fit\tm02_fit\tm03_fit\tm04_fit\tm05_fit\tm06_fit\tm07_fit\tm08_fit\tm09_fit\tm10_fit\tm11_fit\t"; # genotype fitness values (male)
print OUTPUT "f01_f\tf02_f\tf03_f\tf04_f\tf05_f\tf06_f\tf07_f\t"; # final genotype freqs (female)
print OUTPUT "m01_f\tm02_f\tm03_f\tm04_f\tm05_f\tm06_f\tm07_f\tm08_f\tm09_f\tm10_f\tm11_f\t"; # final genotype freqs (male)
print OUTPUT "X1_f\tY1_f\tX2_f\tY2_f\tZ_f\tW_f\t"; # final chromosome freqs
print OUTPUT "parameterization\n"; 


# Run simulations
for (my $i=0; $i<$simulations; $i++){

#	1. Draw random frequencies of Y2, Y1, and W
	my @allele_freq = random_uniform(3, 0, 1); 	

#	2. Perform 5 generations of random mating to determine starting frequencies of each genotype
	my @initial_freq = genotype_freq(@allele_freq); # initial genotype frequencies

#	3. Draw random dominance and fitness values for each chromosome
	# set allele fitness values for each allele in each sex
	#	0 Y1 	female
	#	1 W		female
	#	2 Y1	male
	#	3 Y2   	male
	#	4 W		male
	my @allele_fitness = random_uniform(5, -1, 1); 

#	4. Calculate multi-chromosome genotype fitness value for all 18 genotypes
	# calculate genotype fitness values 
	my @dominance = random_uniform(4, 0, 1); # female & male, Y and W dominance effects
	my @fitness_array = genotype_fitness(\@allele_fitness, \@dominance); # need to pass the function references to the arrays, not the arrays themselves

#	5. Simulate for 1,000 generations to determine freqs after selection

	# Run the population simulation
	my @simulation_results;
	$simulation_results[0] = [ @initial_freq ];
	for(my $a=1; $a <= $generations; $a++){
		$simulation_results[$a] = [ recursion_seln(\@{$simulation_results[$a-1]}, \@fitness_array) ]; # need to pass the function references to the arrays, not the arrays themselves
	}
	


	# Print final results of simulation to file

	# print OUTPUT "Y1_fem\tW_fem\tY1_male\tY2_male\tW_male\t"; # fitness values
	for(my $b=0; $b<5; $b++){
		print OUTPUT "$allele_fitness[$b]\t";
	}

	# print OUTPUT "female_Y_dom\tmale_Y_dom\tfemale_W_dom\tmale_W_dom\t"; # dominance effects of proto-sex chromosomes
	print OUTPUT "$dominance[0]\t"; # female Y dominance
	print OUTPUT "$dominance[1]\t"; # male Y dominance
	print OUTPUT "$dominance[2]\t"; # female W dominance
	print OUTPUT "$dominance[3]\t"; # male W dominance

	# print OUTPUT "X1_i\tY1_i\tX2_i\tY2_i\tZ_i\tW_i\t"; # initial chromosome freqs
	my @allele_freq_initial = allele_freqs(@initial_freq); 
	for(my $b=0; $b<6; $b++){
		print OUTPUT "$allele_freq_initial[$b]\t";
	}

	# initial genotype freqs
	# print OUTPUT "f01_i\tf02_i\tf03_i\tf04_i\tf05_i\tf06_i\tf07_i\t"; # initial genotype freqs (female)
	# print OUTPUT "m01_i\tm02_i\tm03_i\tm04_i\tm05_i\tm06_i\tm07_i\tm08_i\tm09_i\tm10_i\tm11_i\t"; # initial genotype freqs (male)
	for(my $b=0; $b<18; $b++){
		print OUTPUT "$initial_freq[$b]\t"
	}
	
	# fitness values
	# print OUTPUT "f01_fit\tf02_fit\tf03_fit\tf04_fit\tf05_fit\tf06_fit\tf07_fit\t"; # genotype fitness values (female)
	# print OUTPUT "m01_fit\tm02_fit\tm03_fit\tm04_fit\tm05_fit\tm06_fit\tm07_fit\tm08_fit\tm09_fit\tm10_fit\tm11_fit\t"; # genotype fitness values (male)
	for(my $b=0; $b<18; $b++){
		print OUTPUT "$fitness_array[$b]\t";
	}

	# final genotype freqs
	# print OUTPUT "f01_f\tf02_f\tf03_f\tf04_f\tf05_f\tf06_f\tf07_f\t"; # final genotype freqs (female)
	# print OUTPUT "m01_f\tm02_f\tm03_f\tm04_f\tm05_f\tm06_f\tm07_f\tm08_f\tm09_f\tm10_f\tm11_f\t"; # final genotype freqs (male)
	for(my $b=0; $b<18; $b++){
		print OUTPUT "$simulation_results[$generations][$b]\t";
	}

	# final chromosome freqs
	# print OUTPUT "X1_f\tY1_f\tX2_f\tY2_f\tZ_f\tW_f\t"; # final chromosome freqs
	my @allele_freq_final = allele_freqs(@{$simulation_results[$generations]});
	for(my $b=0; $b<6; $b++){
		print OUTPUT "$allele_freq_final[$b]\t";
	}
	
	my $parameterization = $i+1;
	print OUTPUT "$parameterization\n";
													
}
close OUTPUT;



# Function performs 5 generations of random mating to determine initial genotype frequencies
sub genotype_freq{

	# store allele frequencies into array
	# print OUTPUT "X1_i\tY1_i\tX2_i\tY2_i\tZ_i\tW_i\t"; # initial chromosome freqs
	my @alleles = @_;
	my $Y1 = $alleles[0];
	my $X1 = 1 - $alleles[0];
	my $Y2 = $alleles[1];
	my $X2 = 1 - $alleles[1];
	my $W = $alleles[2];
	my $Z = 1 - $alleles[2];

	# use HW to determine genotype frequencies
	my @genotypes = (		
		($X1*$X1)*($Z*$Z)*($X2*$X2),		#female	f01	X1X1	ZZ	X2X2
		($X1*$X1)*(2*$W*$Z)*($X2*$X2),		#female	f02	X1X1	ZW	X2X2
		($X1*$X1)*($W*$W)*($X2*$X2),		#female	f03	X1X1	WW	X2X2
		(2*$X1*$Y1)*(2*$W*$Z)*($X2*$X2),	#female	f04	X1Y1	ZW	X2X2
		(2*$X1*$Y1)*($W*$W)*($X2*$X2),		#female	f05	X1Y1	WW	X2X2
		($Y1*$Y1)*(2*$W*$Z)*($X2*$X2),		#female	f06	Y1Y1	ZW	X2X2
		($Y1*$Y1)*($W*$W)*($X2*$X2),		#female	f07	Y1Y1	WW	X2X2
	 	($X1*$X1)*($Z*$Z)*(2*$X2*$Y2),		#male	m01	X1X1	ZZ	X2Y2
	 	($X1*$X1)*(2*$W*$Z)*(2*$X2*$Y2),	#male	m02	X1X1	ZW	X2Y2
	 	($X1*$X1)*($W*$W)*(2*$X2*$Y2),		#male	m03	X1X1	WW	X2Y2
	 	(2*$X1*$Y1)*($Z*$Z)*($X2*$X2),		#male	m04	X1Y1	ZZ	X2X2
	 	(2*$X1*$Y1)*($Z*$Z)*(2*$X2*$Y2),	#male	m05	X1Y1	ZZ	X2Y2
	 	(2*$X1*$Y1)*(2*$Z*$W)*(2*$X2*$Y2),	#male	m06	X1Y1	ZW	X2Y2
	 	(2*$X1*$Y1)*($W*$W)*(2*$X2*$Y2),	#male	m07	X1Y1	WW	X2Y2
	 	($Y1*$Y1)*($Z*$Z)*($X2*$X2),		#male	m08	Y1Y1	ZZ	X2X2
	 	($Y1*$Y1)*($Z*$Z)*(2*$X2*$Y2),		#male	m09	Y1Y1	ZZ	X2Y2
	 	($Y1*$Y1)*(2*$Z*$W)*(2*$X2*$Y2),	#male	m10	Y1Y1	ZW	X2Y2
	 	($Y1*$Y1)*($W*$W)*(2*$X2*$Y2)		#male	m11	Y1Y1	WW	X2Y2
	);
	my @genotypes_std;
	for (my $i=0; $i<18; $i++) {
		$genotypes_std[$i] = $genotypes[$i] / sum(@genotypes);
	}
	my @genotypes_1 = random_mating(@genotypes_std);
	my @genotypes_2 = random_mating(@genotypes_1);
	my @genotypes_3 = random_mating(@genotypes_2);
	my @genotypes_4 = random_mating(@genotypes_3);
	my @genotypes_5 = random_mating(@genotypes_4);
	@genotypes_5;
}


# function calculates new genotype frequencies after one generation of random mating
sub random_mating{
	
	my @counts = @_; 
	
	my $f01 = $counts[0];
	my $f02 = $counts[1];
	my $f03 = $counts[2];
	my $f04 = $counts[3];
	my $f05 = $counts[4];
	my $f06 = $counts[5];
	my $f07 = $counts[6];

	my $m01 = $counts[7];
	my $m02 = $counts[8];
	my $m03 = $counts[9];
	my $m04 = $counts[10];
	my $m05 = $counts[11];
	my $m06 = $counts[12];
	my $m07 = $counts[13];
	my $m08 = $counts[14];
	my $m09 = $counts[15];
	my $m10 = $counts[16];
	my $m11 = $counts[17]; 
	
	# Calculate genotype frequencies after random mating
	my @after_mating = (
		0.5*$f01*$m01 + 0.25*$f02*$m01 + 0.125*$f04*$m01 + 0.25*$f01*$m02 + 0.125*$f02*$m02 + 0.0625*$f04*$m02 + 0.5*$f01*$m04 + 0.25*$f02*$m04 + 0.125*$f04*$m04 + 0.25*$f01*$m05 + 0.125*$f02*$m05 + 0.0625*$f04*$m05 + 0.125*$f01*$m06 + 0.0625*$f02*$m06 + 0.03125*$f04*$m06,
		0.25*$f02*$m01 + 0.5*$f03*$m01 + 0.125*$f04*$m01 + 0.25*$f05*$m01 + 0.25*$f01*$m02 + 0.25*$f02*$m02 + 0.25*$f03*$m02 + 0.125*$f04*$m02 + 0.125*$f05*$m02 + 0.5*$f01*$m03 + 0.25*$f02*$m03 + 0.125*$f04*$m03 + 0.25*$f02*$m04 + 0.5*$f03*$m04 + 0.125*$f04*$m04 + 0.25*$f05*$m04 + 0.125*$f02*$m05 + 0.25*$f03*$m05 + 0.0625*$f04*$m05 + 0.125*$f05*$m05 + 0.125*$f01*$m06 + 0.125*$f02*$m06 + 0.125*$f03*$m06 + 0.0625*$f04*$m06 + 0.0625*$f05*$m06 + 0.25*$f01*$m07 + 0.125*$f02*$m07 + 0.0625*$f04*$m07,
		0.125*$f02*$m02 + 0.25*$f03*$m02 + 0.0625*$f04*$m02 + 0.125*$f05*$m02 + 0.25*$f02*$m03 + 0.5*$f03*$m03 + 0.125*$f04*$m03 + 0.25*$f05*$m03 + 0.0625*$f02*$m06 + 0.125*$f03*$m06 + 0.03125*$f04*$m06 + 0.0625*$f05*$m06 + 0.125*$f02*$m07 + 0.25*$f03*$m07 + 0.0625*$f04*$m07 + 0.125*$f05*$m07,
		0.125*$f04*$m01 + 0.25*$f05*$m01 + 0.25*$f06*$m01 + 0.5*$f07*$m01 + 0.125*$f04*$m02 + 0.125*$f05*$m02 + 0.25*$f06*$m02 + 0.25*$f07*$m02 + 0.125*$f04*$m03 + 0.25*$f06*$m03 + 0.25*$f02*$m04 + 0.5*$f03*$m04 + 0.25*$f04*$m04 + 0.5*$f05*$m04 + 0.25*$f06*$m04 + 0.5*$f07*$m04 + 0.125*$f02*$m05 + 0.25*$f03*$m05 + 0.125*$f04*$m05 + 0.25*$f05*$m05 + 0.125*$f06*$m05 + 0.25*$f07*$m05 + 0.125*$f01*$m06 + 0.125*$f02*$m06 + 0.125*$f03*$m06 + 0.125*$f04*$m06 + 0.125*$f05*$m06 + 0.125*$f06*$m06 + 0.125*$f07*$m06 + 0.25*$f01*$m07 + 0.125*$f02*$m07 + 0.125*$f04*$m07 + 0.125*$f06*$m07 + 0.5*$f02*$m08 + 1*$f03*$m08 + 0.25*$f04*$m08 + 0.5*$f05*$m08 + 0.25*$f02*$m09 + 0.5*$f03*$m09 + 0.125*$f04*$m09 + 0.25*$f05*$m09 + 0.25*$f01*$m10 + 0.25*$f02*$m10 + 0.25*$f03*$m10 + 0.125*$f04*$m10 + 0.125*$f05*$m10 + 0.5*$f01*$m11 + 0.25*$f02*$m11 + 0.125*$f04*$m11,
		0.0625*$f04*$m02 + 0.125*$f05*$m02 + 0.125*$f06*$m02 + 0.25*$f07*$m02 + 0.125*$f04*$m03 + 0.25*$f05*$m03 + 0.25*$f06*$m03 + 0.5*$f07*$m03 + 0.0625*$f02*$m06 + 0.125*$f03*$m06 + 0.0625*$f04*$m06 + 0.125*$f05*$m06 + 0.0625*$f06*$m06 + 0.125*$f07*$m06 + 0.125*$f02*$m07 + 0.25*$f03*$m07 + 0.125*$f04*$m07 + 0.25*$f05*$m07 + 0.125*$f06*$m07 + 0.25*$f07*$m07 + 0.125*$f02*$m10 + 0.25*$f03*$m10 + 0.0625*$f04*$m10 + 0.125*$f05*$m10 + 0.25*$f02*$m11 + 0.5*$f03*$m11 + 0.125*$f04*$m11 + 0.25*$f05*$m11,
		0.125*$f04*$m04 + 0.25*$f05*$m04 + 0.25*$f06*$m04 + 0.5*$f07*$m04 + 0.0625*$f04*$m05 + 0.125*$f05*$m05 + 0.125*$f06*$m05 + 0.25*$f07*$m05 + 0.0625*$f04*$m06 + 0.0625*$f05*$m06 + 0.125*$f06*$m06 + 0.125*$f07*$m06 + 0.0625*$f04*$m07 + 0.125*$f06*$m07 + 0.25*$f04*$m08 + 0.5*$f05*$m08 + 0.5*$f06*$m08 + 1*$f07*$m08 + 0.125*$f04*$m09 + 0.25*$f05*$m09 + 0.25*$f06*$m09 + 0.5*$f07*$m09 + 0.125*$f04*$m10 + 0.125*$f05*$m10 + 0.25*$f06*$m10 + 0.25*$f07*$m10 + 0.125*$f04*$m11 + 0.25*$f06*$m11,
		0.03125*$f04*$m06 + 0.0625*$f05*$m06 + 0.0625*$f06*$m06 + 0.125*$f07*$m06 + 0.0625*$f04*$m07 + 0.125*$f05*$m07 + 0.125*$f06*$m07 + 0.25*$f07*$m07 + 0.0625*$f04*$m10 + 0.125*$f05*$m10 + 0.125*$f06*$m10 + 0.25*$f07*$m10 + 0.125*$f04*$m11 + 0.25*$f05*$m11 + 0.25*$f06*$m11 + 0.5*$f07*$m11,
		0.5*$f01*$m01 + 0.25*$f02*$m01 + 0.125*$f04*$m01 + 0.25*$f01*$m02 + 0.125*$f02*$m02 + 0.0625*$f04*$m02 + 0.25*$f01*$m05 + 0.125*$f02*$m05 + 0.0625*$f04*$m05 + 0.125*$f01*$m06 + 0.0625*$f02*$m06 + 0.03125*$f04*$m06,
		0.25*$f02*$m01 + 0.5*$f03*$m01 + 0.125*$f04*$m01 + 0.25*$f05*$m01 + 0.25*$f01*$m02 + 0.25*$f02*$m02 + 0.25*$f03*$m02 + 0.125*$f04*$m02 + 0.125*$f05*$m02 + 0.5*$f01*$m03 + 0.25*$f02*$m03 + 0.125*$f04*$m03 + 0.125*$f02*$m05 + 0.25*$f03*$m05 + 0.0625*$f04*$m05 + 0.125*$f05*$m05 + 0.125*$f01*$m06 + 0.125*$f02*$m06 + 0.125*$f03*$m06 + 0.0625*$f04*$m06 + 0.0625*$f05*$m06 + 0.25*$f01*$m07 + 0.125*$f02*$m07 + 0.0625*$f04*$m07,
		0.125*$f02*$m02 + 0.25*$f03*$m02 + 0.0625*$f04*$m02 + 0.125*$f05*$m02 + 0.25*$f02*$m03 + 0.5*$f03*$m03 + 0.125*$f04*$m03 + 0.25*$f05*$m03 + 0.0625*$f02*$m06 + 0.125*$f03*$m06 + 0.03125*$f04*$m06 + 0.0625*$f05*$m06 + 0.125*$f02*$m07 + 0.25*$f03*$m07 + 0.0625*$f04*$m07 + 0.125*$f05*$m07,
		0.125*$f04*$m01 + 0.25*$f06*$m01 + 0.0625*$f04*$m02 + 0.125*$f06*$m02 + 0.5*$f01*$m04 + 0.25*$f02*$m04 + 0.25*$f04*$m04 + 0.25*$f06*$m04 + 0.25*$f01*$m05 + 0.125*$f02*$m05 + 0.125*$f04*$m05 + 0.125*$f06*$m05 + 0.125*$f01*$m06 + 0.0625*$f02*$m06 + 0.0625*$f04*$m06 + 0.0625*$f06*$m06 + 1*$f01*$m08 + 0.5*$f02*$m08 + 0.25*$f04*$m08 + 0.5*$f01*$m09 + 0.25*$f02*$m09 + 0.125*$f04*$m09 + 0.25*$f01*$m10 + 0.125*$f02*$m10 + 0.0625*$f04*$m10,
		0.125*$f04*$m01 + 0.25*$f06*$m01 + 0.0625*$f04*$m02 + 0.125*$f06*$m02 + 0.25*$f01*$m05 + 0.125*$f02*$m05 + 0.125*$f04*$m05 + 0.125*$f06*$m05 + 0.125*$f01*$m06 + 0.0625*$f02*$m06 + 0.0625*$f04*$m06 + 0.0625*$f06*$m06 + 0.5*$f01*$m09 + 0.25*$f02*$m09 + 0.125*$f04*$m09 + 0.25*$f01*$m10 + 0.125*$f02*$m10 + 0.0625*$f04*$m10,
		0.125*$f04*$m01 + 0.25*$f05*$m01 + 0.25*$f06*$m01 + 0.5*$f07*$m01 + 0.125*$f04*$m02 + 0.125*$f05*$m02 + 0.25*$f06*$m02 + 0.25*$f07*$m02 + 0.125*$f04*$m03 + 0.25*$f06*$m03 + 0.125*$f02*$m05 + 0.25*$f03*$m05 + 0.125*$f04*$m05 + 0.25*$f05*$m05 + 0.125*$f06*$m05 + 0.25*$f07*$m05 + 0.125*$f01*$m06 + 0.125*$f02*$m06 + 0.125*$f03*$m06 + 0.125*$f04*$m06 + 0.125*$f05*$m06 + 0.125*$f06*$m06 + 0.125*$f07*$m06 + 0.25*$f01*$m07 + 0.125*$f02*$m07 + 0.125*$f04*$m07 + 0.125*$f06*$m07 + 0.25*$f02*$m09 + 0.5*$f03*$m09 + 0.125*$f04*$m09 + 0.25*$f05*$m09 + 0.25*$f01*$m10 + 0.25*$f02*$m10 + 0.25*$f03*$m10 + 0.125*$f04*$m10 + 0.125*$f05*$m10 + 0.5*$f01*$m11 + 0.25*$f02*$m11 + 0.125*$f04*$m11,
		0.0625*$f04*$m02 + 0.125*$f05*$m02 + 0.125*$f06*$m02 + 0.25*$f07*$m02 + 0.125*$f04*$m03 + 0.25*$f05*$m03 + 0.25*$f06*$m03 + 0.5*$f07*$m03 + 0.0625*$f02*$m06 + 0.125*$f03*$m06 + 0.0625*$f04*$m06 + 0.125*$f05*$m06 + 0.0625*$f06*$m06 + 0.125*$f07*$m06 + 0.125*$f02*$m07 + 0.25*$f03*$m07 + 0.125*$f04*$m07 + 0.25*$f05*$m07 + 0.125*$f06*$m07 + 0.25*$f07*$m07 + 0.125*$f02*$m10 + 0.25*$f03*$m10 + 0.0625*$f04*$m10 + 0.125*$f05*$m10 + 0.25*$f02*$m11 + 0.5*$f03*$m11 + 0.125*$f04*$m11 + 0.25*$f05*$m11,
		0.125*$f04*$m04 + 0.25*$f06*$m04 + 0.0625*$f04*$m05 + 0.125*$f06*$m05 + 0.03125*$f04*$m06 + 0.0625*$f06*$m06 + 0.25*$f04*$m08 + 0.5*$f06*$m08 + 0.125*$f04*$m09 + 0.25*$f06*$m09 + 0.0625*$f04*$m10 + 0.125*$f06*$m10,
		0.0625*$f04*$m05 + 0.125*$f06*$m05 + 0.03125*$f04*$m06 + 0.0625*$f06*$m06 + 0.125*$f04*$m09 + 0.25*$f06*$m09 + 0.0625*$f04*$m10 + 0.125*$f06*$m10,
		0.0625*$f04*$m05 + 0.125*$f05*$m05 + 0.125*$f06*$m05 + 0.25*$f07*$m05 + 0.0625*$f04*$m06 + 0.0625*$f05*$m06 + 0.125*$f06*$m06 + 0.125*$f07*$m06 + 0.0625*$f04*$m07 + 0.125*$f06*$m07 + 0.125*$f04*$m09 + 0.25*$f05*$m09 + 0.25*$f06*$m09 + 0.5*$f07*$m09 + 0.125*$f04*$m10 + 0.125*$f05*$m10 + 0.25*$f06*$m10 + 0.25*$f07*$m10 + 0.125*$f04*$m11 + 0.25*$f06*$m11,
		0.03125*$f04*$m06 + 0.0625*$f05*$m06 + 0.0625*$f06*$m06 + 0.125*$f07*$m06 + 0.0625*$f04*$m07 + 0.125*$f05*$m07 + 0.125*$f06*$m07 + 0.25*$f07*$m07 + 0.0625*$f04*$m10 + 0.125*$f05*$m10 + 0.125*$f06*$m10 + 0.25*$f07*$m10 + 0.125*$f04*$m11 + 0.25*$f05*$m11 + 0.25*$f06*$m11 + 0.5*$f07*$m11
	);
	
	# Calculate standardized genotype frequencies to sum to 1
	my @after_mating_std;
	for (my $i=0; $i<18; $i++) {
		$after_mating_std[$i] = $after_mating[$i] / sum(@after_mating);
	}

	return @after_mating_std;
}

# Function performs 1 generation of selection, followed by random mating
# Send function 2 arrays:
#	1. counts of the 18 genotypes
#	2. fitness values for the 18 genotypes
sub recursion_seln {
	my($raw_counts_ref, $fitness_ref) = @_; # pass the function references to the arrays
	
	# Dereference into arrays
	my @raw_counts = @{$raw_counts_ref}; 
	my @fitness = @{$fitness_ref};

	my @counts;
	
	# impose selection pressures
	for(my $i=0; $i<18; $i++) {
		$counts[$i] = $raw_counts[$i] * $fitness[$i];
	}
	
	# Calculate standardized genotype frequencies to sum to 1
	my @counts_norm;
	for (my $i=0; $i<18; $i++) {
		$counts_norm[$i] = $counts[$i] / sum(@counts);
	}	

	my $f01 = $counts_norm[0];
	my $f02 = $counts_norm[1];
	my $f03 = $counts_norm[2];
	my $f04 = $counts_norm[3];
	my $f05 = $counts_norm[4];
	my $f06 = $counts_norm[5];
	my $f07 = $counts_norm[6];

	my $m01 = $counts_norm[7];
	my $m02 = $counts_norm[8];
	my $m03 = $counts_norm[9];
	my $m04 = $counts_norm[10];
	my $m05 = $counts_norm[11];
	my $m06 = $counts_norm[12];
	my $m07 = $counts_norm[13];
	my $m08 = $counts_norm[14];
	my $m09 = $counts_norm[15];
	my $m10 = $counts_norm[16];
	my $m11 = $counts_norm[17]; 
	
	# Calculate genotype frequencies after random mating
	my @after_mating = (
		0.5*$f01*$m01 + 0.25*$f02*$m01 + 0.125*$f04*$m01 + 0.25*$f01*$m02 + 0.125*$f02*$m02 + 0.0625*$f04*$m02 + 0.5*$f01*$m04 + 0.25*$f02*$m04 + 0.125*$f04*$m04 + 0.25*$f01*$m05 + 0.125*$f02*$m05 + 0.0625*$f04*$m05 + 0.125*$f01*$m06 + 0.0625*$f02*$m06 + 0.03125*$f04*$m06,
		0.25*$f02*$m01 + 0.5*$f03*$m01 + 0.125*$f04*$m01 + 0.25*$f05*$m01 + 0.25*$f01*$m02 + 0.25*$f02*$m02 + 0.25*$f03*$m02 + 0.125*$f04*$m02 + 0.125*$f05*$m02 + 0.5*$f01*$m03 + 0.25*$f02*$m03 + 0.125*$f04*$m03 + 0.25*$f02*$m04 + 0.5*$f03*$m04 + 0.125*$f04*$m04 + 0.25*$f05*$m04 + 0.125*$f02*$m05 + 0.25*$f03*$m05 + 0.0625*$f04*$m05 + 0.125*$f05*$m05 + 0.125*$f01*$m06 + 0.125*$f02*$m06 + 0.125*$f03*$m06 + 0.0625*$f04*$m06 + 0.0625*$f05*$m06 + 0.25*$f01*$m07 + 0.125*$f02*$m07 + 0.0625*$f04*$m07,
		0.125*$f02*$m02 + 0.25*$f03*$m02 + 0.0625*$f04*$m02 + 0.125*$f05*$m02 + 0.25*$f02*$m03 + 0.5*$f03*$m03 + 0.125*$f04*$m03 + 0.25*$f05*$m03 + 0.0625*$f02*$m06 + 0.125*$f03*$m06 + 0.03125*$f04*$m06 + 0.0625*$f05*$m06 + 0.125*$f02*$m07 + 0.25*$f03*$m07 + 0.0625*$f04*$m07 + 0.125*$f05*$m07,
		0.125*$f04*$m01 + 0.25*$f05*$m01 + 0.25*$f06*$m01 + 0.5*$f07*$m01 + 0.125*$f04*$m02 + 0.125*$f05*$m02 + 0.25*$f06*$m02 + 0.25*$f07*$m02 + 0.125*$f04*$m03 + 0.25*$f06*$m03 + 0.25*$f02*$m04 + 0.5*$f03*$m04 + 0.25*$f04*$m04 + 0.5*$f05*$m04 + 0.25*$f06*$m04 + 0.5*$f07*$m04 + 0.125*$f02*$m05 + 0.25*$f03*$m05 + 0.125*$f04*$m05 + 0.25*$f05*$m05 + 0.125*$f06*$m05 + 0.25*$f07*$m05 + 0.125*$f01*$m06 + 0.125*$f02*$m06 + 0.125*$f03*$m06 + 0.125*$f04*$m06 + 0.125*$f05*$m06 + 0.125*$f06*$m06 + 0.125*$f07*$m06 + 0.25*$f01*$m07 + 0.125*$f02*$m07 + 0.125*$f04*$m07 + 0.125*$f06*$m07 + 0.5*$f02*$m08 + 1*$f03*$m08 + 0.25*$f04*$m08 + 0.5*$f05*$m08 + 0.25*$f02*$m09 + 0.5*$f03*$m09 + 0.125*$f04*$m09 + 0.25*$f05*$m09 + 0.25*$f01*$m10 + 0.25*$f02*$m10 + 0.25*$f03*$m10 + 0.125*$f04*$m10 + 0.125*$f05*$m10 + 0.5*$f01*$m11 + 0.25*$f02*$m11 + 0.125*$f04*$m11,
		0.0625*$f04*$m02 + 0.125*$f05*$m02 + 0.125*$f06*$m02 + 0.25*$f07*$m02 + 0.125*$f04*$m03 + 0.25*$f05*$m03 + 0.25*$f06*$m03 + 0.5*$f07*$m03 + 0.0625*$f02*$m06 + 0.125*$f03*$m06 + 0.0625*$f04*$m06 + 0.125*$f05*$m06 + 0.0625*$f06*$m06 + 0.125*$f07*$m06 + 0.125*$f02*$m07 + 0.25*$f03*$m07 + 0.125*$f04*$m07 + 0.25*$f05*$m07 + 0.125*$f06*$m07 + 0.25*$f07*$m07 + 0.125*$f02*$m10 + 0.25*$f03*$m10 + 0.0625*$f04*$m10 + 0.125*$f05*$m10 + 0.25*$f02*$m11 + 0.5*$f03*$m11 + 0.125*$f04*$m11 + 0.25*$f05*$m11,
		0.125*$f04*$m04 + 0.25*$f05*$m04 + 0.25*$f06*$m04 + 0.5*$f07*$m04 + 0.0625*$f04*$m05 + 0.125*$f05*$m05 + 0.125*$f06*$m05 + 0.25*$f07*$m05 + 0.0625*$f04*$m06 + 0.0625*$f05*$m06 + 0.125*$f06*$m06 + 0.125*$f07*$m06 + 0.0625*$f04*$m07 + 0.125*$f06*$m07 + 0.25*$f04*$m08 + 0.5*$f05*$m08 + 0.5*$f06*$m08 + 1*$f07*$m08 + 0.125*$f04*$m09 + 0.25*$f05*$m09 + 0.25*$f06*$m09 + 0.5*$f07*$m09 + 0.125*$f04*$m10 + 0.125*$f05*$m10 + 0.25*$f06*$m10 + 0.25*$f07*$m10 + 0.125*$f04*$m11 + 0.25*$f06*$m11,
		0.03125*$f04*$m06 + 0.0625*$f05*$m06 + 0.0625*$f06*$m06 + 0.125*$f07*$m06 + 0.0625*$f04*$m07 + 0.125*$f05*$m07 + 0.125*$f06*$m07 + 0.25*$f07*$m07 + 0.0625*$f04*$m10 + 0.125*$f05*$m10 + 0.125*$f06*$m10 + 0.25*$f07*$m10 + 0.125*$f04*$m11 + 0.25*$f05*$m11 + 0.25*$f06*$m11 + 0.5*$f07*$m11,
		0.5*$f01*$m01 + 0.25*$f02*$m01 + 0.125*$f04*$m01 + 0.25*$f01*$m02 + 0.125*$f02*$m02 + 0.0625*$f04*$m02 + 0.25*$f01*$m05 + 0.125*$f02*$m05 + 0.0625*$f04*$m05 + 0.125*$f01*$m06 + 0.0625*$f02*$m06 + 0.03125*$f04*$m06,
		0.25*$f02*$m01 + 0.5*$f03*$m01 + 0.125*$f04*$m01 + 0.25*$f05*$m01 + 0.25*$f01*$m02 + 0.25*$f02*$m02 + 0.25*$f03*$m02 + 0.125*$f04*$m02 + 0.125*$f05*$m02 + 0.5*$f01*$m03 + 0.25*$f02*$m03 + 0.125*$f04*$m03 + 0.125*$f02*$m05 + 0.25*$f03*$m05 + 0.0625*$f04*$m05 + 0.125*$f05*$m05 + 0.125*$f01*$m06 + 0.125*$f02*$m06 + 0.125*$f03*$m06 + 0.0625*$f04*$m06 + 0.0625*$f05*$m06 + 0.25*$f01*$m07 + 0.125*$f02*$m07 + 0.0625*$f04*$m07,
		0.125*$f02*$m02 + 0.25*$f03*$m02 + 0.0625*$f04*$m02 + 0.125*$f05*$m02 + 0.25*$f02*$m03 + 0.5*$f03*$m03 + 0.125*$f04*$m03 + 0.25*$f05*$m03 + 0.0625*$f02*$m06 + 0.125*$f03*$m06 + 0.03125*$f04*$m06 + 0.0625*$f05*$m06 + 0.125*$f02*$m07 + 0.25*$f03*$m07 + 0.0625*$f04*$m07 + 0.125*$f05*$m07,
		0.125*$f04*$m01 + 0.25*$f06*$m01 + 0.0625*$f04*$m02 + 0.125*$f06*$m02 + 0.5*$f01*$m04 + 0.25*$f02*$m04 + 0.25*$f04*$m04 + 0.25*$f06*$m04 + 0.25*$f01*$m05 + 0.125*$f02*$m05 + 0.125*$f04*$m05 + 0.125*$f06*$m05 + 0.125*$f01*$m06 + 0.0625*$f02*$m06 + 0.0625*$f04*$m06 + 0.0625*$f06*$m06 + 1*$f01*$m08 + 0.5*$f02*$m08 + 0.25*$f04*$m08 + 0.5*$f01*$m09 + 0.25*$f02*$m09 + 0.125*$f04*$m09 + 0.25*$f01*$m10 + 0.125*$f02*$m10 + 0.0625*$f04*$m10,
		0.125*$f04*$m01 + 0.25*$f06*$m01 + 0.0625*$f04*$m02 + 0.125*$f06*$m02 + 0.25*$f01*$m05 + 0.125*$f02*$m05 + 0.125*$f04*$m05 + 0.125*$f06*$m05 + 0.125*$f01*$m06 + 0.0625*$f02*$m06 + 0.0625*$f04*$m06 + 0.0625*$f06*$m06 + 0.5*$f01*$m09 + 0.25*$f02*$m09 + 0.125*$f04*$m09 + 0.25*$f01*$m10 + 0.125*$f02*$m10 + 0.0625*$f04*$m10,
		0.125*$f04*$m01 + 0.25*$f05*$m01 + 0.25*$f06*$m01 + 0.5*$f07*$m01 + 0.125*$f04*$m02 + 0.125*$f05*$m02 + 0.25*$f06*$m02 + 0.25*$f07*$m02 + 0.125*$f04*$m03 + 0.25*$f06*$m03 + 0.125*$f02*$m05 + 0.25*$f03*$m05 + 0.125*$f04*$m05 + 0.25*$f05*$m05 + 0.125*$f06*$m05 + 0.25*$f07*$m05 + 0.125*$f01*$m06 + 0.125*$f02*$m06 + 0.125*$f03*$m06 + 0.125*$f04*$m06 + 0.125*$f05*$m06 + 0.125*$f06*$m06 + 0.125*$f07*$m06 + 0.25*$f01*$m07 + 0.125*$f02*$m07 + 0.125*$f04*$m07 + 0.125*$f06*$m07 + 0.25*$f02*$m09 + 0.5*$f03*$m09 + 0.125*$f04*$m09 + 0.25*$f05*$m09 + 0.25*$f01*$m10 + 0.25*$f02*$m10 + 0.25*$f03*$m10 + 0.125*$f04*$m10 + 0.125*$f05*$m10 + 0.5*$f01*$m11 + 0.25*$f02*$m11 + 0.125*$f04*$m11,
		0.0625*$f04*$m02 + 0.125*$f05*$m02 + 0.125*$f06*$m02 + 0.25*$f07*$m02 + 0.125*$f04*$m03 + 0.25*$f05*$m03 + 0.25*$f06*$m03 + 0.5*$f07*$m03 + 0.0625*$f02*$m06 + 0.125*$f03*$m06 + 0.0625*$f04*$m06 + 0.125*$f05*$m06 + 0.0625*$f06*$m06 + 0.125*$f07*$m06 + 0.125*$f02*$m07 + 0.25*$f03*$m07 + 0.125*$f04*$m07 + 0.25*$f05*$m07 + 0.125*$f06*$m07 + 0.25*$f07*$m07 + 0.125*$f02*$m10 + 0.25*$f03*$m10 + 0.0625*$f04*$m10 + 0.125*$f05*$m10 + 0.25*$f02*$m11 + 0.5*$f03*$m11 + 0.125*$f04*$m11 + 0.25*$f05*$m11,
		0.125*$f04*$m04 + 0.25*$f06*$m04 + 0.0625*$f04*$m05 + 0.125*$f06*$m05 + 0.03125*$f04*$m06 + 0.0625*$f06*$m06 + 0.25*$f04*$m08 + 0.5*$f06*$m08 + 0.125*$f04*$m09 + 0.25*$f06*$m09 + 0.0625*$f04*$m10 + 0.125*$f06*$m10,
		0.0625*$f04*$m05 + 0.125*$f06*$m05 + 0.03125*$f04*$m06 + 0.0625*$f06*$m06 + 0.125*$f04*$m09 + 0.25*$f06*$m09 + 0.0625*$f04*$m10 + 0.125*$f06*$m10,
		0.0625*$f04*$m05 + 0.125*$f05*$m05 + 0.125*$f06*$m05 + 0.25*$f07*$m05 + 0.0625*$f04*$m06 + 0.0625*$f05*$m06 + 0.125*$f06*$m06 + 0.125*$f07*$m06 + 0.0625*$f04*$m07 + 0.125*$f06*$m07 + 0.125*$f04*$m09 + 0.25*$f05*$m09 + 0.25*$f06*$m09 + 0.5*$f07*$m09 + 0.125*$f04*$m10 + 0.125*$f05*$m10 + 0.25*$f06*$m10 + 0.25*$f07*$m10 + 0.125*$f04*$m11 + 0.25*$f06*$m11,
		0.03125*$f04*$m06 + 0.0625*$f05*$m06 + 0.0625*$f06*$m06 + 0.125*$f07*$m06 + 0.0625*$f04*$m07 + 0.125*$f05*$m07 + 0.125*$f06*$m07 + 0.25*$f07*$m07 + 0.0625*$f04*$m10 + 0.125*$f05*$m10 + 0.125*$f06*$m10 + 0.25*$f07*$m10 + 0.125*$f04*$m11 + 0.25*$f05*$m11 + 0.25*$f06*$m11 + 0.5*$f07*$m11
	);
	
	# Calculate standardized genotype frequencies to sum to 1
	my @after_mating_std;
	for (my $i=0; $i<18; $i++) {
		$after_mating_std[$i] = $after_mating[$i] / sum(@after_mating);
	}

	return @after_mating_std;
}


# Function calculates the allele frequencies 
# given an array of all possible genotypes
sub allele_freqs{
	my @freqs = @_;

#0 female	f01	X1X1	ZZ	X2X2
#1 female	f02	X1X1	ZW	X2X2
#2 female	f03	X1X1	WW	X2X2
#3 female	f04	X1Y1	ZW	X2X2
#4 female	f05	X1Y1	WW	X2X2
#5 female	f06	Y1Y1	ZW	X2X2
#6 female	f07	Y1Y1	WW	X2X2
#7 male		m01	X1X1	ZZ	X2Y2
#8 male		m02	X1X1	ZW	X2Y2
#9 male		m03	X1X1	WW	X2Y2
#10 male	m04	X1Y1	ZZ	X2X2
#11 male	m05	X1Y1	ZZ	X2Y2
#12 male	m06	X1Y1	ZW	X2Y2
#13 male	m07	X1Y1	WW	X2Y2
#14 male	m08	Y1Y1	ZZ	X2X2
#15 male	m09	Y1Y1	ZZ	X2Y2
#16 male	m10	Y1Y1	ZW	X2Y2
#17 male	m11	Y1Y1	WW	X2Y2

# print OUTPUT "X1\tY1\tX2\tY2\tZ\tW\t"; # final chromosome freqs
	
	my $X1_freq = ($freqs[0] + $freqs[1] + $freqs[2] + $freqs[7] + $freqs[8] + $freqs[9]
				+ 0.5*($freqs[3] + $freqs[4] + $freqs[10] + $freqs[11] + $freqs[12] + $freqs[13])) / sum(@freqs);
	my $Y1_freq = 1 - $X1_freq;

	my $X2_freq = ($freqs[0] + $freqs[1] + $freqs[2] + $freqs[3] + $freqs[4] + $freqs[5] + $freqs[6] + $freqs[10] + $freqs[14]
				+ 0.5*($freqs[7] + $freqs[8] + $freqs[9] + $freqs[11] + $freqs[12] + $freqs[13] + $freqs[15] + $freqs[16] + $freqs[17])) / sum(@freqs);
	my $Y2_freq = 1 - $X2_freq;
	
	my $Z_freq = ($freqs[0] + $freqs[7] + $freqs[10] + $freqs[11] + $freqs[14] + $freqs[15] 
				+ 0.5*($freqs[1] + $freqs[3] + $freqs[5] + $freqs[8] + $freqs[12] + $freqs[16])) / sum(@freqs);
	my $W_freq = 1 - $Z_freq;

	return ($X1_freq, $Y1_freq, $X2_freq, $Y2_freq, $Z_freq, $W_freq);
}


# Function determines the genotype fitness from allele fitness values based on: 
# 1. array with 5 fitness values
# 2. array with dominance effects
sub genotype_fitness{
	my($allele_ref, $dom_ref) = @_; # pass the function references to the arrays

	
	# Dereference into arrays
	my @allele_fits = @{$allele_ref}; 
	my @doms = @{$dom_ref};

	#	0 Y1 	female
	#	1 W		female
	#	2 Y1	male
	#	3 Y2   	male
	#	4 W		male	
	my($Y1_fem, $W_fem, $Y1_male, $Y2_male, $W_male) = @allele_fits;

	# 0 female Y dominance
	# 1 male Y dominance
	# 2 female W dominance
	# 3 male W dominance
	my($fem_Y_dom, $male_Y_dom, $fem_W_dom, $male_W_dom) = @doms;
	
	## Calculate single locus fitness effects of genotypes:	
	
	# Y1 female
	my $X1X1_fem; 
	my $X1Y1_fem; 
	my $Y1Y1_fem; 
	if($Y1_fem < 0){ #if Y1 is deleterious to females
		$X1X1_fem  = 1;
		$X1Y1_fem = 1 + $fem_Y_dom*$Y1_fem;
		$Y1Y1_fem = 1 + $Y1_fem;
	}else{ # if Y1 is beneficial to females
		$X1X1_fem  = 1 - $Y1_fem;
		$X1Y1_fem = 1 - $fem_Y_dom*$Y1_fem;
		$Y1Y1_fem = 1; 
	}
		
	# W female
	my $ZZ_fem; 
	my $ZW_fem; 
	my $WW_fem;
	if($W_fem < 0){ #if W is deleterious to females
		$ZZ_fem  = 1;
		$ZW_fem = 1 + $fem_W_dom*$W_fem;
		$WW_fem = 1 + $W_fem;
	}else{ # if W is beneficial to females
		$ZZ_fem  = 1 - $W_fem;
		$ZW_fem = 1 - $fem_W_dom*$W_fem;
		$WW_fem = 1; 
	}
	
	# Y1 male
	my $X1X1_male; 
	my $X1Y1_male; 
	my $Y1Y1_male; 
	if($Y1_male < 0){ #if Y1 is deleterious to males
		$X1X1_male  = 1;
		$X1Y1_male = 1 + $male_Y_dom*$Y1_male;
		$Y1Y1_male = 1 + $Y1_male;
	}else{ # if Y1 is beneficial to males
		$X1X1_male  = 1 - $Y1_male;
		$X1Y1_male = 1 - $male_Y_dom*$Y1_male;
		$Y1Y1_male = 1; 
	}

	# Y2 male
	my $X2X2_male; 
	my $X2Y2_male; 
	if($Y2_male < 0){ #if Y2 is deleterious to males
		$X2X2_male = 1;
		$X2Y2_male = 1 + $Y2_male;
	}else{ # if Y2 is beneficial to males
		$X2X2_male = 1;
		$X2Y2_male = 1 - $Y2_male;
	}

	# W male
	my $ZZ_male; 
	my $ZW_male; 
	my $WW_male;
	if($W_male < 0){ #if W is deleterious to males
		$ZZ_male  = 1;
		$ZW_male = 1 + $male_W_dom*$W_male;
		$WW_male = 1 + $W_male;
	}else{ # if W is beneficial to males
		$ZZ_male  = 1 - $W_male;
		$ZW_male = 1 - $male_W_dom*$W_male;
		$WW_male = 1; 
	}	
	
	# Calculate the genotype fitness 
	my @genotype_fit = (
		$X1X1_fem * $ZZ_fem ,					#0  female	f01	X1X1	ZZ	X2X2
		$X1X1_fem * $ZW_fem ,					#1  female	f02	X1X1	ZW	X2X2
		$X1X1_fem * $WW_fem ,					#2  female	f03	X1X1	WW	X2X2
		$Y1Y1_fem * $ZW_fem ,					#3  female	f04	X1Y1	ZW	X2X2
		$X1Y1_fem * $WW_fem , 					#4  female	f05	X1Y1	WW	X2X2
		$Y1Y1_fem * $ZW_fem ,					#5  female	f06	Y1Y1	ZW	X2X2
		$Y1Y1_fem * $WW_fem , 					#6  female	f07	Y1Y1	WW	X2X2
		$X1X1_male * $ZZ_male * $X2Y2_male ,	#7  male	m01	X1X1	ZZ	X2Y2
		$X1X1_male * $ZW_male * $X2Y2_male ,	#8  male	m02	X1X1	ZW	X2Y2
		$X1X1_male * $WW_male * $X2Y2_male ,	#9  male	m03	X1X1	WW	X2Y2
		$Y1Y1_male * $ZZ_male * $X2X2_male ,	#10 male	m04	X1Y1	ZZ	X2X2
		$X1Y1_male * $ZZ_male * $X2Y2_male ,	#11 male	m05	X1Y1	ZZ	X2Y2
		$X1Y1_male * $ZW_male * $X2Y2_male ,	#12 male	m06	X1Y1	ZW	X2Y2
		$X1Y1_male * $WW_male * $X2Y2_male ,	#13 male	m07	X1Y1	WW	X2Y2
		$Y1Y1_male * $ZZ_male * $X2X2_male ,	#14 male	m08	Y1Y1	ZZ	X2X2
		$Y1Y1_male * $ZZ_male * $X2Y2_male ,	#15 male	m09	Y1Y1	ZZ	X2Y2
		$Y1Y1_male * $ZW_male * $X2Y2_male ,	#16 male	m10	Y1Y1	ZW	X2Y2
		$Y1Y1_male * $WW_male * $X2Y2_male		#17 male	m11	Y1Y1	WW	X2Y2		
	);

	my @genotype_norm = @genotype_fit;
	
	# Normalize male fitness to ensure that max = 1
	for(my $i=7; $i < 18; $i++){
		$genotype_norm[$i] = $genotype_fit[$i] / max(@genotype_fit[7..17]);
	}

	# Normalize female fitness to ensure that max = 1
	for(my $i=0; $i < 7; $i++){
		$genotype_norm[$i] = $genotype_fit[$i] / max(@genotype_fit[0..6]);
	}

	# return normalize genotype fitness values
	return(@genotype_norm);
}



