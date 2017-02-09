# This file illustrates how results from the citizens panel questionnaires can be used to 
# calculate inequaltity aversion parameters for the common constant relative and absolute 
# inequality aversion style social welfare functions the Atkinson and Kolm funtions.
#
# Given the health distributions that respondents have classed as being equally good this code
# numerically solves these social welfare functions to reveal the level of inequality aversion
# implied by the respondents choice.
#
# Author: Miqdad Asaria
# Date: August 2014
#####################################################################################################


#####################################################################################################
# define Atkinson and Kolm social welfare functions that when given a health distribution and a 
# level of inequality aversion return the equally distributed equivalent level of health
#####################################################################################################

# this function calculates the EDE of a health distribution defined in programme
# using the atkinson index using an inequality aversion level given by e
calculate_atkinson_ede = function(programme, e){
	 mean(programme)*mean((programme/mean(programme))^(1-e))^(1/(1-e))
}

# this function calculates the EDE of a health distribution defined in programme
# using the kolm index using an inequality aversion level given by alpha
calculate_kolm_ede = function(programme, alpha){
	mean(programme)-(1/alpha)*log(mean(exp(alpha*(mean(programme)-programme))))
}

#####################################################################################################
# define health distributions in terms of the health achieved by each quintile in the population
# under programmes A and B
#####################################################################################################

programme_a = c(65,68,70,72,81)
programme_b = c(70,68,70,72,77)


#####################################################################################################
# solve numerically by finding the inequality aversion level where the difference in EDE between the 
# health distributions defined by the programmes are zero
######################################################################################################

atkinson_e = uniroot(function(e, prog_a, prog_b) calculate_atkinson_ede(prog_a, e) - calculate_atkinson_ede(prog_b, e),
		interval=c(-50,50),
		prog_a=programme_a,
		prog_b=programme_b)

kolm_alpha = uniroot(function(alpha, prog_a, prog_b) calculate_kolm_ede(prog_a, alpha) - calculate_kolm_ede(prog_b, alpha),
		interval=c(-50,50),
		prog_a=programme_a,
		prog_b=programme_b)

#####################################################################################################
# print out the implied inequality aversion parameters for the two social welfare functions
#####################################################################################################

print(paste("Atkinson Inequality Aversion Level (e) making A and B equivalent is: ", round(atkinson_e$root,5),sep=""))
print(paste("Kolm Inequality Aversion Level (alpha) making A and B equivalent is: ", round(kolm_alpha$root,5),sep=""))