# This file calculates the relative weights of health gains to the most and least healthy fifths 
# of the population based on inequality aversion parameters found from the citizens panel.
#
# Author: Miqdad Asaria
# Date: May 2015
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
# solve numerically by finding the health sacrifice in the most affluent which equates to the social 
# welfare of a health sacrifice of 1 QALY in the most deprived
#####################################################################################################

baseline_health_distribution = c(62,68,70,72,74)
atkinson_e = 10.946
kolm_alpha = 0.152

atkinson_ratio = function(increment){
  ratio = uniroot(function(x,e,baseline) calculate_atkinson_ede(c(baseline[1]+increment,baseline[c(2:5)]),e) - calculate_atkinson_ede(c(baseline[c(1:4)],c(baseline[5])+x),e),
                     interval=c(-20,20),
                     e=atkinson_e,
                     tol=abs(increment)*.Machine$double.eps^0.25,
                     baseline=baseline_health_distribution)
  return(ratio$root)
}

kolm_ratio = function(increment){
  ratio = uniroot(function(x,alpha,baseline) calculate_kolm_ede(c(baseline[1]+increment,baseline[c(2:5)]),alpha) - calculate_kolm_ede(c(baseline[c(1:4)],c(baseline[5])+x),alpha),
                       interval=c(-20,20),
                       alpha=kolm_alpha,
                       tol=abs(increment)*.Machine$double.eps^0.25,
                       baseline=baseline_health_distribution)
  return(ratio$root)
}

print_relative_weights = function(increment){
  cat(paste("Using an increment of ",increment," QALYs",sep=""))
  cat(paste("Atkinson Social Welfare Index with inequality aversion level ", atkinson_e, " implies investment in the health of the most deprived fifth of the population should be given ", round(atkinson_ratio(increment)*1/increment,2) ," times as much priority as investment in the health of the most affluent fifth of the population.",sep=""))
  cat(paste("Kolm Social Welfare Index with inequality aversion level ", kolm_alpha, " implies investment in the health of the most deprived fifth of the population should be given ", round(kolm_ratio(increment)*1/increment,2) ," times as much priority as investment in the health of the most affluent fifth of the population.",sep=""))  
  cat("\n\n")
}
#####################################################################################################
# print out the implied weights given by the two social welfare functions
#####################################################################################################
sink("relative_weights.txt")
print_relative_weights(increment=1)
print_relative_weights(increment=-1)

print_relative_weights(increment=0.1)
print_relative_weights(increment=-0.1)

print_relative_weights(increment=0.01)
print_relative_weights(increment=-0.01)

print_relative_weights(increment=0.001)
print_relative_weights(increment=-0.001)

print_relative_weights(increment=0.0001)
print_relative_weights(increment=-0.0001)

print_relative_weights(increment=0.00001)
print_relative_weights(increment=-0.00001)
sink()