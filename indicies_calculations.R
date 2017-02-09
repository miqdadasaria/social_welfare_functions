###########################################################################
# This code creates stylised examples demonstrating the use of stochastic
# dominance of Lorenz and Generalised Lorenz curves to rank health policies
# in terms of their effect on inequalities in the univariate health case. It
# also calculates Atkinson indices for the various examples.
#
# PHRC Inequalities Project
# CHE - University of York
# Miqdad Asaria 
# August 2011
###########################################################################

# calculate the atkinson index for a health distribution 
# using inequality aversion parameter e
calculate_atkinson_index = function(health,e){
	mean_health = mean(health)
	n = length(health)
	health_ede = NULL
	if(e==1){
		health_ede = prod(health^(1/n))
	} else {
		health_ede = (sum(health^(1-e))/n)^(1/(1-e))
	}
	atkinson = 1-(health_ede/mean_health)
	return(atkinson)
}

# outputs the atkinson index across a range of inequality aversion values
atkinson_index = function(health_baseline, health_intervention, baseline_label, intervention_label, N=10, increment = 0.25){
	atkinson = array(NA, c(N,5),list(1:N,c("e","index_baseline","index_intervention","ede_baseline","ede_intervention")))
	for(i in 1:N){
		e = i*increment
		atkinson_baseline = calculate_atkinson_index(health_baseline,e)
		ede_baseline = mean(health_baseline) *(1-atkinson_baseline)
		atkinson_intervention = calculate_atkinson_index(health_intervention,e)
		ede_intervention = mean(health_intervention)*(1-atkinson_intervention)
		atkinson[i,] = c(e,atkinson_baseline, atkinson_intervention, ede_baseline, ede_intervention)
	}
	return(atkinson)
}

# calculated kolm index for a health vector with alpha as the inequality aversion parameter
calculate_kolm_index <- function(health, alpha)
{
	kolm <- (1/alpha)*log(mean(exp(alpha * (mean(health)-health))))
	return(kolm)
}

# outputs the atkinson index across a range of inequality aversion values
kolm_index = function(health_baseline, health_intervention, baseline_label, intervention_label, N=10, increment = 0.0125){
	kolm = array(NA, c(N,5),list(1:N,c("alpha","index_baseline","index_intervention","ede_baseline","ede_intervention")))
	for(i in 1:N){
		k = i*increment
		kolm_baseline = calculate_kolm_index(health_baseline,k)
		ede_baseline = (kolm_baseline * -1) + mean(health_baseline)
		kolm_intervention = calculate_kolm_index(health_intervention,k)
		ede_intervention = (kolm_intervention * -1) + mean(health_intervention)
		kolm[i,] = c(k,kolm_baseline,kolm_intervention,ede_baseline,ede_intervention)
	}
	return(kolm)
}

# calculate slope index
calculate_slope_index = function(health, cum_population){
	res = lm(health ~ cum_population[-1])
	return(res$coefficients[2])
}

# calculate 20:20 ratio index
calculate_ratio_index = function(health, health_frequencies){
	quintiles = get_quantiles(health, health_frequencies)
	return((quintiles[5,"health"]/quintiles[1,"health"])-1)
}

# calculate 20:20 gap index
calculate_gap_index = function(health, health_frequencies){
	quintiles = get_quantiles(health, health_frequencies)
	return(quintiles[5,"health"]-quintiles[1,"health"])
}

# calculate gini index assuming individual level data (Angus Deaton Formula 1997)
calculate_gini_index = function(health){
	health = sort(health, decreasing = TRUE)
	N = length(health)
	gini = ((N+1)/(N-1)) - ((2*sum(health * 1:N))/(N*(N-1)*mean(health)))
	return(gini)
}

# calculates the concentration index for individual level data v is inequality aversion parameter
# v=2 gives standard concentrarion index
calculate_concentration_index = function(health, v=2){
	rank = cumulative(rep(1,length(health)))[-1]
	conc = 1 - (v/(length(health)*mean(health))*sum(health*(1-rank)^(v-1)))
	return(conc)
}


# poverty measures
calculate_headcount_ratio = function(health, fair_innings){
	health_poor = health[health<fair_innings]
	H = 0
	if(length(health_poor)>1){
		H = sum(health<fair_innings)/length(health)
	}
	return(H)
}

calculate_health_gap_ratio = function(health, fair_innings){
	health_poor = health[health<fair_innings]
	I = 0
	if(length(health_poor)>1){
		I = sum((fair_innings-health_poor)/fair_innings)/length(health_poor)
	}
	return(I)
}

calculate_sen_index = function(health, fair_innings){
	health_poor = health[health<fair_innings]
	S = 0
	if(length(health_poor)>1){
		H = calculate_headcount_ratio(health, fair_innings)
		I = calculate_health_gap_ratio(health, fair_innings)
		G = calculate_gini_index(health_poor)
		S = H*(I+(1-I)*G)
	}
	return(S)
}

