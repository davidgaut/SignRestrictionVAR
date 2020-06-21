@doc doc"""

Compute optimal lag criterion based on AIC, HQIC and SIC.
Takes arguments TS::TimeArray and options_set::Options and fills options_set with optimal criterion.
"""


function lag_crit!(TS :: TimeArray, options_set :: Options)

data  = values(TS);
(T,N) = size(data);

# Lag Criteria
lagmax  = 10;

AICC  = HQIC = SICC = LR = pval = zeros(lagmax + 1,1);

s2_d  = s2  = zeros(N,N)
s2_v  = zeros(N,N,lagmax+1)
s2_det= zeros(lagmax+1)

for λ in 0 : lagmax

	options_set.options_var["lags"] = λ

# Generate lag indiferent dataset
	yLagSelect = data[lagmax-λ+1:end,:];
	(Tt,K)     = size(data[lagmax+1:end,:]);

# Posterior with mean and variance parameters corresponding to OLS estimates
	var_results = var_OLS(TS,options_set.options_var);

	# for ii in 1:length(var_results); eval(:($(Symbol(var_results[ii])) = $var_results[ii])); end
     u = var_results["Res"]
# Calculate Residual and bfo for each draw
	s2_d  		  = u'*u./Tt;  # Estimated error variance
    s2_v[:,:,λ+1] = (s2_d);
    s2_det[λ+1]   = det(s2_d);

# Kilian and Lutkepohl (p55)
	regNbr 	      = (λ * K^2 + K);
    AICC[lags+1]  = log(s2_det[λ+1]) + regNbr * 2 / Tt;
    HQIC[lags+1]  = log(s2_det[λ+1]) + regNbr * 2 / Tt * log(log(Tt));
    SICC[lags+1]  = log(s2_det[λ+1]) + regNbr * 1 / Tt * log(Tt);

end

Crit     = [argmin(AICC)[1] argmin(HQIC)[1] argmin(SICC)[1]];
LagTable = DataFrame(Crit,[:AICC, :HQIC, :SICC])

println(LagTable)

if mode(Crit) != 1
Crit_Val = mode(Crit)
else
Crit_Val = minimum(Crit);
end

options_set.options_var["lags"] = Crit_Val
options_set.options_var["T"]    = T
println("\nCriterion Selection:  ", string(Crit_Val))

return options_set

end
