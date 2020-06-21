@doc doc"""

Load former results to compute graphs
"""

struct resps
    f			:: Array{Float64,4}
	ratio_alpha	:: Int64
	logPost  	:: Array{Float64,1}
	eps_big     :: Array{Float64,3}
	HDstore     :: Array{Float64,4}
end

function load_results(options_set :: Options)

	options_set = get_draws!(options_set)

	unpack(options_set.paths)
	unpack(options_set.options_var)
	unpack(options_set.options_sr)

Response    = resps(zeros(nvar,nvar,horizon,length(DrawDone)*alpha),0,NaN*Array{Float64}(undef,length(DrawDone)*alpha),NaN*Array{Float64}(undef,nvar,T-lags,length(DrawDone)*alpha),NaN*Array{Float64,4}(undef,T,nvar,nvar,length(DrawDone)));
HDstore 	= Array{Float64,4}(undef,T,nvar,nvar,length(DrawDone))
resp_low 	= Array{Float64,3}(undef,nvar,nvar,horizon)
resp_high 	= Array{Float64,3}(undef,nvar,nvar,horizon)
resp_50 	= Array{Float64,3}(undef,nvar,nvar,horizon)
resp_mean	= Array{Float64,3}(undef,nvar,nvar,horizon)

patternFile = r"Draw_(.+)_(.+).jld";
fileirf 	= readdir(join(save_path_IRF));

print("Loading draws:    ")
if isempty(options_set.options_sr["DrawDone"]) 
    print("No Draws Done")
    filename 		                          = join(save_path_IRF)*"Imprest.jld";
    results_mat                               = load(filename)

    for Shock in 1 : results_mat["nvar"]
        for Reac in 1 : results_mat["nvar"]
            for i in 1 : results_mat["horizon"]
                resp_low[Reac,Shock,i]   =   percentile(results_mat["Response"].f[Reac,Shock,i,:],low_band);
                resp_high[Reac,Shock,i]  =   percentile(results_mat["Response"].f[Reac,Shock,i,:],high_band);
                resp_50[Reac,Shock,i]    =   percentile(results_mat["Response"].f[Reac,Shock,i,:],50);
                resp_mean[Reac,Shock,i]  =   mean(results_mat["Response"].f[Reac,Shock,i,:]);
            end # for i
        end # for ReacSh
    end # for Shock 
    
# Variance Decomposition
HDM = (median(results_mat["HDstore"],dims=4));

# Median Structural Shocks
SHK = (results_mat["Response"].eps_big);

# Variance Decomposition
VD = (vardec_sr_new(resp_50));

# IRFs
IRF = Dict("resp_low"=>resp_low,"resp_high"=>resp_high,"resp_50"=>resp_50,"resp_mean"=>resp_mean)


else  

for sim in 1 : length(DrawDone)
    print([repeat("\b",length(string.(sim))) * string.(sim)][1])

    filename 		                            = join(save_path_IRF)*fileirf[sim];
    OutputLoad                                  = load(filename)["Draw"]
    
    Response.logPost[(sim-1)*alpha+1]           = OutputLoad["logPost"];
    Response.eps_big[:,:,(sim-1)*alpha+1]       = OutputLoad["Eps"];
    Response.HDstore[:,:,:,sim]                 = OutputLoad["HD_save_Shock"];
    Response.f[:,:,:,(sim-1)*alpha+1:sim*alpha] = OutputLoad["IRF"][:,:,1:horizon,1];


end


for Shock in 1 : nvar
    for Reac in 1 : nvar
        for i in 1 : horizon
            resp_low[Reac,Shock,i]   =   percentile(Response.f[Reac,Shock,i,:],low_band);
            resp_high[Reac,Shock,i]  =   percentile(Response.f[Reac,Shock,i,:],high_band);
            resp_50[Reac,Shock,i]    =   percentile(Response.f[Reac,Shock,i,:],50);
            resp_mean[Reac,Shock,i]  =   mean(Response.f[Reac,Shock,i,:]);
        end # for i
    end # for ReacSh
end # for Shock 
draws = options_set.options_sr["draws"]*alpha;

# Variance Decomposition
HDM = (median(Response.HDstore,dims=4));

# Median Structural Shocks
SHK = (Response.eps_big);

# Variance Decomposition
VD = (vardec_sr_new(resp_50));

# IRFs
IRF = Dict("resp_low"=>resp_low,"resp_high"=>resp_high,"resp_50"=>resp_50,"resp_mean"=>resp_mean)

# Save final output in Imprest
# JLD format
save(join(save_path_IRF)*"Imprest.jld","Response",Response,"VD",VD,"horizon",horizon,"nvar",nvar,"draws",draws,"HDM",HDM,"eps_big",SHK)

# Matlab format
matwrite(join(save_path_IRF)*"Imprest.mat",Dict("Response"=>Response,"VD"=>VD,"horizon"=>horizon,"nvar"=>nvar,"draws"=>draws,"HDM"=>HDM,"eps_big"=>SHK))

end

results_mat = Dict("HDM"=>HDM,"VD"=>VD,"IRF"=>IRF,"Response"=>Response,"SHK"=>SHK)

return results_mat
end



function vardec_sr_new(Response_burn)
#       vardec2 takes as input structural impulse Responses
#       from the Uhlig sign-restriction identification
#       and the associated Choleski impulse Responses
#       and computes the fraction of the variance attributable to the structural shock

(nvar,nshocks,nstep) = size(Response_burn);

forvar       = zeros(nvar,nshocks,nstep-1);
forvar1      = zeros(nvar,nshocks,nstep-1);
forvarm      = zeros(nvar,nshocks,nstep-1);
forvarm1     = zeros(nvar,nshocks,nstep-1);
diag_forvarm = zeros(nvar,nstep-1);

for h in 1:nstep-1
   forvar1[:,:,h]    = Response_burn[:,:,h].^2;
   forvar[:,:,h]     = sum(forvar1,dims=3);
   forvarm1[:,:,h]   = Response_burn[:,:,h]*Response_burn[:,:,h]';
   forvarm[:,:,h]    = sum(forvarm1,dims=3);
   diag_forvarm[:,h] = diag(forvarm[:,:,h]);
end

vardec = Array{Float64,3}(undef,nvar,nshocks,nstep-1);
#for each nstep j, vardec(:,:,j) has size nvar, nshocks
#in row i, vardec(:,:,j) decomposes forecast error variance in variable i
#into fraction arising from shocks columnwise
for j in 1 : nstep-1
    for m in 1 : nshocks
        for n in 1 : nvar
        vardec[n,m,j] = forvar[n,m,j] / diag_forvarm[n,j];
        end
    end
end

return vardec
end
