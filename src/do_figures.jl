@doc doc"""
Generate figures

Inputs:
options
save_path: location of the folder
save_name: name     of the folder

Outputs:
return Options with paths
"""


function do_figures(results_mat :: Dict ,TS :: TimeArray, options_set :: Options)

    PyCall.PyDict(matplotlib.rcParams)["font.family"] = "STIXGeneral"
    PyCall.PyDict(matplotlib.rcParams)["font.size"]   = 12

    options_set.plot["hatch_plot"] ? nothing : (options_set.plot["hatch_map"] = split(" /"^length(options_set.plot["hatch_map"]),"/")) 

    
if options_set.plot["HD"] & options_set.plot["plot"]
    HSD_plot(results_mat["HDM"],TS,options_set);
    HSD_plot_OneShock(results_mat["HDM"],TS,options_set);
end
if options_set.plot["SHK"] & options_set.plot["plot"]
    SHK_plot(results_mat["SHK"],TS,options_set)
end
if options_set.plot["IRF"] & options_set.plot["plot"]
    IRF_plot(results_mat["IRF"],options_set)
end
if options_set.plot["VD"] & options_set.plot["plot"]
    idx_HD = Array{Int64,2}(undef,options_set.options_var["nvar"],options_set.options_var["nvar"])
    for ii in 1 : options_set.options_var["nvar"]
        idx_HD[ii,:] = sortperm(dropdims(sum(abs.(results_mat["HDM"][options_set.options_var["lags"]+1:end,:,ii]),dims=1),dims=1),rev=true);
    end
    VD_plot(results_mat["VD"],idx_HD,options_set)
end
if options_set.plot["MM"] & options_set.plot["plot"]
    IRF_modal_plot(results_mat["Response"],options_set)
end

end

function HSD_plot(HDM :: Array, TS :: TimeArray, options_set :: Options)
    TS = TS[1:133]
    HDM = HDM[1:133,:,:]
    nvar   = nshocks = options_set.options_var["nvar"]
    lags   =           options_set.options_var["lags"]
    HD     = Array{Float64}(undef,size(HDM,1),size(HDM,2),size(HDM,3))
    idx_HD = Array{Int64,2}(undef,nvar,nvar)
    
    for ii in 1 : nvar
        idx_HD[ii,:] = sortperm(dropdims(sum(abs.(HDM[lags+1:end,:,ii]),dims=1),dims=1),rev=true);
        idx_HD[ii,:] = [idx_HD[ii,findall(idx_HD[ii,:].!=6)] ; 6];
        idx_HD[ii,:] = [5;idx_HD[ii,findall(idx_HD[ii,:].!=5)]];
        HD[:,:,ii]   = HDM[:,idx_HD[1,:],ii];
    end

    #  Pick Series to filter
    idx_filter =  indexin(options_set.plot["filter"], options_set.options_var["series"])

    # Replace NaN with zero
    replace_nan(v) = map(x -> isnan(x) ? zero(x) : x, v)
    HD = replace_nan(HD);

    # Annual Growth
    if options_set.plot["annnualize"]
        fun2 = var  ->  var[4:end,:] + var[3:(end-1),:] + var[2:(end-2),:] + var[1:(end-3),:]
    else
        fun2 = var  ->  var
    end

    (nsteps,)  = size(HD);
    filename   = options_set.paths["save_path_FIG"]*"HSD_";

    sizbar 	   = 50
    snamesPlot = options_set.options_sr["shocks"]
    vnamesPlot = colnames(TS)
    xpos 	   = Array{Any,1}(undef,size(HD,3));
    xneg 	   = Array{Any,1}(undef,size(HD,3));

    for ii in 1 : nvar
        if in(ii,idx_filter)
            HD_loop = fun2(diff(HD[:,:,ii],dims=1)); dates = timestamp(TS)[5:end];
        else
            HD_loop = HD[:,:,ii]; dates = timestamp(TS)[1:end];
        end

        fig,axesData  = subplots()
        for ij in idx_HD[1,:]
            xpos[ij] = axesData.bar(dates,Int.(HD_loop[:,ij].>0).*HD_loop[:,ij],sizbar,facecolor = options_set.plot["color_map"][idx_HD[1,:],:][ij,:],
            bottom = (ij > 1 ? dropdims(sum(HD_loop[:,1:ij-1].*Int.(HD_loop[:,1:ij-1].>0),dims=2),dims=2) : 0), hatch = options_set.plot["hatch_map"][idx_HD[1,:],:][ij])
            xpos[ij] = axesData.bar(dates,Int.(HD_loop[:,ij].<0).*HD_loop[:,ij],sizbar,facecolor = options_set.plot["color_map"][idx_HD[1,:],:][ij,:],
            bottom = (ij > 1 ? dropdims(sum(HD_loop[:,1:ij-1].*Int.(HD_loop[:,1:ij-1].<0),dims=2),dims=2) : 0), hatch = options_set.plot["hatch_map"][idx_HD[1,:],:][ij])
        end        

        axesData.legend(xpos[1:end],snamesPlot[idx_HD[1,:]],labelspacing=0.05)
        plot(dates,sum(HD_loop,dims=2),color = "k",linewidth=1.5);
        axesData.set_xlim(dates[2],dates[end])
        axesData.set_title(vnamesPlot[ii])
        ylims = axis()[3:4];
        for iiii=1:1:size(NBER_recessions,1)
            axesData.fill([NBER_recessions[iiii,1],NBER_recessions[iiii,1],NBER_recessions[iiii,2],NBER_recessions[iiii,2]],[ylims[1],ylims[2],ylims[2],ylims[1]],facecolor="grey",alpha=0.45)
        end
        for ij in idx_HD[1,:]
            xpos[ij] = axesData.bar(dates,Int.(HD_loop[:,ij].>0).*HD_loop[:,ij],sizbar,facecolor = options_set.plot["color_map"][idx_HD[1,:],:][ij,:],
            bottom = (ij > 1 ? dropdims(sum(HD_loop[:,1:ij-1].*Int.(HD_loop[:,1:ij-1].>0),dims=2),dims=2) : 0), hatch = options_set.plot["hatch_map"][idx_HD[1,:],:][ij])
            xpos[ij] = axesData.bar(dates,Int.(HD_loop[:,ij].<0).*HD_loop[:,ij],sizbar,facecolor = options_set.plot["color_map"][idx_HD[1,:],:][ij,:],
            bottom = (ij > 1 ? dropdims(sum(HD_loop[:,1:ij-1].*Int.(HD_loop[:,1:ij-1].<0),dims=2),dims=2) : 0), hatch = options_set.plot["hatch_map"][idx_HD[1,:],:][ij])
        end         
      
        axesData.set_ylim(ylims[1],ylims[end])
        axesData.set_xlim(dates[3],dates[end])

        if in(ii,idx_filter)
            options_set.plot["annnualize"] ? axesData.set_ylabel("Annual Growth Rate (%)") : axesData.set_ylabel("Quarterly Growth Rate (%)")
        end
        if options_set.plot["disp_figure"]; 
            display(fig)
        end
        if options_set.plot["save_figure"]; 
            savefig(join([filename,vnamesPlot[ii],".pdf"]),format = "pdf",dpi = options_set.plot["precision"]); 
        end
     end

end


function HSD_plot_OneShock(HDM :: Array, TS :: TimeArray, options_set :: Options)
    TS = TS[1:133]
    HDM = HDM[1:133,:,:]
 
    nvar   = nshocks = options_set.options_var["nvar"]
    lags   =           options_set.options_var["lags"]
    HD     = Array{Float64}(undef,size(HDM,1),size(HDM,2),size(HDM,3))
    idx_HD = Array{Int64,2}(undef,nvar,nvar)
    
    for ii in 1 : nvar
        idx_HD[ii,:] = sortperm(dropdims(sum(abs.(HDM[lags+1:end,:,ii]),dims=1),dims=1),rev=true);
        idx_HD[ii,:] = [idx_HD[ii,findall(idx_HD[ii,:].!=6)] ; 6];
        idx_HD[ii,:] = [5;idx_HD[ii,findall(idx_HD[ii,:].!=5)]];
        HD[:,:,ii]   = HDM[:,idx_HD[1,:],ii];
    end

    HD              = [HD[:,1:1,:] sum(HD[:,2:end,:],dims=2)]
    HD[:,:,[1 5 6]] = HD[:,:,[1 5 6]]*100

    # Replace NaN with zero
    replace_nan(v) = map(x -> isnan(x) ? zero(x) : x, v)
    HD = replace_nan(HD);

    #  Pick Series to filter
    idx_filter =  indexin(options_set.plot["filter"], options_set.options_var["series"])

    # Annual Growth
    if options_set.plot["annnualize"]
        fun2 = var  ->  var[4:end,:] + var[3:(end-1),:] + var[2:(end-2),:] + var[1:(end-3),:]
    else
        fun2 = var  ->  var
    end

    (nsteps,)  = size(HD);
    filename   = options_set.paths["save_path_FIG"]*"HSD_";
    color_map  = [purple;[0.5 0.5 0.5]]
    snamesPlot = ["Financial" ; "Non-Financial"]
    sizbar 	   = 50
    # snamesPlot = options_set.options_sr["shocks"]
    vnamesPlot = colnames(TS)
    xpos 	   = Array{Any,1}(undef,size(HD,3));
    xneg 	   = Array{Any,1}(undef,size(HD,3));

    for ii in 1 : nvar
        if in(ii,idx_filter)
            HD_loop = fun2(diff(HD[:,:,ii],dims=1)); dates = timestamp(TS)[5:end];
        else
            HD_loop = HD[:,:,ii]; dates = timestamp(TS)[1:end];
        end

        fig,axesData  = subplots();
        for ij in 1:size(HD,2)
            xpos[ij] = axesData.bar(dates,Int.(HD_loop[:,ij].>0).*HD_loop[:,ij],sizbar,facecolor = color_map[ij,:],
            bottom = (ij > 1 ? dropdims(sum(HD_loop[:,1:ij-1].*Int.(HD_loop[:,1:ij-1].>0),dims=2),dims=2) : 0), hatch = options_set.plot["hatch_map"][ij])
            xpos[ij] = axesData.bar(dates,Int.(HD_loop[:,ij].<0).*HD_loop[:,ij],sizbar,facecolor = color_map[ij,:],
            bottom = (ij > 1 ? dropdims(sum(HD_loop[:,1:ij-1].*Int.(HD_loop[:,1:ij-1].<0),dims=2),dims=2) : 0), hatch = options_set.plot["hatch_map"][ij])
        end        
        if ii == 1
            axesData.legend(xpos[1:2],snamesPlot,labelspacing=0.075,loc="lower left")
        end
        plot(dates,sum(HD_loop,dims=2),color = "k",linewidth=1.5);
        axesData.set_xlim(dates[2],dates[end])
        axesData.set_title(vnamesPlot[ii])
        ylims = axis()[3:4];
        for iiii=1:1:size(NBER_recessions,1)
            axesData.fill([NBER_recessions[iiii,1],NBER_recessions[iiii,1],NBER_recessions[iiii,2],NBER_recessions[iiii,2]],[ylims[1],ylims[2],ylims[2],ylims[1]],facecolor="grey",alpha=0.45)
        end
        for ij in  1:size(HD,2)
            xpos[ij] = axesData.bar(dates,Int.(HD_loop[:,ij].>0).*HD_loop[:,ij],sizbar,facecolor = color_map[ij,:],
            bottom = (ij > 1 ? dropdims(sum(HD_loop[:,1:ij-1].*Int.(HD_loop[:,1:ij-1].>0),dims=2),dims=2) : 0), hatch = options_set.plot["hatch_map"][ij])
            xpos[ij] = axesData.bar(dates,Int.(HD_loop[:,ij].<0).*HD_loop[:,ij],sizbar,facecolor = color_map[ij,:],
            bottom = (ij > 1 ? dropdims(sum(HD_loop[:,1:ij-1].*Int.(HD_loop[:,1:ij-1].<0),dims=2),dims=2) : 0), hatch = options_set.plot["hatch_map"][ij])
        end         
      
        axesData.set_ylim(ylims[1],ylims[end])
        axesData.set_xlim(dates[3],dates[end])

        if in(ii,idx_filter)
            # options_set.plot["annnualize"] ? axesData.set_ylabel("Annual Growth Rate (%)") : axesData.set_ylabel("Quarterly Growth Rate (%)")
        end
        if options_set.plot["disp_figure"]; 
            display(fig)
        end
        if options_set.plot["save_figure"]; 
            savefig(join([filename,vnamesPlot[ii],"_Sum",".pdf"]),format = "pdf",dpi = options_set.plot["precision"]); 
        end
     end

end

function IRF_plot(IRF,options_set)

    series      = options_set.options_var["series"]
    horizon     = options_set.options_sr["horizon"]
    shocks      = options_set.options_sr["shocks"]
    resp_low    = IRF["resp_low"]
    resp_high   = IRF["resp_high"]
    resp_50     = IRF["resp_50"]
    resp_mean   = IRF["resp_mean"]

    o               = zeros(1,horizon);        # zero line
    tt              = collect(0:horizon-1);
    patch_time      = [tt ; reverse(tt,dims=1)];   # runs through time and back
    filename        = options_set.paths["save_path_FIG"]*"IRF_";

    for idx_shock  in 1 : nvar
        figIRF = figure(shocks[idx_shock],figsize=(8,8));# suptitle(F_names_imp(idx_shock,1:F_namesize_imp(idx_shock)))
        for idx_var in 1 : nvar

            rlow    =   reshape(resp_low[idx_var,idx_shock,:]  ,1,horizon) 	* 1;
            rhigh   =   reshape(resp_high[idx_var,idx_shock,:] ,1,horizon) 	* 1;
            r50     =   reshape(resp_50[idx_var,idx_shock,:]   ,1,horizon) 	* 1;
            rM      =   reshape(resp_mean[idx_var,idx_shock,:] ,1,horizon)  * 1;

            subplot(Int(ceil(sqrt(nvar))), Int(ceil(sqrt(nvar))), [idx_var])
            Hpatch  = [patch_time [rlow  reverse(rhigh,dims=2)]'];
            fill(Hpatch[:,1],Hpatch[:,2],"grey",alpha=0.85)
            # Median and zero line:
            if false
            plot(tt, rM[:],"--",  color = (0.6,0.,1.),  linewidth= 1.5);
            else
            plot(tt, r50[:],"--", color = (0.,0.5,0.9), linewidth= 1.5);
            end
            plot(tt, o[:],color = "k", linewidth = 0.75)
			ticklabel_format(axis="y", style="sci",scilimits = (0,0))
            xlim(tt[1],tt[end]); ylim(minimum(rlow).*1.3,maximum(rhigh).*1.3)
            grid(); 
            title(series[idx_var])
            tight_layout(); 
        end
        if options_set.plot["disp_figure"]; 
            display(figIRF)
        end
        if options_set.plot["save_figure"]; 
            savefig(join([filename,shocks[idx_shock],".pdf"]),format = "pdf",dpi = options_set.plot["precision"]); 
        end
    end

end


function SHK_plot(SHK,TS,options_set)

    SHK_low 	= Array{Float64,2}(undef,nvar,options_set.options_var["T"])
    SHK_high 	= Array{Float64,2}(undef,nvar,options_set.options_var["T"])
    SHK_50   	= Array{Float64,2}(undef,nvar,options_set.options_var["T"])
    SHK_mean	= Array{Float64,2}(undef,nvar,options_set.options_var["T"])
    
    for Shock in 1 : nvar
            for Time in 1 : options_set.options_var["T"]-lags
                SHK_low[Shock,Time]    =   percentile(SHK[Shock,:,:][Time,:],options_set.options_sr["low_band"])
                SHK_high[Shock,Time]   =   percentile(SHK[Shock,:,:][Time,:],options_set.options_sr["high_band"])
                SHK_50[Shock,Time]     =   percentile(SHK[Shock,:,:][Time,:],50);
            end # for Time
    end # for Shock
    
    shocks       = options_set.options_sr["shocks"]

    o               = zeros(1,options_set.options_var["T"]);        # zero line
    tt              = timestamp(TS)#collect(0:options_set.options_var["T"]-1);
    patch_time      = [tt ; reverse(tt,dims=1)];   # runs through time and back
    filename        = options_set.paths["save_path_FIG"]*"MedianShocks";

    fig = figure(figsize=(8,8));# suptitle(F_names_imp(idx_shock,1:F_namesize_imp(idx_shock)))
        for idx_var in 1 : nvar

            rlow    =   SHK_low[idx_var,:];
            rhigh   =   SHK_high[idx_var,:];
            r50     =   SHK_50[idx_var,:];

            subplot(Int(ceil(sqrt(nvar))), Int(ceil(sqrt(nvar))), [idx_var])
            Hpatch  = [patch_time [rlow ; reverse(rhigh,dims=1)]];
            fill(Hpatch[:,1],Hpatch[:,2],"grey",alpha=0.85)
            # Median and zero line:
            if false
            plot(tt, rM[:],"--",  color = (0.6,0.,1.),  linewidth= 1.5);
            else
            plot(tt, r50[:],"--", color = (0.,0.5,0.9), linewidth= 1.5);
            end
            plot(tt, o[:],color = "k", linewidth = 0.75)
			xlim(tt[1],tt[end])
			ticklabel_format(axis="y", style="sci",scilimits = (0,0))
            ylim(axis()[3:4].*1.3)
            grid(); 
            title(shocks[idx_var])
            tight_layout(); 
        end

        if options_set.plot["disp_figure"]; 
            display(fig)
        end
        if options_set.plot["save_figure"]; 
            savefig(join([filename,".pdf"]),format = "pdf",dpi = options_set.plot["precision"]); 
        end

end


function VD_plot(VD,idx_HD,options_set)

    series      = options_set.options_var["series"]
    horizon     = options_set.options_sr["horizon"]
    shocks      = options_set.options_sr["shocks"]
    
    nvar = nshocks = options_set.options_var["nvar"]
    # orderPlot   = [1 5 3 ; 2 4 6]
    orderPlot   = 1:nvar
    
    tt              =   collect(1:horizon-1);
    patch_time      =   [tt ; reverse(tt,dims=1)];   # runs through time and back
    o               =   zeros(1,horizon-1);    # zero line
    v50				=   Array{Float64,3}(undef,horizon-1,nvar,nvar)
    band			=   Array{Any,1}(undef,nvar)

    for idx_var in 1 : nvar
    	for idx_shock in 1 : nvar
            v50[:,idx_shock,idx_var] = reshape(VD[idx_var,idx_HD[1,idx_shock],:],1,horizon-1);
    	end
    end


    fig  = figure("Variance Decomposition",figsize=(12,8))
    for idx_var in 1 : nvar
    	subplot(Int(ceil(sqrt(nvar))), Int(ceil(sqrt(nvar))), orderPlot[idx_var])
        for idx_shock in 1 : nvar
    		if idx_shock == 1; vdist = 0; vpl = v50[:,1,idx_var]; else vdist = dropdims(sum(v50[:,1:idx_shock-1,idx_var],dims=2),dims=2);
    			global vpl = vpl + v50[:,idx_shock,idx_var];
            end
            
        # Plot stacked variance decompositions
    	global band[idx_shock] = fill_between(collect(0:size(v50)[1]-1),vpl,y2 = vdist,facecolor = options_set.plot["color_map"][idx_HD[1,:],:][idx_shock,:], hatch =  options_set.plot["hatch_map"][idx_HD[1,:],:][idx_shock])
        end
    	xlim(0,horizon-2)
    	ylim(0,1)
        title(series[idx_var],fontsize=27)
    	xticks(fontsize=18)
    	yticks(fontsize=18)
    end
    # idxx = [5; 1; 2; 3; 4; 6]
    idxx = idx_HD[1,:]# [5; 1; 2; 3; 4; 6]
    legend(labels=shocks[idxx],handles=[band[1:nvar]... ],fontsize=18)
    tight_layout()

    if options_set.plot["disp_figure"]; 
        display(fig)
    end
    if options_set.plot["save_figure"]; 
         savefig(join([options_set.paths["save_path_FIG"],"VD_Median.pdf"]),format = "pdf",dpi = options_set.plot["precision"]); 
    end

end

function IRF_modal_plot(Response,options_set)

    series      = options_set.options_var["series"]
    horizon     = options_set.options_sr["horizon"]
    shocks      = options_set.options_sr["shocks"]
    alphaHPD    = options_set.plot["alphaHPD"]

    # Modal Model
    logPostVec  = (Response.logPost);
    Ids         = sortperm(logPostVec);
    logPostVecS = logPostVec[Ids];
    idm         = Ids[end];
    modalHD     = Response.HDstore[:,:,:,idm];
    modalResp   = Response.f[:,:,:,idm];
    respBand    = Response.f[:,:,:,Ids[end-Int(ceil((1-alphaHPD)*size(logPostVec)[1])):end]];

    # Impulse Response Function
    tt          = 1:horizon;
    o           = zeros(horizon,1);

for idx_shock in 1 : nvar
    fig1 = figure(shocks[idx_shock],figsize=(8,8)); # suptitle(F_names_imp(idx_shock,1:F_namesize_imp(idx_shock)))
	subplots()
        for idx_var in 1 : nvar
            rMM          = dropdims(reshape(modalResp[idx_var,idx_shock,:],  (1, horizon)), dims = 1);
            respBandPlot = dropdims(reshape(respBand[idx_var,idx_shock,:,:], (1, horizon, 1 + Int(ceil((1-alphaHPD)*size(logPostVec)[1])))),dims=1);
            subplot(Int(ceil(sqrt(nvar))), Int(ceil(sqrt(nvar))), idx_var)
            # Median and zero line:
            for i in 1:1:size(respBand,4)
            plot(tt,respBandPlot[:,i],"r-");
            end
            plot(tt, rMM, color=(0,0,0), linewidth= 1.5);
            plot(tt, o,   color=(0,0,0), linewidth= 1)

            title(series[idx_var])
        end

        if options_set.plot["disp_figure"]; 
            display(gcf())
        end
        if options_set.plot["save_figure"]; 
            savefig(join([options_set.paths["save_path_FIG"],"IRF_MM_",shocks[idx_shock],".pdf"]),format = "pdf", dpi = options_set.plot["precision"]); 
        end
end
return
end
