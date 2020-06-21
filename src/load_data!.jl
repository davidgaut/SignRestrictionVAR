@doc doc"""

Load data located in path_data and sheet in a TimeSeries object.

By default uses first Excel sheet

Returns TimeSeries, and options for VAR
"""

function load_data!(options_set :: Options, sheet = 1)

path_data = options_set.paths["data_path"]
data_loc  = join([options_set.paths["save_path_FIG"],"DataPlot.pdf"])

database   = readxlsheet(path_data,sheet);
Time       = map(Date,database[2:end,1])
Data       = (database[2:end,2:end])
Names      = (database[1,2:end])
TS		   = TimeArray(Time,float.(Data),Symbol.(Names))
date_VAR   = collect(options_set.options_var["date_sample"][1]:Month(3):options_set.options_var["date_sample"][2])
date_VAR   = intersect(Time,date_VAR);
TS         = TS[date_VAR]

if options_set.plot["plot"]
    fig, axesData = subplots(Int(round(size(TS,2)/2)),2,figsize = (6,6));
    for ii = 1 : size(TS,2)
    	ax = axesData[ii]
    	ax.plot(timestamp(TS),values(TS)[:,ii]);
    	# for iiii=1:1:size(NBER_recessions,1)
    	# ax[:fill]([NBER_recessions[iiii,1],NBER_recessions[iiii,1],NBER_recessions[iiii,2],NBER_recessions[iiii,2]],[ylims[1],ylims[2],ylims[2],ylims[1]],facecolor="grey",alpha=0.45)
    	# end

    	ax.set_title(colnames(TS)[ii])
    	ax.set_xticks(timestamp(TS)[1:30:end])
    	ax.set_xlim(timestamp(TS)[1],timestamp(TS)[end])
		dates_str = map(t->string.(t)[1:4],timestamp(TS)[1:30:end])
		ax.set_xticklabels(dates_str)
    	ax.grid(1)
    end
    if size(TS,2) < 6; delaxes(axesData[6]); end
    tight_layout()
    savefig(data_loc,format = "pdf",dpi = options_set.plot["precision"])

end
options_set.options_var["nvar"]   = (size(TS)[2])
options_set.options_var["series"] = string.(colnames(TS))

ColorMaps     = get_cmaps()
ColorMaps     = ColorMaps[39]
ColorMapPlot  = Array{Float64}(undef,options_set.options_var["nvar"],3)
for ii  in 1 : options_set.options_var["nvar"]
    global ColorMapPlot[ii,:] = collect(ColorMaps(ii)[1:3])
end
options_set.plot["color_map"] = ColorMapPlot

TableDS = convert(DataFrame,[options_set.options_sr["shocks"] options_set.options_var["series"]])
DataFrames.rename!(TableDS, :x1 => :Shocks)
DataFrames.rename!(TableDS, :x2 => :Variables)
println(TableDS)

println("---------------------------");
println("Starting Date: ", (timestamp(TS)[1]));
println("Ending Date:   ", (timestamp(TS)[end]));
println("---------------------------");


RS = options_set.options_sr["RS"]

	  # Impact Restrictions
	  R1 = zeros(size(RS,1),3,size(RS,2));

	  # Set durations
				  for ij in 1 : size(RS,2)
					  R1[:,3,ij] = RS[:,ij];
					  R1[findall(RS[:,ij].!=0),1:2,ij] = ones(sum(.!(RS[:,ij] .== 0)),2);

					  # Restrictions for loans and bonds
			#         if (ij == 2)
			# #            R1(5:6,1,ij) = 1 .* ones(2,1);
			# #            R1(5:6,2,ij) = 4 .* ones(2,1);
			#            R1[4:5,1,ij] = 1 .* ones(2,1);
			#            R1[4:5,2,ij] = 1 .* ones(2,1);
			#         end

					# Restrictions for loans and bonds
			#         if (ij == 2) || (ij == 4)
			# #            R1((RS(:,ij) == 2),1,ij) = 1 .* ones(sum((RS(:,ij) == 2)),1);
			# #            R1((RS(:,ij) == 2),2,ij) = 4 .* ones(sum((RS(:,ij) == 2)),1);
			#
			#            R1(5:6,1,ij) = 1 .* ones(2,1);
			#            R1(5:6,2,ij) = 4 .* ones(2,1);
			#         end
				  end

			# Final restriction matrix
			options_set.options_sr["R"] = Int.(R1);

			# options_set.options_sr["R"] = string.(collect(1:size(R)[1]))
			TSR = convert(DataFrame, [options_set.options_var["series"] RS]);
			DataFrames.rename!(TSR,(Symbol.(vcat("Variables", options_set.options_sr["shocks"]))));
			println(TSR)


return TS :: TimeArray, options_set :: Options

end
   