@doc doc"""

Load previous successful draws
"""

function get_draws!(options_set :: Options)

unpack(options_set.options_sr)
unpack(options_set.paths)

# patternFile = r"Draw_(.+).jld";
patternFile = r"Draw_(.+)_(.+).jld";
fileirf 	= readdir(join(save_path_IRF));
n1 			= map(x->match(patternFile,x),fileirf)
filter!(x -> nothing!= x, n1)
n1 			= [(n1[ij].captures[1]) for ij in 1 : size(n1)[1]]
n1 		    = parse.(Int,n1);
DrawDone    = sort(n1);


# DrawLeft = findall(indexin(collect(1:draws),DrawDone) .== nothing);
DrawLeft = options_set.options_sr["draws"] - length(DrawDone)
options_set.options_sr["DrawDone"] = DrawDone
options_set.options_sr["DrawLeft"] = DrawLeft

println(string(length(DrawDone))*" draws completed in a previous session, " * string(DrawLeft) *" draws to go.")

return options_set
end
