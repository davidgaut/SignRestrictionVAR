

@doc doc"""
Create folders for output from the VAR, retained IRFs and generated graphs

Inputs:
options
save_path: location of the folder
save_name: name     of the folder

Outputs:
return Options with paths
"""

function initialize!(options_set :: Options)

    
     VAR_path 	= "VAR\\";
     IRF_path 	= "DRAWS\\";
     FIG_path	= "FIG\\";

     all_path      = [VAR_path;IRF_path;FIG_path];
     save_path     = options_set.paths["save_path"]*options_set.paths["save_name"]*"\\";
     save_path_VAR = [save_path,VAR_path];
     save_path_IRF = [save_path,IRF_path];
     save_path_FIG = [save_path,FIG_path];

    for ii = 1 : length(all_path)
    	mkpath(save_path*all_path[ii]);
    end

    println("Creating folder in " * save_path)

    options_set.paths["save_path_VAR"] = join(save_path_VAR)
    options_set.paths["save_path_IRF"] = join(save_path_IRF)
    options_set.paths["save_path_FIG"] = join(save_path_FIG)

                
    return options_set

end
