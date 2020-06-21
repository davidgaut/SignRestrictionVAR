
using MAT, Distributions, DataFrames, Distributed
# addprocs() # need to open VSC as admin
# addprocs([2, 3, 4, 5])
@everywhere cd("D:\\Julia_File\\CodeVAR_Jl\\")
@everywhere using Pkg, MAT, Distributions, DataFrames, Distributed, ProgressMeter
@everywhere Pkg.activate("SignRestrictionVAR") # Make sure the right path to SignRestrictionVAR is provided to instantiate @everywhere
Pkg.instantiate()
@everywhere using SignRestrictionVAR



close("all")
using Dates, JLD

options_set = SignRestrictionVAR.Options()

FL = 1; FL1 = 0; FL2 = 1; FL3 = 0;
    # Supply     Demand    Invest   Monetary  Financial    Residual
RS = [  FL         FL         FL       FL        1           FL1       # Output
       -FL2        FL2        FL2      FL2       0           0         # Inflation
        0          FL2        FL2     -FL2      FL3          0         # IR
        0         -FL2        FL2       0        0           0         # Inv
        FL         FL         FL     1-FL3       1           FL1       # Loan
        FL         FL         FL     1-FL3      -1           FL1       # Bond
];



options_set.options_sr["RS"] = RS
options_set.options_var["date_sample"] = [Date(1985,01,01) Date(2019,12,01)];

options_set.plot["plot"]            = true 
options_set.paths["save_path"]      = "D:\\Julia_File\\CodeVAR_Jl\\"
options_set.paths["save_name"]      = "Baseline_NW_2019_1"
options_set.paths["data_path"]      = "D:\\Julia_File\\CodeVAR_Jl\\BLD_BaseData.xlsx"
# options_set.paths["data_path"]      = "D:\\Codes\\Matlab_File\\BLD_Codes\\BLD_data_19Q4_SR.xls"
options_set.paths["data_path"]      = "D:\\Codes\\Matlab_File\\BLD_Codes\\BLD_data_19Q4.xls"
# options_set.options_sr["shocks"]    = ["Supply","Demand","Investment","Monetary","Financial","Residual","Residual2"]
options_set.plot["filter"]          = ["Output" "Loans" "Bonds"]
options_set.plot["annnualize"]      = true
 
# options_set.paths["data_path"]      = "D:\\Codes\\Matlab_File\\BLD_Codes\\BLD_data_19Q4_Spread.xls"
# options_set.paths["save_name"]      = "Baseline_NW_2019_SpreadInvest_ShortSampleNoResid"
# options_set.options_sr["shocks"]    = ["Supply","Demand","Investment","Monetary","Financial","Residual","Residual2"]
options_set.options_sr["shocks"]    = ["Supply","Demand","Investment","Monetary","Financial","Residual"]
# FL = 1; FL1 = 0; FL2 = 0; FL3 = 0; FL4 = 1;
#     # Supply     Demand    Invest   Monetary  Financial    Residual
# RS = [  FL         FL         FL       FL        1      FL4   FL4       # Output
#        -FL2        FL2        FL2      FL2       0      0     0         # Inflation
#         0          FL2        FL2     -FL2      FL3     0     0         # IR
#         0          FL2        FL2     -FL2      FL3     0     0         # IR
#         0         -FL2        FL2       0        0      0     0         # Inv
#         FL         FL         FL     1-FL3       1      FL4   FL4     # Loan
#         FL         FL         FL     1-FL3      -1      FL4   FL4     # Bond
# ];

# # RS = [  FL         FL         FL       FL        1          0        # Output
# #        -FL2        FL2        FL2      FL2       0          0        # Inflation
# #         0          FL2        FL2     -FL2      FL3         0        # IR
# #         0         -FL2        FL2       0        0          0        # Inv
# #         FL         FL         FL     1-FL3       1          0        # Loan
#         FL         FL         FL     1-FL3      -1          0        # Bond
# ];
options_set.options_sr["RS"] = RS

options_set.options_sr["draws"] = 0000

options_set      = SignRestrictionVAR.initialize!(options_set)
TS, options_set  = SignRestrictionVAR.load_data!(options_set)
options_set      = SignRestrictionVAR.lag_crit!(TS,options_set)

var_results      = SignRestrictionVAR.var_OLS(TS,options_set.options_var)
bvar_results     = SignRestrictionVAR.bvar_base(var_results,options_set.options_var)

options_set      = SignRestrictionVAR.get_draws!(options_set)
SRvar_results    = SignRestrictionVAR.SR_var(bvar_results,var_results,options_set,TS) 
Results          = SignRestrictionVAR.load_results(options_set)

if false
A1 =  cholesky(dropdims(mean(bvar_results["b"][:,:,:],dims=3),dims=3)).U'
BetaCoeff = mean(bvar_results["PhiVAR_a"],dims=3)[2:end,:]

vmat = dropdims(mean(Response.f[:,:,1,:],dims=3),dims=3) # A1 * Rotation
vmat = dropdims(median(Response.f[:,:,1,:],dims=3),dims=3) # A1 * Rotation

nvars       = 6
nlags       = nlags
filter_data = true
frequency   = [6,32]


VD_bc   = VAR_VD(BetaCoeff,vmat,nvars,nlags,filter_data,frequency)
TableDS = convert(DataFrame,[vnames VD_bc])
DataFrames.rename!(TableDS,Symbol.(["Variables" ; snames]))
Results          = SignRestrictionVAR.VAR_VD(options_set) 
end

SignRestrictionVAR.do_figures(Results,TS,options_set)


if false

var_results      = SignRestrictionVAR.var_OLS(TS,options_set.options_var)
bvar_results     = SignRestrictionVAR.bvar_base(var_results,options_set.options_var)
w = 1


randmatrix  = rand(MvNormal(zeros(nvar),Matrix{Float64}(I,nvar,nvar)),nvar);
randvec     = rand(Normal(0,1), var_results["nvar"]*var_results["ncoeffs"]);

Beta,b,Q,ut,bfo,Fcomp,s2,Constant  = bvar_draw(randmatrix,randvec,bvar_results["b"],var_results["xxu"],bvar_results["sxx"],
bvar_results["sinv"],bvar_results["vres"],var_results["y"],var_results["T"],var_results["lags"],var_results["nvar"],var_results["ncoeffs"],options_set.options_var["constant"]);


# rfvarconstant(values(TS),2,[0 0],1)[3] - var_results["xxu"]

Beta1, b1, Q1,nx, ut1, bfo1, Fcomp1, B1,s21, Constant1 = VARmin(values(TS),2,1,[0 0]);




Bnorm,Q = niw(b,sxx,sinv,T,randmatrix,randvec);




(y_0,q,constV,prior) = (values(TS),2,1,[0 0]);

# Data
y             = y_0[q+1:end,:];
(T,K)         = size(y);

# Computation of hyperparameters
# Posterior with mean and variance parameters corresponding to OLS estimates
(B, u, xxu,xxi, nx) = rfvarconstant(y_0,q,prior,constV);

# Normal–Inverse Wishart Priors
M           = u'*u;
Q           = rand(InverseWishart(T-K*q-1,M));
Q           = M/(T-K*q-1);
XX          = kron(Q,xxi);
# Bnorm       = rand(MvNormal(vec(B), XX));<
Bnorm       = vec(B)
Bnorm       = Bnorm';
b           = reshape(Bnorm,K*q + nx,K);
bb          = b[1:K*q,:];
bbb         = reshape(bb,K,q,K);        # variables, q, equations
bbb         = permutedims(bbb,[3,1,2]); # equations, variables, q
Constant    = b[K*q + 1,:];

# Calculate Residual and bfo for each draw
ut    = y - xxu*b; # Residuals
bfo   = xxu*b;     # Forcasted variables
F     = b';
s2    = ut'*ut./T;  # Estimated error variance

# Companion Matrix
Fcomp = [F[:,1:K*(q)+0];  Matrix{Float64}(I,K*(q-1),K*(q-1)) zeros(K*(q-1),K)];



randmatrix  = Matrix{Float64}(I,nvar,nvar);
randvec     = 0*rand(Normal(0,1), var_results["nvar"]*var_results["ncoeffs"]);

(randmatrix,randvec,b,xxu,sxx,sinv,vres,y,T,lags,nvar,ncoeffs,constant) = (randmatrix,randvec,bvar_results["b"],var_results["xxu"],bvar_results["sxx"],
bvar_results["sinv"],bvar_results["vres"],var_results["y"],var_results["T"],var_results["lags"],var_results["nvar"],var_results["ncoeffs"],options_set.options_var["constant"]);

# Minnesota Method
Bnorm,Q = niw(b,sxx,sinv,T,randmatrix,randvec);

# NIW
Q           = rand(InverseWishart(T-nvar*lags-1,(T-nvar*lags-1)*Matrix(bvar_results["vres"])));
xxi         = inv(X'*X)
XX          = kron(Q,xxi);
Bnorm       = rand(MvNormal(vec(B), XX));

b           = reshape(Bnorm',nvar*lags + constant,nvar);
bb          = b[1:nvar*lags,:];
bbb         = reshape(bb,nvar,lags,nvar);        # variables, lags, equations
bbb         = permutedims(bbb,[3,1,2]);          # equations, variables, lags

if constant
    Constant = b[nvar*lags + 1,:];
else
    Constant = Matrix(undef,nvar,1)
end

# Residual and bfo for each draw
ut    = y - xxu*b; # Residuals
bfo   = xxu*b;     # Forcasted variables
F     = b';
s2    = ut'*ut./T;  # Estimated error variance

# Companion Matrix
Fcomp = [F[:,1:nvar*(lags)];  Matrix{Float64}(I,nvar*(lags-1),nvar*(lags-1)) zeros(nvar*(lags-1),nvar)];

# Beta 
PhiVAR  = bbb;
end







