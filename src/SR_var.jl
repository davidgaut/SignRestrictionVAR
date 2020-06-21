@doc doc"""
Takes draws from the Bvar results, retain IRFs satisfying the sign restrictions (files are saved progressively)

Inputs:
bvar_results :: Dict
options_set  :: Options
TS           :: TimeArray

Outputs:
SR_output :: Dict
"""

function SR_var(bvar_results :: Dict, var_results :: Dict, options_set :: Options,TS :: TimeArray)

options_sr  = options_set.options_sr
options_var = options_set.options_var

unpack(options_set.paths)
unpack(options_sr)
unpack(options_var)

ratio_alpha_s       = SharedArray{Int64}(length(DrawLeft))
ratio_alpha_Tot_s   = SharedArray{Int64}(length(DrawLeft))
Eps_s               = SharedArray{Float64}(length(DrawLeft),nvar,T-lags);
impres_alpha_s   	= SharedArray{Float64}(length(DrawLeft),nvar,nvar,horizon,1)
HD_save_Shock_s 	= SharedArray{Float64}(length(DrawLeft),T,nvar,nvar)
HD_save_Const_s 	= SharedArray{Float64}(length(DrawLeft),T,nvar)
HD_save_Init_s 		= SharedArray{Float64}(length(DrawLeft),T,nvar)
logPost_s 			= SharedArray{Float64}(length(DrawLeft))

start_time = Dates.unix2datetime(time())
println("Starting draws, number of procs is "*string(nprocs()))

TotalDraws  = length(readdir(join(save_path_IRF)));
patternFile = r"Draw_(.+)_.+.jld";
fileirf 	= readdir(join(save_path_IRF));
n1 			= map(x->match(patternFile,x),fileirf)
filter!(x -> nothing!= x, n1)
n1 			= [(n1[ij].captures[1]) for ij in 1 : size(n1)[1]]

if isempty(n1)
    n1 = 0
else
    n1 = maximum(parse.(Int,n1));
end

# randmatrix  = rand(MvNormal(zeros(nvar),Matrix{Float64}(I,nvar,nvar)),nvar);
# randvec     = rand(Normal(0,1), var_results["nvar"]*var_results["ncoeffs"]);
# SignRestrictionVAR.SR_draw(randmatrix,randvec,bvar_results,var_results,options_set,TS,1)


global draw_num = n1
while TotalDraws < draws
    
    randmatrix  = rand(MvNormal(zeros(nvar),Matrix{Float64}(I,nvar,nvar)),nvar);
    randvec     = rand(Normal(0,1), var_results["nvar"]*var_results["ncoeffs"]);
    
    # @showprogress 
    pmap(x -> SR_draw(randmatrix,randvec,bvar_results,var_results,options_set,TS,x), collect(draw_num: draw_num + nprocs() - 1));

    TotalDraws = length(readdir(join(save_path_IRF)));
    draw_num   = draw_num + nprocs()
end

finish = convert(Int, Dates.value(Dates.unix2datetime(time())- start_time))/1000/60/60;
println("Total elapsed time: ", finish, " hours. \n")
end


function SR_draw(randmatrix,randvec,bvar_results::Dict,var_results::Dict,options_set :: Options, TS :: TimeArray, w :: Int)
  
    unpack(options_set.options_var)
    unpack(options_set.paths)
    unpack(options_set.options_sr)

    
    horizon_max  = maximum(R[:,2,:]);
    impres_alpha = NaN*Array{Float64,4}(undef,nvar,nvar,horizon,alpha);
    HT_inverse   = NaN*Array{Float64,3}(undef,nvar,nvar,alpha);
    index_alpha  = zeros(Int64,alpha);

    # bvar_results = bvar_base(TS,options_set.options_var)
    PhiVAR, b, Q, ut, bfo, Fcomp, s2, Constant =
    bvar_draw(randmatrix,randvec,bvar_results["b"],var_results["xxu"],bvar_results["sxx"],bvar_results["sinv"],bvar_results["vres"],var_results["y"],var_results["X"],var_results["T"],var_results["lags"],var_results["nvar"],var_results["ncoeffs"],options_set.options_var["constant"]);
    
    A1           = Ref((cholesky(Q).U)');       # Chol decomposition of reduced shocks VAR-COVAR matrice
    A2           = PhiVAR;
    Constant     = Constant;
    RConditions  = reshape(R[:,1:2,:],nvar,2*nvar); 
    
    if all(Bool.((RConditions .== 1) + (RConditions .== 0))); 
        FastSign  = true; 
        Rsigns    = reshape(R[:,3,:],nvar,nvar); 
        idx_signs = findall(.!(Bool.(iszero.(Rsigns))))
        idx_sum   = sum(abs.(Rsigns),dims=1)
        RsignsV   = Rsigns[idx_signs]
        idx_nzero = findall(.!(cumprod(iszero.(Rsigns),dims=1)[end,:]))
    end
    

    # function rotate()

    ratio_alpha_Tot = ratio_alpha = maxdraw = success_count = 0
    # global success, bqr
        # for ix in 1:alpha
            while maxdraw <= drawcount
                ix = 1
                # global ratio_alpha_Tot, maxdraw
                    maxdraw += 1
                    success = false
                    (Qr,)   = qr_frs(rand(MvNormal(zeros(nvar), Matrix{Float64}(I,nvar,nvar)),nvar)); #check qr vs qr_frs
                    bqr     = similar(A1[])
                    mul!(bqr,(A1)[],Matrix(Qr));
                
                if FastSign & !(horizon_max == 1)
                    (horizon_max == 1) ? (IrfD = bqr) : (IrfD = impulsdtrf(A2,bqr,horizon_max));
                    SRR     = sign.(IrfD[:][idx_signs]).== RsignsV
                elseif FastSign & (horizon_max == 1)
                    crit = (Int.(sum(sign.(Rsigns./bqr),dims=1)) ./ idx_sum)[idx_nzero]
                    if any(1.0 .> crit .> -1.0)
                        SRR = 0
                    else
                        bqr[:,idx_nzero] = bqr[:,idx_nzero] .* crit'
                        SRR   = 1
                    end  
                else
                    (horizon_max == 1) ? (IrfD = bqr) : (IrfD = impulsdtrf(A2,bqr,horizon_max));
                    SRR     = BL_Condition(IrfD,R);
                end
            
            if SRR  == 1
                impres_alpha_check = impulsdtrf(A2,bqr,horizon);
                ratio_alpha        = ratio_alpha + 1;
                success            = true;
            elseif SRR  == -1
                bqr                = - bqr;
                impres_alpha_check = impulsdtrf(A2,bqr,horizon);
                ratio_alpha        = ratio_alpha + 1;
                success            = true;
            end
            
         ratio_alpha_Tot = ratio_alpha_Tot + 1;
           

if success

    success_count += 1
    impres_alpha[:,:,:,ix] = impres_alpha_check;

    # Computation of the Logposterior
vecB        = vec([Constant' ; reshape(vec(permutedims(A2,(2,3,1))),lags*nvar,nvar)]);
(Y, X)      = VARmakexy(values(TS),lags,constant);
Ft          = (X'*X)\(X'*Y);
F           = Ft';
Residuals   = Y - X*Ft;
FcompU      = [F[:,1+constant:nvar*lags+constant]; Matrix{Float64}(I,nvar*(lags-1),nvar*(lags-1)) zeros(nvar*(lags-1),nvar)];
sigma       = (1 / T)*(Y-X*Ft)'*(Y-X*Ft); # adjusted for # of estimated coeff per equation
EvecB       = vec(Ft);

logPost     = irfpdf(bqr,EvecB,Symmetric(X'*X),T,sigma,vecB);


# Computation of the Historical Variance Decomposition for accepted IRFs
#--------------------------------------------------------------------------

invA        = bqr;
Eps         = invA\transpose(Residuals);  # structural errors

# Contribution of each shock
invA_big           = zeros(nvar*lags,nvar);
invA_big[1:nvar,:] = invA;
Icomp              = [Matrix{Float64}(I,nvar,nvar) zeros(nvar,(lags-1)*nvar)];
HDshock_big        = zeros(lags*nvar,T,nvar);
HDshock            = zeros(nvar,T,nvar);

    for j in 1:nvar; # for each variable
        global eps_big         = zeros(nvar,T); # matrix of shocks conformable with companion
        eps_big[j,lags+1:end]  = Eps[j,:];
        for i in lags+1:T
            HDshock_big[:,i,j] = invA_big * eps_big[:,i] + FcompU * HDshock_big[:,i-1,j];
            HDshock[:,i,j]     = Icomp    * HDshock_big[:,i,j];
        end
    end


# Initial value
    HDinit_big      = zeros(lags*nvar,T);
    HDinit          = zeros(nvar,T);
    HDinit_big[:,1] = X[1,1+constant:nvar*lags+constant]';
    HDinit[:,1]     = Icomp * HDinit_big[:,1];
    for i in 2:T
        HDinit_big[:,i] = FcompU * HDinit_big[:,i-1];
        HDinit[:,i]     = Icomp  * HDinit_big[:,i];
    end

# Constant
    HDconst_big = zeros(lags*nvar,T);
    HDconst     = zeros(nvar, T);
    CC          = zeros(lags*nvar,1);
    if constant > 0
        CC[1:nvar,:] = F[:,1];
        for i in 2:T
            HDconst_big[:,i] = CC + FcompU * HDconst_big[:,i-1];
            HDconst[:,i]     = Icomp * HDconst_big[:,i];
        end
    end

# Reshape all HDs
HD_save_Shock  = zeros(T,nvar,nvar);  # [Nobs x shock x var]
    for i in 1 : nvar
        for j in 1 : nvar
            # println( [NaN*Array{Float64,1}(undef,lags); HDshock[i,1:end,j]])
            # println([HDshock[i,1:end,j]][1])
            HD_save_Shock[:,j,i]  = [NaN*Array{Float64,1}(undef,1); HDshock[i,2:end,j]];
            # HD_save_Shock[:,j,i]  =  [HDshock[i,:,j]][1];
        end
    end

HD_save_Const = HDconst[:,1:end]'[:];
HD_save_Init  = HDinit[:,1:end]'[:];

SR_output = Dict(
"ratio_alpha"       => ratio_alpha     ,
"ratio_alpha_Tot"   => ratio_alpha_Tot ,
"Eps"               => Eps             ,
"IRF"               => impres_alpha    ,
"A1"                => A1   ,
"A2"                => A2   ,
"HD_save_Shock"     => HD_save_Shock   ,
"HD_save_Const"     => HD_save_Const   ,
"HD_save_Init"      => HD_save_Init    ,
"logPost"           => logPost);

savefile = join(save_path_IRF)*"Draw_$(w)_$(success_count).jld";
save(savefile,"Draw",SR_output);

# TotalDraws = length(readdir(join(save_path_IRF)));

            
        end
    end
    return success_count./ratio_alpha_Tot
end

# return SR_output
# end

function impulsdtrf(B::Array{S,3},smat::Matrix{S},nstep::Int) where {S<:AbstractFloat}
    # function Response=impulsdtrf(B,smat,nstep)
    # Assumes the same model as in rfvar, except here only the By part is used.
    # smat is a square matrix of initial shock vectors.  To produce "orthogonalized
    # impulse Responses" it should have the property that smat'*smat=sigma, where sigma
    # is the Var(u(t)) matrix and u(t) is the residual vector.  One way to get such a smat
    # is to set smat=chol(sigma).  To get the smat corresponding to a different ordering,
    # use smat=chol(P*Sigma*P')*P, where P is a permutation matrix.
    # B is a neq x nvar x nlags matrix.  neq=nvar, of course, but the first index runs over
    # equations.  In Response, the first index runs over variables, the second over
    # shocks (in effect, equations), the third over time.
    # Code written by Christopher Sims.  This version 6/15/03.
    (neq,nvar,lags) = size(B);
    Response        = zeros(neq,nvar,nstep);
    Response[:,:,1] = smat;
    for it=2:nstep
       for ilag in 1:min(lags,it-1)
          Response[:,:,it] = Response[:,:,it] + B[:,:,ilag]*Response[:,:,it-ilag];
       end
    end
    return Response
end



function BL_Condition(candd::Array{S},R::Array{Int,3}) where {S<:AbstractFloat}


    nvar     = size(R,1);
    nshocks  = size(R,3);

    checkall      = (zeros(nvar,nshocks));
    checkall_flip = (zeros(nvar,nshocks));

   for ss in 1:nshocks
        for ii in 1:nvar
            if R[ii,1,ss] != 0
                if R[ii,3,ss] == 1
                    check = candd[ii,ss,R[ii,1,ss]:R[ii,2,ss]] .> 0;
                    checkall[ii,ss] = minimum(check);
                    # Check flipped signs
                    check_flip = candd[ii,ss,R[ii,1,ss]:R[ii,2,ss]] .< 0;
                    checkall_flip[ii,ss] = minimum((check_flip));
                elseif R[ii,3,ss] == -1
                    check = candd[ii,ss,R[ii,1,ss]:R[ii,2,ss]] .< 0;
                    checkall[ii,ss] = minimum(check);
                    # Check flipped signs
                    check_flip = candd[ii,ss,R[ii,1,ss]:R[ii,2,ss]] .> 0;
                    checkall_flip[ii,ss] = minimum(check_flip);
                elseif R[ii,3,ss] == 2
                    check = sum(prod(candd[Int.(R[ii,3,ss].==2),ss,R[ii,1,ss]:R[ii,2,ss]]) .> 0) .== R[ii,2,ss];
                    checkall[ii,ss] = minimum(check);
#                     if check && (ss==4)
#                         figure; plot(squeeze(candd((R(:,3,ss)==2),ss,R(ii,1,ss):R(ii,2,ss)))');
#                     end
                    # return

                end
                else
                    checkall[ii,ss] = 1;
                    # Check flipped signs
                    checkall_flip[ii,ss] = 1;
            end
        end
   end

   if minimum(minimum(checkall)) == 1
       logical = 1;
   elseif minimum(minimum(checkall_flip)) == 1
       logical = - 1;
   else
       logical = 0;
   end
   return logical
end


function qr_norm(WW :: Array{Float64,2})
x = qr!(WW); x = sign.(diag(x.R))'.*x.Q
return x
end

# QR factorization with non-negative diagonals in R
function qr_frs(WW :: Array{Float64,2})
M  = qr!(WW);
Qr = M.Q[:,:];
Rr = M.R;
for rx in 1:size(WW,1)
    if Rr[rx,rx] < 0
        Rr[rx,:]  = -Rr[rx,:];
        Qr[:,rx]  = -Qr[:,rx];
    end
end
return Qr, Rr
end

# Indepedent  and dependent variables preparation
# function var_companion(TS)
# (Y, X) = VARmakexy(values(TS),lags,constant);
# Nobs   = size(Y,1);
# d1     = data;
#
# ncoeff        = nvar * lags;
# nvarXeq       = nvar * lags;
# ntotcoeff     = ncoeff + constant;
# Nobse         = size(d1,1) - lags;
#
# # Compute companion form
# Ft          = (X'*X)\(X'*Y);
# SIGMA       = (1/(Nobse))*(Y-X*Ft)'*(Y-X*Ft); # adjusted for # of estimated coeff per equation
# sigma       = SIGMA;
# Residuals   = Y - X*Ft;
# EvecB       = vec(Ft);
# F           = Ft';
# FcompU      = [F[:,1+constant:nvar*lags+constant]; Matrix{Float64}(I,nvar*(lags-1),nvar*(lags-1)) zeros(nvar*(lags-1),nvar)];
# maxEig      = maximum(abs.(real(eigvals(FcompU))));
# yy          = Array{Float64,2}(undef,lags * nvar,Nobse-lags)
# end

function  VARmakexy(DATA,lags,constant)
# =======================================================================
# Builds the VAR process from the data-matrix DATA. It orders the data into
# the Y and X matrix --> Example: [x y] = [x(-1) y(-1) x(-2) y(-2)]
# =======================================================================
# [Y, X] = VARmakexy(DATA, lags, constant)
# -----------------------------------------------------------------------
# INPUT
#   DATA: matrix containing the original data
#   lags: lag order of the VAR
#   constant : 0, no constantant, no trend
#           1, constant, no trend
#           2, constant, trend
#           3, constant, trend, trend^2
# -----------------------------------------------------------------------
# OUTPUT
#   Y: dependent variable
#   X: independent variable
# =======================================================================



# Get dimesion of DATA
(nobs,) = size(DATA);

# Y matrix
Y = DATA[lags+1:end,:];

# X-matrix
if constant==0
    X = Array{Float64,2}(undef,nobs-lags,0);
    for jj in 0:lags-1
        X = hcat(DATA[jj+1:nobs-lags+jj,:],X);
    end

elseif constant==1 #constant
        X = Array{Float64,2}(undef,nobs-lags,0);
        for jj in 0:lags-1
            X = hcat(DATA[jj+1:nobs-lags+jj,:],X);
        end
       X = hcat(ones(nobs-lags,1), X);

elseif constant==2 # time trend and constant
    X = Array{Float64,2}(undef,nobs-lags,0);
    for jj in 0:lags-1
        X = hcat(DATA[jj+1:nobs-lags+jj,:],X);
    end
        trend=1:size(X,1);
        X = hcat(ones(nobs-lags,1),collect(trend), X);

elseif constant==3 # squared time trend, linear time trend, and constant
    X = Array{Float64,2}(undef,nobs-lags,0);
    for jj in 0:lags-1
        X = hcat(DATA[jj+1:nobs-lags+jj,:],X);
    end
        trend=1:size(X,1);
        X = hcat(ones(nobs-lags,1),collect(trend),collect(trend).^2, X);
end

return (Y, X)
end


function irfpdf(Atilde,EvecB,NT,nuT,ST,vecB)

# Purpose:
# This code computes the log of the joint posterior density of structural impulse responses up to constant
# input :
# Atilde: Obtained from rotating matrix A by matrix U, i.e., Atilde=A*U
# EvecB : Posterior mean of vec(B)
# NT    : N0+X'*X
# nuT   : degrees of freedom for the inverse-Wishart distribution
# ST    : Scale matrix for the inverse-Wishart distribution
# vecB  : Posterior draw of vec(B) = vec([c Phi1 ... Phip]')
# output:
# y     : Log of the joint posterior density up to constant
#
# Record of revisions:
# Date        Programmer            Description of change
# 02/15/2011  Atsushi Inoue         Original code
# 05/11/2011  Atsushi Inoue         Modified
# 02/26/2014  Jonas Arias' corrections incorporated
# 06/22/2018  Tom Doan's correction incorporated

n     = size(Atilde,1);           # The number of variables
p     = Int((size(vecB)[1]/n-1)/n);     # The order of the VAR model
index = kron(ones(n,1),[0;ones(n*p,1)]);


# Communication Matrix
Kn=Array{Float64}(undef,n*n,0);
for i in 1:n
    for j in 1:n
        E      = zeros(n,n);
        E[i,j] = 1;
        Kn     = [Kn reshape(E,n*n,1)];
    end
end

# Duplication Matrix
Dn=Array{Float64}(undef,n*n,0);
for j in 1:n
    for i in j:n
        E      = zeros(n,n);
        E[i,j] = 1;
        E[j,i] = 1;
        Dn     = [Dn reshape(E,n*n,1)];
    end
end
Dnplus = inv(Dn'*Dn)*Dn';

# Duplication Matrix" for Non-Symmetric Matrices (such as A)
Dbarn=Array{Float64}(undef,n*n,0);
for j in 1:n
    for i in j:n
        E      = zeros(n,n);
        E[i,j] = 1;
        Dbarn  = [Dbarn reshape(E,n*n,1)];
    end
end

# En Matrices
En=Array{Float64}(undef,n*n,0);
for j in 2:n
    for i in 1:j-1
        E      = zeros(n,n);
        E[i,j] =  1;
        E[j,i] = -1;
        En     = [En reshape(E,n*n,1)];
    end
end
Enplus = inv(En'*En)*En';

# Compute reduced-form impulse responses
B     = (reshape(vecB,1+n*p,n))';
Phi   = [B[:,1+1:end] ; ones(n*(p-1),n*(p-1)) zeros(n*(p-1),n)];
Theta = Array{Float64}(I,n,n);
for i  in  1:p
    Phii  = Phi^i;
    Theta = [Theta Phii[1:n,1:n]];
end

# Compute the Jacobian of vec(Theta)
Sigma       = Atilde*Atilde';
A           = cholesky(Sigma).U';
U           = inv(A)*Atilde;
v           = eigen(U).values;
j=1;
while j<=n
    j
     if all((abs.(real.(v[j])+1).<sqrt.(eps())).&(abs.(imag.(v[j]).<sqrt.(eps()))))
         W=diagm(0=>dropdims([-1; ones(n-1,1)],dims=2));
        A=A*W;
        U=W*U;
         #Atilde=(A*W)*(W*U)=A*U is unchanged
         j=n+1; # Terminate while statement
     else
         j=j+1;
     end;
end;
Idn         = Array{Float64}(I,n,n);
S           = Idn-2*inv(Idn+U);
s           = Enplus*reshape(S,n*n,1);
Ju          = kron(inv(Idn-S)',inv(Idn-S))*En;
J1          = [kron(Idn,A)*Ju kron(U',Idn)*Dbarn];
J3          = Dnplus*(kron(A,Idn)+kron(Idn,A)*Kn)*Dbarn;

idx         = LinearIndices(findall(index.==1));
logdetJ     = -log(abs(det(J1)))-0.5*n*p*log(det(Sigma))+log(abs(det(J3)));
VvecB       = kron(Sigma,inv(NT));
VvecBinv    = inv(VvecB[idx,idx]);
logdetVvecB = logdet(VvecB[idx,idx]);
f           = logdetJ-0.5*logdetVvecB-0.5*(vecB[idx,1]-EvecB[idx,1])'*VvecBinv*(vecB[idx,1]-EvecB[idx,1]);
f           = f-0.5*(nuT+n+1)*log(det(Sigma))-0.5*tr(nuT*ST*inv(Sigma));
f           = f-(n-1)*log(det(Idn+S));
end



function qr_fast(nvar)
    Q = Ref(qr(rand(MvNormal(zeros(nvar), Matrix{Float64}(I,nvar,nvar)),nvar)).Q[:,:])
    return Q
end

