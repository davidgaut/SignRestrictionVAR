@doc doc"""

Compute Bayesian VAR model.
Takes arguments TS::TimeArray and options_set::Options and returns a dict containing results.
"""

# function bvar(TS :: TimeSeries.TimeArray,options_set :: Options)

#     unpack(options_set.options_var)

#     #--------------------------------------------------------------------------
#     b_a     	 = SharedArray{Float64}( lags * nvar + constant, nvar, IterationsNbr + 1)
#     Q_a     	 = SharedArray{Float64}( nvar , nvar, IterationsNbr + 1)
#     PhiVAR_a	 = SharedArray{Float64}( nvar, nvar, lags,  IterationsNbr + 1)
#     ut_a    	 = SharedArray{Float64}( T - lags, nvar, IterationsNbr + 1)
#     bfo_a   	 = SharedArray{Float64}( T - lags, nvar, IterationsNbr + 1)
#     Fcomp_a 	 = SharedArray{Float64}( lags * nvar , lags * nvar , IterationsNbr + 1)
#     B_a          = SharedArray{Float64}( lags * nvar + constant, nvar, IterationsNbr + 1)
#     Constant_a   = SharedArray{Float64}( 1, nvar, IterationsNbr + 1)


#     var_results  = var_OLS(TS,options_set.options_var)
#     bvar_results = bvar_base(var_results,options_set.options_var)
    
  
#     randmatrix  = rand(MvNormal(zeros(nvar),Matrix{Float64}(I,nvar,nvar)),nvar);
#     randvec     = rand(Normal(0,1), nvar*ncoeffs);
    
# # @sync @distributed
#  for w in 1 : IterationsNbr + 1
        
#     PhiVAR_a[:,:,:,w], b_a[:,:,w], Q_a[:,:,w], ut_a[:,:,w], bfo_a[:,:,w], Fcomp_a[:,:,w], s2, Constant_a[:,:,w] = bvar_draw(bvar_results["b"],var_results["xxu"],bvar_results["sxx"],bvar_results["sinv"],bvar_results["vres"],var_results["T"],var_results["lags"],var_results["nvar"],var_results["ncoeffs"],options_set.options_var["constant"])
#     # unpack(bvar_draw_it)                                                                                       
#     # Q_a[:,:,w]        = Q;
#     # PhiVAR_a[:,:,:,w] = PhiVAR;
#     # ut_a[:,:,w]       = ut;
#     # bfo_a[:,:,w]      = bfo;
#     # Fcomp_a[:,:,w]    = Fcomp;
#     # B_a[:,:,w]        = B;
#     # Constant_a[:,:,w] = Constant;
# end

# b_a        = b_a[:,:,b_sel];
# Q_a        = Q_a[:,:,b_sel];
# PhiVAR_a   = PhiVAR_a[:,:,:,b_sel];
# ut_a       = ut_a[:,:,b_sel];
# bfo_a      = bfo_a[:,:,b_sel];
# Fcomp_a    = Fcomp_a[:,:,b_sel];
# # B_a        = B_a[:,:,b_sel];
# Constant_a = Constant_a[:,:,b_sel];

# bvar_a = Dict("b_a"=>b_a,"Q_a"=>Q_a,"PhiVAR_a"=>PhiVAR_a,"ut_a"=>ut_a,"bfo_a"=>bfo_a,"Fcomp_a"=>Fcomp_a,"Constant_a"=>Constant_a,"T"=>T)

# return bvar_a
# end

function var_OLS(TS::TimeSeries.TimeArray,options_var::Dict)

    unpack(options_var)

    ydata         = values(TS)  
    (T,nvar)      = size(ydata);
    # var_results   = var_base(TS,options_var)
    
    y = ydata[lags+1:end,:];
    X = zeros(T - lags,nvar * lags);
    for lags_count in 1:lags
       X[:,1+(lags_count-1)*nvar:lags_count*nvar] = ydata[(lags+1-lags_count):(T-lags_count),:];
    end
    if constant; X = [X ones(T - lags,1)]; end
    
    # OLS VAR
    Beta = X \ y;
    Res  = y - X * Beta;     
    xxu  = X;    
    vres = (1 / (T - lags)) * (Res' * Res);
    stds = sqrt.(diag(vres));
    
    ncoeffs   = nvar * lags + constant; 
    

    vres = Symmetric(vres)
    return var_results = Dict("Beta"=>Beta,"Res"=>Res,"vres"=>vres,"stds"=>stds,"y"=>y,"X"=>X,"xxu"=>xxu,"ncoeffs"=>ncoeffs,"T"=>T,"nvar"=>nvar,"lags"=>lags)

end



function bvar_base(var_results::Dict,options_var::Dict)

    unpack(var_results)

    # Options
    quarterly   = 1
    prior       = 1
    lev         = zeros(nvar,1)
    if prior == "Minnesota" # Minnesota Priors
        println("Minnesota Prior")
        hm, bm       = minneprc(y,X,lags,quarterly,constant,lev,prior);
    else                    # Diffuse Priors
        hm = zeros(ncoeffs*nvar,1);
        bm = zeros(ncoeffs*nvar,1);   
    end    
    
    xxx         = inv(Symmetric(kron(I(nvar),X'*X) + Matrix{Float64}(I, size(hm,1), size(hm,1)) .* hm));
    bb          = xxx * (vec(X'*y) + bm);      # stacked posterior means of regression coeffs: vec(bpost) = inv[inv(H) + kron(I,X'X)] * [vec[X'Y] + inv(H)*bprior]
    b           = reshape(bb,ncoeffs,nvar);    # posterior coefficient matrix 
    sxx         = cholesky(xxx).U';
    sinv        = cholesky(inv(vres)).U;
    
    return bvar_results = Dict("b"=>b,"sxx"=>sxx,"sinv"=>sinv,"vres"=>vres)

end


function bvar_draw(randmatrix,randvec,b,xxu,sxx,sinv,vres,y,X,T,lags,nvar,ncoeffs,constant)


    # Bnorm, Q = niw(b,sxx,sinv,T,randmatrix,randvec);

    Q           = rand(InverseWishart(T-nvar*lags-1,(T-nvar*lags-1)*Matrix(vres)));
    xxi         = inv(Symmetric(X'*X))
    XX          = kron(Q,xxi);
    Bnorm       = rand(MvNormal(vec(b), XX));

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
    PhiVAR  = bbb;# # [equations, variables, lags] to match impulsdt.m

    # Result
    # bvar_results = Dict("PhiVAR"=>PhiVAR,"b"=>b,"Q"=>Q,"nx"=>nx,"ut"=>ut,"bfo"=>bfo,"Fcomp"=>Fcomp,"B"=>B,"s2"=>s2,"Constant"=>Constant)
    # return bvar_results
    return PhiVAR, b, Q , ut, bfo, Fcomp, s2, Constant

end



function minneprc(y,x,nlag,quarterly,constant,lev,prior)

   t,nvar = size(y);      #number of variables
   k = constant+nvar*nlag; #number of regressors per equation   
   
   #values for tightness paramters
   lambda0 = prior;
   lambda1 = .2;       # = gamma in Hamilton pg 361
   lambda2 = .5;       # = w
   lambda3 = 2;        # modified from lambda 3 = 1 to match Angeletos(2020) 
   lambda4 = 1e5;      # prior on intercept
   
   function seqa(a,b,c)
      seq = collect(a:b:(a+b*(c-1)));
      return seq
   end

   #create monthly lag decay to match 1/q decay in quarterly data where q = quarters
   if Bool.(quarterly)
       ld = seqa(1.0,1.0,nlag).^(-lambda3);    #produces a regularly sequence of lag decay [note ^-mu[4]] */
   else
      j = ceil(nlag/3)^(-lambda3);   # last quarter [rounded up] eg. l2=13=>xx2=5
      b = 0;
      if nlag>1;
         b = ( log(1.0)-log(j) ) / (1-nlag);
       end
      a  = exp(-b);
      ld = a*exp(b*seqa(1.0,1.0,nlag));  # Tao's lag decay to match 13th lag 
   end
   
   #add other parameters and take the squared inverse
   #this is H^-1 for own lags
   ld = (lambda0*lambda1*lambda2*ld).^(-2);  #note: the weight lambda2 is taken
                                             #off from the own lag coeffs
                                             #later on (line 107)
   
   #Scale factors from univariate OLS AR's to compute stdevs of priors of lags
   #of other variables
   s = zeros(nvar,1);
   i = 1;
   
   #compute variance of residuals of univariate AR regressions
   for i in 1:nvar;
      #y = Y(nlag+1:rows(Y),i);
      #t = rows(y);
      #x = ones(t,1);
      xi = x[:,k];
      for j in 1:nlag; 
         xi=[xi x[:,i+(j-1)*nvar]];
          #x=[x Y(nlag+1-j:rows(Y)-j,i)];  
      end
      bsh = inv(xi'*xi)*xi'*y[:,i];
      u = y[:,i] - xi*bsh;
      s[i] = (u'*u)/t;
   end
   
   #prepare construction of inv(H) for each equation and add prior for intercept
      if lambda4 > 0
       if constant
           H = [kron(ld,s); ((lambda0*lambda4).^(-2))];
       else
           H = kron(ld,s);
       end
   
   elseif lambda4 == 0;
       if constant
           H = [kron(ld,s); 0];
       else
           H = kron(ld,s);
       end
   end
  
   #construct output vectors
   hm = zeros(k*nvar,1);
   bm = zeros(k*nvar,1);

      # stack the H and distinguish own lag vs other 
   for i in 1:nvar;
       hadd = copy(H);      #NOTE: H/s(i) MEANS WE NORMALIZE BY VARIANCE OF DEPENDENT VARIABLE 
       for j in 0:(nlag-1);
          hadd[i+(j*nvar)] = (lambda2^2)*H[i+(j*nvar)];
       end
       hm[k*(i-1)+1:k*i,1] = hadd;       #stacking for each equation
   
       badd = zeros(k,1);
       
       #XXXXXXXXXXXXXXXXXXXXXXXXXXXX#
       if lev[i]==1;   # if lev=1, put unit root prior on own first lag; otherwise prior = 0 
          bm[k*(i-1)+i,1]=hm[k*(i-1)+i]*1; #note: identical results than with Chris' formula
       end
    end

       return hm, bm
   end


   function niw(BETAHAT,SXX,SINV,DF,RANDMATRIX,RANDVEC)

    NVARS   = size(SINV,1);
    NPARAMS = size(BETAHAT,1);
    RANW    = RANDMATRIX/sqrt(DF); #randn(NVARS,DF)/sqrt(DF);
    RANTR   = RANW'*SINV;
    WISH    = inv(RANTR'*RANTR);
    SWISH   = cholesky(Symmetric(WISH)).U';
    RANC    = RANDVEC; #randn(NPARAMS*NVARS,1);
    
    V=zeros(NVARS*NPARAMS,NVARS*NPARAMS);
    for i in 1:NVARS;
        for j in 1:NVARS;
            V[1+(i-1)*NPARAMS:i*NPARAMS,1+(j-1)*NPARAMS:j*NPARAMS] = SWISH[i,j] * SXX[1+(i-1)*NPARAMS:i*NPARAMS,1+(i-1)*NPARAMS:i*NPARAMS];
        end
    end
            
    SHOCK = V * RANC;
    BETADRAW = vec(BETAHAT) + SHOCK;

    return BETADRAW,WISH
end






# function bvar_base(TS::TimeSeries.TimeArray,options_var::Dict)

#     unpack(var_results)

#     # if false
#     # d         = svd(concatd2(X),full=false);
#     # vl,d,vr   = d.U, d.S, d.V;
#     # di        = one(length(d)) ./ d;
#     # xxu       = concatd2(X);
#     # B         = vl'*y;
#     # B         = (vr.*repeat(di',nvar*lags+nx,1))*B;
#     # u         = y - xxu*B;
#     # xxi       = vr.*repeat(di',nvar*lags+nx,1);
#     # xxi       = xxi*xxi';
#     # B         = reshape(B,nvar*lags+nx,nvar); # rhs [variables, equations]

#     # # Normal–Inverse Wishart Priors
#     # M           = u'*u;
#     # Q           = rand(InverseWishart(T-nvar*lags-1,M));
#     # XX          = kron(Q,xxi);
#     # Bnorm       = rand(MvNormal(vec(B), XX));
#     # Bnorm       = Bnorm';
#     # b           = reshape(Bnorm,nvar*lags + nx,nvar);
#     # bb          = b[1:nvar*lags,:];
#     # bbb         = reshape(bb,nvar,lags,nvar);        # variables, lags, equations
#     # bbb         = permutedims(bbb,[3,1,2]);          # equations, variables, lags
#     # Constant    = b[nvar*lags + 1,:];

#     # else
    
#     # Options
#     quarterly   = 1
#     prior       = 1
#     lev         = zeros(nvar,1)
#     if true # Minnesota Priors
#         hm,bm       = minneprc(y,X,lags,quarterly,constant,lev,prior);
#     else # Diffuse Priors
#         hm = bm = 0
#     end    
    
#     xxx         = inv(kron(I(nvar),X'*X) + Matrix{Float64}(I, size(hm,1), size(hm,1)) .* hm);
#     bb          = xxx * (vec(X'*y) + bm);      # stacked posterior means of regression coeffs: vec(bpost) = inv[inv(H) + kron(I,X'X)] * [vec[X'Y] + inv(H)*bprior]
#     b           = reshape(bb,ncoeffs,nvar);    # posterior coefficient matrix 
#     sxx         = cholesky(xxx,check = false).U';
#     sinv        = cholesky(inv(vres),check = false).U;


#     bvar_results = Dict("b"=>b,"sxx"=>sxx,"sinv"=>sinv)


#     # end

#     # randmatrix  = rand(MvNormal(zeros(nvar),Matrix{Float64}(I,nvar,nvar)),nvar);
#     # randvec     = rand(Normal(0,1),nvar*ncoeffs);

#     # Bnorm,Q     = niw(b,sxx,sinv,T,randmatrix,randvec);
#     # b           = reshape(Bnorm,nvar*lags + nx,nvar);
#     # bb          = b[1:nvar*lags,:];
#     # bbb         = reshape(bb,nvar,lags,nvar);        # variables, lags, equations
#     # bbb         = permutedims(bbb,[3,1,2]);          # equations, variables, lags
#     # Constant    = b[nvar*lags + 1,:];

#     # # Calculate Residual and bfo for each draw
#     # ut    = y - xxu*b; # Residuals
#     # bfo   = xxu*b;     # Forcasted variables
#     # F     = b';
#     # s2    = ut'*ut./T;  # Estimated error variance

#     # # Companion Matrix
#     # Fcomp = [F[:,1:nvar*(lags)];  Matrix{Float64}(I,nvar*(lags-1),nvar*(lags-1)) zeros(nvar*(lags-1),nvar)];

#     # # Beta Draws
#     # PhiVAR  = bbb;# # [equations, variables, lags] to match impulsdt.m

#     # bvar_results = Dict("PhiVAR"=>PhiVAR,"b"=>b,"Q"=>Q,"nx"=>nx,"ut"=>ut,"bfo"=>bfo,"Fcomp"=>Fcomp,"B"=>B,"s2"=>s2,"Constant"=>Constant)

#     return bvar_results

# end
