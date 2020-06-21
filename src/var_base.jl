# THIS IS A COPY OF RFVAR.M WITH THE MINNESOTA PRIOR ADDED
# This algorithm goes for acuturacy without worrying about memory requirements.
# ydata:   dependent variable data matrix
# xdata:   exogenous variable data matrix
# lags:    number of lags

function var_base(TS,options_var)

 # prior,lag,constant = options.options_var["prior"],options.options_var["lag"]

 # for ii in 1:length(options.options_var); eval(:($(Symbol(options.options_var[ii])) = $options.options_var[ii])); end

 function concatd2(XX)
   XX = dropdims(reshape(XX,size(XX)[1],:,1);dims = 3);
   return XX
end

unpack(options_var)



ydata    = values(TS)
(T,nvar) = size(ydata)

y = ydata[lags+1:end,:];
X = zeros(T-lags,nvar,lags);
for i = 1 : lags
   X[:,:,i] = ydata[lags+1-i:end-i,:];
end




if  constant
    xdata           = ones(T,1);
else
    xdata           = [];
end

nox             = isempty(xdata);

if ~ nox
   (T2,nx) = size(xdata);
else
   T2 = T; nx = 0;
end

# [T2,nx]=size(xdata);
# note that x must be same length as y, even though first part of x will not be used.
# This is so that the lags parameter can be changed without reshaping the xdata matrix.
if T2 != T; println("Mismatch of x and y data lengths");end

# Minnesota Prior
#--------------------------------------------------------------------------
pi1, pi3  = prior[1],prior[2]
if pi1 > 0 && pi3 > 0;
   Xp  = zeros(nvar*lags,nvar,lags);
   #pi1 = 0.2; pi3 = 1;                             
   ysd = X[1,:,:]; ysd = diag(std(ysd,1,3));
   for i = 1 : nvar;
      y = [pi1*ysd[nvar+1-i,:]; zeros(lags-1,nvar); y];
      Xp[lags*(i-1)+1,:,1] = pi1*ysd[i,:];
      for j = 1 : lags-1;
         Xp[lags*(i-1)+(j+1),:,j+1]=pi1*ysd[i,:]*(j+1)^pi3;
      end
   end

X = [X ; Xp];

if ~nox
   (T2,nx) = size(xdata);
else
   T2 = T; nx = 0; # xdata=zeros(T2,0);
end


if constant
xdata = [xdata[lags+1:end,:]; zeros(nvar*lags,nx)];
end
X     = [concatd2(X) xdata];
end
#--------------------------------------------------------------------------

if  pi1 == 0 && pi3 == 0;
    X = [concatd2(X) xdata[lags+1:T]];
end

d         = svd(concatd2(X),full=false);
vl,d,vr   = d.U, d.S, d.V;
di        = one(length(d)) ./ d;
xxu       = concatd2(X);
B         = vl'*y;
B         = (vr.*repeat(di',nvar*lags+nx,1))*B;
u         = y - xxu*B;
xxi       = vr.*repeat(di',nvar*lags+nx,1);
xxi       = xxi*xxi';
B         = reshape(B,nvar*lags+nx,nvar); # rhs [variables, equations]

var_results = Dict("B" => B, "u" => u, "xxu" => xxu,"xxi" => xxi, "nx" => nx)

return var_results
end



function var_base_ols(TS,options_var)

ydata    = values(TS)
(T,nvar) = size(ydata)

y = ydata[lags+1:end,:];
X = zeros(T - lags,nvar * lags);
for lags_count in 1:lags
   X[:,1+(lags_count-1)*nvar:lags_count*nvar] = ydata[(lags+1-lags_count):(T-lags_count),:];
end

if constant
   X = [X ones(T - lags,1)];
end

# OLS VAR
Beta = X \ y;
Res  = y - X * Beta;          
vres = (1 / (T - lags)) * (Res' * Res);
stds = sqrt.(diag(vres));

var_results = Dict("Beta" => Beta, "Res" => Res, "vres" => vres, "stds" => stds)

return var_results
end
