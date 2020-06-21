
function VAR_VD(BetaCoeff,vmat,nvars,nlags,filter_data,frequency)

    # David Gauthier
    # 03 - 2020
    # Bank of England
     
     # Companion matrix of demeaned the VAR 
     M                                          = zeros(nvars*nlags,nvars*nlags);			
     M[1:nvars,:]                               = BetaCoeff[1:nvars*nlags,:]';
     M[nvars+1:nvars*nlags,1:nvars*nlags-nvars] = I(nvars*nlags-nvars);
     vmat                                       = [vmat ; zeros(nvars*(nlags-1),nvars)]';
     
     if !filter_data
     # Solve for variance S / [S = F*S*F' + Q (Hamilton (1994) 10.2.13 )]
     VV  = lyapd(M,I(nvars*nlags),vmat'*vmat);
     VV  = diag(VV[1:nvars,1:nvars]);
     global vv2, V
     vv2 = zeros(nvars,1); VD = VS = Array{Float64}(undef,nvars,nvars);
     for i in 1:nvars
         vx1     = lyapd(M,I(nvars*nlags),vmat[i,:]*vmat[i,:]');
         vx2     = abs.(diag(M*vx1*M'+vmat[i,:]*vmat[i,:]'));
         global VS[:,i]  = vx2[1:nvars];
         global vv2     += vx2[1:nvars];
     end
     
     for i in 1:nvars
         VD[:,i] = VS[:,i] ./ VV;
     end
     
     else
        VD = nothing
        # Dens = 2π * real(integral(@(omega) DensVAR(omega,M,Eps),2π / frequency(2),2π / frequency(1),'ArrayValued',true));
        # Dens = diag(Dens);
        # Dens = Dens(1:nvars);
        
        # VD = []; nshocks = 1;
        # for shk = 1 : nshocks
        #     V  = 2π * real(integral(@(omega) DensLoop(omega,M,q1),2*pi / frequency(2),2π / frequency(1),'ArrayValued',true));
        #     VV = diag(V(1:nvars,1:nvars));
        #     VD = [VD (VV)];
        # end
        
        # VD = VD ./ repmat(Dens,1,nshocks);
        
        return VD*100
    end
end
     
# function Dens = DensLoop(omega,M,q1)
#     S     = (eye(size(M)) - M * exp(-1i.*omega)) \ (q1) ;
#  return   Dens  = ((S) * ctranspose(S))/ (2π;
    
# function Dens = DensVAR(omega,M,vmat)
#     S     = (eye(size(M)) - M * exp(-1i.*omega)) \ vmat' ;
#  return   Dens  = ((S) * ctranspose(S))/ (2π);
    
    
     
     # Old, to remove
#         ngrid    = 512;
#      passband = frequency; 
     
#      # Filter Settings 
#      grid_filter = collect(Float64,0 : ((2*pi) / ngrid) : (2*pi*(1 - .5/ngrid)));
#      filter_gain = zeros(1,ngrid);
#      lowest_periodicity=passband[2];
#      highest_periodicity=passband[1];
#      highest_periodicity=max(2,highest_periodicity); # restrict to upper bound of pi
#      filter_gain[(grid_filter.>=2*pi/lowest_periodicity) .* (grid_filter.<=2*pi/highest_periodicity)] .= 1;
#      filter_gain[(grid_filter.<=-2*pi/lowest_periodicity+2*pi) .* (grid_filter.>=-2*pi/highest_periodicity+2*pi)] .= 1;
     
     
#      S_VAR_grid      = Array{ComplexF64}(undef,ngrid,nvars*nlags*nvars*nlags);
#      for ig in 1 : ngrid     
#          if !Bool(filter_gain[ig])
#          tmps             = zeros(nvars*nlags,nvars*nlags);
#          else
#          tmps             = map(omega -> inv(2*pi) * inv(I(nvars*nlags) - M*exp(-im.*omega)) * (vmat)'*I(nvars)*(vmat) * inv(I(nvars*nlags) - M'*exp(im.*omega)), grid_filter[ig]);
#          end
#          S_VAR_grid[ig,:] = tmps[:]';
#      end
 
#      dd         = [2π*real(ifft(S_VAR_grid[:,ig])) for ig in 1:(nvars*nlags)^2];
#      S_VAR_grid = convert(Array,VectorOfArray(dd));
#      vv         = diag(reshape(S_VAR_grid[1,:],nvars*nlags,nvars*nlags));
#      VD         = Array{Float64}(undef,nvars*nlags,nvars*nlags)
 
#      for ii in 1 : nvars
#          vmat_round = vmat;
#          global eye_shocks = zeros(nvars*1,nvars*1); eye_shocks[ii,ii] = 1;         
     
#      S_VAR_grid   = Array{ComplexF64}(undef,ngrid,nvars*nlags*nvars*nlags);
#      for ig in 1:ngrid                
#          if filter_gain[ig] == 0
#          tmps             = zeros(nvars*nlags,nvars*nlags);
#          else
#          tmps             = map(omega -> inv(2*pi) * inv(I(nvars*nlags) - M*exp(-im.*omega)) * (vmat)'*eye_shocks*(vmat) * inv(I(nvars*nlags) - M'*exp(im.*omega)), grid_filter[ig]);
#          end
#          S_VAR_grid[ig,:] = tmps[:]';
#      end
#      dd         = [2π*real(ifft(S_VAR_grid[:,ig])) for ig in 1:(nvars*nlags)^2];
#      S_VAR_grid = convert(Array,VectorOfArray(dd));
#      VD[:,ii]   = abs.(diag(reshape(S_VAR_grid[1,:],nvars*nlags,nvars*nlags)))./vv;
#      end
#      VD = VD[1:nvars,1:nvars];
     
#  end
 
