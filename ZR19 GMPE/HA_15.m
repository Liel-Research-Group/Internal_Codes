%Eric Johnson
%Liel Research Group
%Referenced Empirical Ground-Motion Model for Eastern North America
%Behzad Hassani and Gail M. Atkinson
%
%This GMPE is dependent on the following files:
%   HA_15.mat, BSSA_14.m, and BSSA_14.mat
%
%Inputs required are:
%   siteprop.Rjb, siteprop.VS30, siteprop.T, faultprop.M, and faultprop.d

function [Sa,SD] = HA_15(siteprop,faultprop)

T = siteprop.T;
coeff = load('HA_15.mat');
T_vec = coeff.T_vec;
    
if isempty(find(abs(T_vec-T)<0.0001)) == 0
    
    [Sa,SD] = HA_15_sub1(T,siteprop,faultprop);
    
else
    
    T_low = T_vec(max(find(T_vec<T)));
    T_high = T_vec(min(find(T_vec>T)));
    
    [Sa_low, SD_low] = HA_15_sub1( T_low,siteprop,faultprop );
    [Sa_high, SD_high] = HA_15_sub1( T_high,siteprop,faultprop );
    
    totalSD_low = SD_low(1);
    totalSD_high = SD_high(1);
    
    interSD_low = SD_low(2);
    interSD_high = SD_high(2);
    
    intraSD_low = SD_low(3);
    intraSD_high = SD_high(3);
    
    T_bounds = [log(T_low) log(T_high)];
    Sa_bounds = [log(Sa_low) log(Sa_high)];
    intraSD_bounds = [intraSD_low intraSD_high];
    interSD_bounds = [interSD_low interSD_high];
    totalSD_bounds = [totalSD_low totalSD_high];
    
    Sa = exp(interp1(T_bounds,Sa_bounds,log(T)));
    
    intraSD = interp1(T_bounds,intraSD_bounds,log(T));
    interSD = interp1(T_bounds,interSD_bounds,log(T));
    totalSD = interp1(T_bounds,totalSD_bounds,log(T));
    
    SD = [totalSD;interSD;intraSD];
    
end

end

function [Sa,SD] = HA_15_sub1(T,siteprop,faultprop)

Rjb = siteprop.Rjb;
VS30 = siteprop.VS30;

coeff = load('HA_15.mat');
T_vec_sub = coeff.T_vec;

Sa_bssa = BSSA_14(siteprop,faultprop);

C1 = interp1(T_vec_sub,coeff.C1_vec,T);
C2 = interp1(T_vec_sub,coeff.C2_vec,T);
C3 = interp1(T_vec_sub,coeff.C3_vec,T);
phi = interp1(T_vec_sub,coeff.phi_vec,T);
sigma = interp1(T_vec_sub,coeff.sigma_vec,T);
tau = interp1(T_vec_sub,coeff.tau_vec,T);

logFENA = C1 + C2*Rjb + C3*max(0,log10(min(Rjb,150)/50));

FENA = 10^(logFENA);

Sa = Sa_bssa*FENA;

SD = [sigma;phi;tau];

end