%Eric Johnson
%Liel Research Group
%NGA-West2 Equations for Predicting PGA, PGV, and 5% Damped PSA for Shallow Crustal Earthquakes
%David M. Boore, Jonathan P. Stewart, Emel Seyhan, and Gail M. Atkinson
%
%This GMPE is dependent on the following files:
%   BSSA_14.mat
%
%Inputs required are:
%   siteprop.Rjb, siteprop.VS30, siteprop.T, faultprop.M, and faultprop.d

function [Sa,SD] = BSSA_14(siteprop,faultprop)

T = siteprop.T;

coeff = load('BSSA_14.mat');
T_vec = coeff.T_vec;

if isempty(find(abs(T_vec-T)<0.0001)) == 0
    
    [Sa,SD] = BSSA_14_sub1(T,siteprop,faultprop);
    
else
    
    T_low = T_vec(max(find(T_vec<T)));
    T_high = T_vec(min(find(T_vec>T)));
    
    [Sa_low, SD_low] = BSSA_14_sub1( T_low,siteprop,faultprop );
    [Sa_high, SD_high] = BSSA_14_sub1( T_high,siteprop,faultprop );
    
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

function [Sa,SD] = BSSA_14_sub1(T,siteprop,faultprop)

Rjb = siteprop.Rjb;
VS30 = siteprop.VS30;

M = faultprop.M;
d = faultprop.d;

coeff = load('BSSA_14.mat');
T_vec_sub = coeff.T_vec;

c = interp1(T_vec_sub,coeff.c_vec,T);
c1 = interp1(T_vec_sub,coeff.c1_vec,T);
c2 = interp1(T_vec_sub,coeff.c2_vec,T);
c3 = interp1(T_vec_sub,coeff.c3_vec,T);
Dc3 = interp1(T_vec_sub,coeff.Dc3_vec,T);
DfR = interp1(T_vec_sub,coeff.DfR_vec,T);
DfV = interp1(T_vec_sub,coeff.DfV_vec,T);
e0 = interp1(T_vec_sub,coeff.e0_vec,T);
e1 = interp1(T_vec_sub,coeff.e1_vec,T);
e2 = interp1(T_vec_sub,coeff.e2_vec,T);
e3 = interp1(T_vec_sub,coeff.e3_vec,T);
e4 = interp1(T_vec_sub,coeff.e4_vec,T);
e5 = interp1(T_vec_sub,coeff.e5_vec,T);
e6 = interp1(T_vec_sub,coeff.e6_vec,T);
f1 = interp1(T_vec_sub,coeff.f1_vec,T);
f3 = interp1(T_vec_sub,coeff.f3_vec,T);
f4 = interp1(T_vec_sub,coeff.f4_vec,T);
f5 = interp1(T_vec_sub,coeff.f5_vec,T);
f6 = interp1(T_vec_sub,coeff.f6_vec,T);
f7 = interp1(T_vec_sub,coeff.f7_vec,T);
h = interp1(T_vec_sub,coeff.h_vec,T);
Mh = interp1(T_vec_sub,coeff.Mh_vec,T);
Mref = interp1(T_vec_sub,coeff.Mref_vec,T);
phi1 = interp1(T_vec_sub,coeff.phi1_vec,T);
phi2 = interp1(T_vec_sub,coeff.phi2_vec,T);
R1 = interp1(T_vec_sub,coeff.R1_vec,T);
R2 = interp1(T_vec_sub,coeff.R2_vec,T);
Rref = interp1(T_vec_sub,coeff.Rref_vec,T);
tau1 = interp1(T_vec_sub,coeff.tau1_vec,T);
tau2 = interp1(T_vec_sub,coeff.tau2_vec,T);
V1 = interp1(T_vec_sub,coeff.V1_vec,T);
V2 = interp1(T_vec_sub,coeff.V2_vec,T);
Vc = interp1(T_vec_sub,coeff.Vc_vec,T);
Vref = interp1(T_vec_sub,coeff.Vref_vec,T);

%Source Function
if M <= Mh
    FE = e0 + e4*(M - Mh) + e5*(M - Mh)^2;
else
    FE = e0 + e6*(M - Mh);
end

%Path Function
R = sqrt(Rjb^2 + h^2);
FP = (c1 + c2*(M - Mref)) * log(R/Rref) + (c3 + Dc3)*(R - Rref);

%Site Function
if VS30 <= Vc
    F_lin = c*log(VS30/Vref);
else
    F_lin = c*log(Vc/Vref);
end

siteprop_r = siteprop;
faultprop_r = faultprop;

if T ~= 0 || VS30 ~= 760
    siteprop_r.T = 0;
    siteprop_r.VS30 = 760;
    [PGAr,~] = BSSA_14_sub1(0,siteprop_r,faultprop_r);
    f2 = f4*(exp(f5*(min(VS30,760)-360)-360)-exp(f5*(760-360)));
    F_nl = f1 + f2*log((PGAr + f3)/f3);
else
    F_nl = 0;
end

FS = F_lin + F_nl;

%Aleatory Uncertainty
if M <= 4.5
    PhiM = phi1;
elseif M > 4.5 && M < 5.5
    PhiM = phi1 + (phi2 - phi1)*(M - 4.5);
elseif M >= 5.5
    PhiM = phi2;
end

if Rjb <= R1
    Phi_MRjb = PhiM;
elseif Rjb > R1 && Rjb <= R2
    Phi_MRjb = PhiM + DfR*(log(Rjb/R1)/log(R2/R1));
elseif Rjb > R2
    Phi_MRjb = PhiM + DfR;
end

if VS30 >= V2
    phi = Phi_MRjb;
elseif VS30 >= V1 && VS30 <=V2
    phi = Phi_MRjb + DfV*(log(V2/VS30)/log(V2/V1));
elseif VS30 <= V1
    phi = Phi_MRjb + DfV;
end

if M <= 4.5
    tau = tau1;
elseif M > 4.5 && M < 5.5
    tau = tau1 + (tau2 - tau1)*(M - 4.5);
elseif M >= 5.5
    tau = tau2;
end 

sigma = sqrt(phi^2 + tau^2);

lnSa = FE + FP + FS;
Sa = exp(lnSa);

SD = [sigma;phi;tau];

end