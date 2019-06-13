%Eric Johnson
%Liel Research Group
%Ground Motion Model for Small-to-Moderate Earthquakes in Texas, Oklahoma, and Kansas
%Georgios Zalachoris and Ellen M. Rathje
%
%This GMPE is dependent on the following files:
%   ZR_19.mat, HA_15.m, HA_15.mat, BSSA_14.m, and BSSA_14.mat
%
%Inputs required are:
%   siteprop.Rjb, siteprop.VS30, siteprop.T, faultprop.M, and faultprop.d



function [Sa,SD] = ZR_19(siteprop,faultprop)

T = siteprop.T;
coeff = load('ZR_19.mat');
T_vec = coeff.T_vec;
    
if isempty(find(abs(T_vec-T)<=0.000001)) == 0

    [Sa,SD] = ZR_19_sub1(T,siteprop,faultprop);

else

    T_low = T_vec(max(find(T_vec<T)));
    T_high = T_vec(min(find(T_vec>T)));

    [Sa_low, SD_low] = ZR_19_sub1( T_low,siteprop,faultprop );
    [Sa_high, SD_high] = ZR_19_sub1( T_high,siteprop,faultprop );

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

    SD = [intraSD;interSD;totalSD];

end
end

function [Sa,SD] = ZR_19_sub1(T,siteprop,faultprop)

M = faultprop.M;
d = faultprop.d;

Rjb = siteprop.Rjb;
Rhyp = sqrt(Rjb^2 + d^2);
VS30 = siteprop.VS30;


coeff = load('ZR_19.mat');
T_vec_sub = coeff.T_vec;

Sa_ha15 = HA_15(siteprop,faultprop);

Alpha = interp1(T_vec_sub,coeff.Alpha_vec,T);
Rb = interp1(T_vec_sub,coeff.Rb_vec,T);
Mb = interp1(T_vec_sub,coeff.Mb_vec,T);
b0 = interp1(T_vec_sub,coeff.b0_vec,T);
b1 = interp1(T_vec_sub,coeff.b1_vec,T);
c = interp1(T_vec_sub,coeff.c_vec,T);
Vc = interp1(T_vec_sub,coeff.Vc_vec,T);
Cadj = interp1(T_vec_sub,coeff.Cadj_vec,T);
Tau = interp1(T_vec_sub,coeff.Tau_vec,T);
Phi = interp1(T_vec_sub,coeff.Phi_vec,T);
Sigma = interp1(T_vec_sub,coeff.Sigma_vec,T);


%Magnitude Adjustment Factor
if M >= 3.0 && M < Mb
    F_M = b0;
elseif M >= Mb && M < 5.8
    F_M = b0 + b1 * (M - Mb);
else
    fprintf('This GMPE is only valid for 3.0 <= Mw <= 5.8')
end

%Near Source Attenuation Correction
if Rhyp < 4 %km
    F_R = Alpha * log(4 / Rb);
elseif Rhyp >= 4 && Rhyp < Rb
    F_R = Alpha * log(Rhyp / Rb);
else
    F_R = 0;
end

%Site Amplification Correction
if VS30 < Vc
    F_S = c * log(VS30 / Vc);
else
    F_S = 0;
end

%Overall Adjustment
Cadj = Cadj;

lnF = Cadj + F_M + F_R + F_S;
F = exp(lnF);

Sa = Sa_ha15*F;


SD = [Phi;Tau;Sigma];

end