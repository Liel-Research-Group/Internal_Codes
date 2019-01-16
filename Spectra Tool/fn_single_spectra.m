function [ spectra ] = fn_single_spectra( dt, ag )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Define Period Range that get more coarse as period increases
T = [0.01:0.01:1,1.02:.02:2,2.05:.05:3,3.1:.1:5];
damp_ratio = [0.01, 0.02, 0.03, 0.05];

%% Run Single Spectra
psa = zeros(length(T),length(damp_ratio));
for i = 1:length(damp_ratio)
    for j = 1:length(T)
        [psuedoAccelerationTH, ~, ~] = fn_sdof_th(T(j), damp_ratio(i), ag, dt);
        psa(j,i) = max(abs(psuedoAccelerationTH));    
    end
end

% Save Spectra
id = 1:length(T);
spectra = table(id',T',psa(:,1),psa(:,2),psa(:,3),psa(:,4),'VariableNames',{'id','period','psa_1','psa_2','psa_3','psa_5'});

end

