function [ stress, direction ] = fn_principle_stress_3d( sig_x, sig_y, sig_z, tau_xy, tau_yz, tau_zx )
% FUNCTION DESCRIPTION
% This Function perfroms an eigen analysis of a 3x3 matrix in order to
% solve for the principle stresses on a loaded body.

%   Created By; Dustin Cook
%   Date Created: 09/22/2017

% FUNCTION INPUTS 
%   sig_ - Normal stresses in the x y and z directions
%   tau_ - Shear stresses in the x y and z directions

% FUNCTION OUTPUTS
%   stress - Magnitudes of three priciple stress in decending order. Same
%            units as input stress.
%   direction - Unit vectoru directions of three principle stresses in the
%               same order as the stess values.

% ASSOCIATED PACKAGES
% 

% NOTES
% 

%% Begin Method
% Assmeble Matrix
A = [sig_x, tau_xy, tau_zx; tau_xy, sig_y, tau_yz; tau_zx, tau_yz, sig_z];

% Run analysis
[V,D] = eig(A);
priciple_stress = diag(D);

% Arrange Principle Stresses
stress = cell(3,1);
direction = cell(3,1);
for i = 1:length(priciple_stress)
    if priciple_stress(i) == max(priciple_stress)
        stress{1} = priciple_stress(i);
        direction{1} = V(:,i);
    elseif priciple_stress(i) == min(priciple_stress)
        stress{3} = priciple_stress(i);
        direction{3} = V(:,i);
    else
        stress{2} = priciple_stress(i);
        direction{2} = V(:,i);
    end
end

end

