function [ dir ] = fn_make_dir( dir )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% Check to see if directory exists and if not, make one
if ~exist(dir,'dir')
    mkdir(dir)
end

end

