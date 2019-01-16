function [ output_txt ] = fn_reformat_NGA_gm( file_name, input_dir, output_dir )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Load GM File
NGA_data = fileread([input_dir filesep file_name]);

% Validate the units are in G
if ~contains(NGA_data,'ACCELERATION TIME SERIES IN UNITS OF G')
    error('Units found to be something other than G')
end

% Find dt
dt_local_start = strfind(NGA_data,'DT=   ');
dt_local_end = strfind(NGA_data,' SEC');
dt_str = NGA_data((dt_local_start+6):(dt_local_end-1));
% Validate dt is 5 chararters long
if length(dt_str)==5
    dt = str2double(dt_str);
else
    error('Could not find dt')
end

% Crop data to just Earthquake and translate to data
eq_local_start = strfind(NGA_data,'E-');
eq_data_str = NGA_data((eq_local_start(1)-9):end);
eq_data = str2double(strsplit(eq_data_str,' '));
eq_data_g = eq_data(~isnan(eq_data));

% Loop through data to create single line signal for running in opensees
fileID = fopen([output_dir filesep erase(file_name,'.AT2') '.tcl'],'w');
for j = 1:length(eq_data_g)
    fprintf(fileID,'%d \n',eq_data_g(j));
end
fclose(fileID);

% Loop through data to create single line signal with g-factor and dt as txt doc
output_txt = [1,dt,eq_data_g];
fileID = fopen([output_dir filesep erase(file_name,'.AT2') '.txt'],'w');
for j = 1:length(output_txt)
    fprintf(fileID,'%d \n',output_txt(j));
end
fclose(fileID);


end

