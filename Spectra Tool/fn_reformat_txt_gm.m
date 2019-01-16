function [ output_txt ] = fn_reformat_txt_gm( file_name, input_dir, output_dir )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Load GM File
txt_data = fileread([input_dir filesep file_name]);

data_raw = str2double(strsplit(txt_data,' '));
data = data_raw(~isnan(data_raw));

% Find G Conversion
g_factor = 980.665;

% Find dt
dt = 0.01;

% Crop data to just Earthquake and translate to data
eq_data = data(1:end);
eq_data_g = eq_data/g_factor;

% Loop through data to create single line signal for running in opensees
fileID = fopen([output_dir filesep erase(file_name,'.txt') '.tcl'],'w');
for j = 1:length(eq_data_g)
    fprintf(fileID,'%d \n',eq_data_g(j));
end
fclose(fileID);

% Loop through data to create single line signal with g-factor and dt as txt doc
output_txt = [g_factor,dt,eq_data_g];
fileID = fopen([output_dir filesep file_name],'w');
for j = 1:length(output_txt)
    fprintf(fileID,'%d \n',output_txt(j));
end
fclose(fileID);

end

