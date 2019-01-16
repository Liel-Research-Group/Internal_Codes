clear 
close all
clc

%% Define Inputs Directory
input_dir = [pwd filesep 'inputs'];

%% Run Function for each GM in inputs folder
% Load all files in Input folder and convert them to single column signals
% with dt as the first entry
files = [dir([input_dir filesep '*.AT2']),dir([input_dir filesep '*.tcl']),dir([input_dir filesep '*.txt'])];
for i = 1:length(files)
    % Grab GM signal
    if contains(files(i).name,'.AT2') % NGA file
        [ output_dir ] = fn_make_dir( [pwd filesep 'outputs' filesep erase(files(i).name,'.AT2')] );
        [ gm_data ] = fn_reformat_NGA_gm( files(i).name, input_dir, output_dir );
    elseif contains(files(i).name,'.txt') % Preformatted TCL signal
        [ output_dir ] = fn_make_dir( [pwd filesep 'outputs' filesep erase(files(i).name,'.txt')] );
        [ gm_data ] = fn_reformat_txt_gm( files(i).name, input_dir, output_dir );
    else
        error('Ground Motion Data Structure Not Recognized')
    end
    dt = gm_data(2);
    ag = gm_data(3:end);

    % Calculate Spectra
    [ spectra ] = fn_single_spectra( dt, ag );

    % Save Spectra
    writetable(spectra,[output_dir filesep 'spectra.csv'])

    % Plot Spectra
    figure
    hold on
    plot(spectra.period,spectra.psa_1,'LineWidth',1.5,'DisplayName','1% Damping') 
    plot(spectra.period,spectra.psa_2,'LineWidth',1.5,'DisplayName','2% Damping') 
    plot(spectra.period,spectra.psa_3,'LineWidth',1.5,'DisplayName','3% Damping') 
    plot(spectra.period,spectra.psa_5,'LineWidth',1.5,'DisplayName','5% Damping') 
    xlabel('Period (s)')
    ylabel('PSa (g)')
    plot_name = 'spectra';
    fn_format_and_save_plot( output_dir, plot_name, 1 )
end


