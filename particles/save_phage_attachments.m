function save_phage_attachments(n_phages_vs_time, num_steps, dt, outputFolder)
% savePhageAttachments
% ---------------------------------------------------------------
% Saves the time vector and number of attached phages vs. time
% into ASCII .txt files inside a specified output folder.
%
% Usage:
%   savePhageAttachments(n_phages_vs_time, num_steps, dt, outputFolder)
%
% Inputs:
%   n_phages_vs_time : matrix (num_steps x N) 
%       first column = number of attached phages over time
%   num_steps        : number of time steps
%   dt               : time step size
%   outputFolder     : path to output folder (string)
%
% Example:
%   savePhageAttachments(n_phages_vs_time, num_steps, dt, 'output/')
%
% ---------------------------------------------------------------

    % Ensure output folder exists
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % Construct time vector
    time_vec = (0:num_steps-1) * dt;

    % Extract first column (number of attached phages)
    n_phages_vs_time_vec = n_phages_vs_time(:, 1);

    % Define file paths
    time_file = fullfile(outputFolder, 'time_vec.txt');
    nphage_file = fullfile(outputFolder, 'n_phages_vs_time_vec.txt');

    % Save to ASCII text files
    save(time_file, 'time_vec', '-ascii');
    save(nphage_file, 'n_phages_vs_time_vec', '-ascii');

    % Optional: print confirmation
    fprintf('Saved:\n');
    fprintf('  %s\n', time_file);
    fprintf('  %s\n', nphage_file);
end
