[aud, Fs] = audioread('SAR_Test_File.m4a');
aud = -aud;

sync = aud(:, 2);
back = aud(:, 1);

Tp = 25e-2; 
N = Tp*Fs; 


% Define Trp (minimum range profile time duration) in seconds
Trp = 0.250; % 250 milliseconds

% Calculate Nrp (guard band samples) based on Trp and sampling frequency (Fs)
Nrp = Trp * Fs;

% Determine the number of positions based on the length of the radar_data
num_positions = floor(length(back) / Nrp);

% Initialize a cell array to store data matrices for each position
data_matrices = cell(num_positions, 1);

% Iterate through positions and parse data matrices
for i = 1:num_positions
    % Calculate the range profile start and end indices based on Nrp
    range_profile_start = (i - 1) * Nrp + 1;
    range_profile_end = i * Nrp;
    
    % Extract the range profile data for this position
    range_profile_data = back(range_profile_start:range_profile_end);
    
    % Store the range profile data in the data matrix
    data_matrices{i} = range_profile_data;
end

% Assuming you have already generated data_matrices as described before

% % Define time axis for the plots (assuming uniform time intervals)
% time_axis = linspace(0, Trp, Nrp);
% 
% % Iterate through data_matrices and create plots for each position
% for i = 1:length(data_matrices)
%     % Create a new figure for each position
% 
%     
%     % Plot the radar data for this position
%     plot(time_axis, data_matrices{i});
%     
%     % Customize the plot labels and title
%     xlabel('Time (s)');
%     ylabel('Amplitude');
%     title(['Converted Data for Position ' num2str(i)]);
%     
%     % You can add more customization as needed (e.g., axis limits, legends)
% end

% Assuming data_matrices is your cell array containing sync pulse data

% Assuming data_matrix is your sync pulse data for one position, stored in a variable

% Threshold to detect gaps (adjust this threshold based on your data)
threshold = 0.1; % Change this value according to your signal amplitude

% Apply median filter to smooth the data
smoothed_data = medfilt1(data_matrices, 5); % Adjust filter size as needed

% Apply thresholding operation element-wise to the smoothed data
below_threshold = smoothed_data < threshold;

% Find the start and end indices of the signal
signal_start_index = find(~below_threshold, 1, 'first');
signal_end_index = find(~below_threshold, 1, 'last');

% Extract the sync signal data within the defined range
sync_signal_data = data_matrix(signal_start_index:signal_end_index);

% Now, sync_signal_data contains the thresholded sync pulse data for one position

