[aud, Fs] = audioread('SAR_Test_File.m4a');
aud = -aud;

sync = aud(:, 2);
back = aud(:, 1);

Tp = 25e-2; 
N = Tp*Fs; 

% Load your sampled down-converted radar data (radar_data)

% Define Trp (minimum range profile time duration) in seconds
Trp = 0.250; % 250 milliseconds

% Calculate Nrp (guard band samples) based on Trp and sampling frequency (Fs)
Nrp = Trp * Fs;

% Determine the number of segments based on the length of radar_data
num_segments = floor(length(back) / Nrp);

% Initialize a cell array to store data matrices for each segment
data_matrices = cell(num_segments, 1);

% Iterate through segments and parse data matrices
for i = 1:num_segments
    % Calculate the range profile start and end indices based on Nrp
    range_profile_start = (i - 1) * Nrp + 1;
    range_profile_end = i * Nrp;
    
    % Extract the range profile data for this segment
    range_profile_data = back(range_profile_start:range_profile_end);
    
    % Store the range profile data in the data matrix
    data_matrices{i} = range_profile_data;
end

Tp1=20e-3
N=Tp1*Fs

time_axis = linspace(0, Trp, N);

% Initialize a cell array to store integrated data for each segment
integrated_data = cell(length(data_matrices), 1);

% Iterate through the data_matrices and perform integration for each segment
for i = 1:length(data_matrices)
    % Check if the data_matrices{i} is empty or has less than Nrp/2 elements
    if ~isempty(data_matrices{i}) && size(data_matrices{i}, 2) >= N/2
        % Sum all the up-chirp data for this segment
        integrated_segment = sum(data_matrices{i}(:, 1:N/2), 2);
        
        % Divide by the number of up-chirps (assuming you know this number)
        % Replace 'num_up_chirps' with the actual number of up-chirps
        num_up_chirps = 2221; 
        integrated_segment = integrated_segment / num_up_chirps;
        
        % Store the integrated data in the cell array
        integrated_data{i} = integrated_segment;
        
       
    end
end





% Define constants (replace with actual values)
Tp = 20; % Range profile time in seconds
fstart = 1e9; % Start frequency in Hz
fstop = 2e9; % Stop frequency in Hz
c = 3e8; % Speed of light in m/s

% Calculate the frequency range corresponding to Kr values
Kr_start = (4 * pi / c) * fstart;
Kr_stop = (4 * pi / c) * fstop;

% Iterate through the data_matrices and perform Hilbert transform and other operations
for i = 1:length(data_matrices)
    % Check if the data_matrices{i} is not empty and has data
    if ~isempty(data_matrices{i})
        % Step a: Perform IFFT on each row
        ifft_result = ifft(data_matrices{i}, 4*N, 2);

        % Step b: IFFT on the positive frequencies of each row
        % Assuming that positive frequencies correspond to the first N/2 elements
        ifft_result(:, Nrp/2+1:end) = 0; % Set the negative frequencies to zero
        inverse_fft_result = ifft(ifft_result, [], 2);

        % Step c: Replace NaN values with 1e-30
        inverse_fft_result(isnan(inverse_fft_result)) = 1e-30;

        
       
    end
end

plot (inverse_fft_result(2000));


% Define constants (replace with actual values)
fstart = 1e9; % Start frequency in Hz
fstop = 2e9; % Stop frequency in Hz
c = 3e8; % Speed of light in m/s


% Calculate Kr matrix
Kr_start = (4 * pi / c) * fstart;
Kr_stop = (4 * pi / c) * fstop;
Kr = linspace(Kr_start, Kr_stop, N/2); 

% Calculate Xa matrix
num_positions = length(data_matrices);
L = N/4;
Xa = linspace(-L / 2, L / 2, L); 








