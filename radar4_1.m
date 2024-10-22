clear;
[aud, Fs] = audioread('SAR_Test_File.m4a');
aud = -aud;

sync = aud(:, 2);
back = aud(:, 1);

Trp = 2;
Nrp = Trp * Fs;

M = floor(length(back) / Nrp);
data_matrices = cell(M, 1);

for i = 1:M
    rstart = (i - 1) * Nrp + 1;
    rend = i * Nrp;
    
    data = back(rstart:rend);
    
    data_matrices{i} = data;
end

plot(data_matrices{5});
xlabel('Data Sample Number');
ylabel('Amplitude');
title('Position parsed down-converted data');

figure;

sync_matrices = cell(M, 1);

for i = 1:M
    rstart = (i - 1) * Nrp + 1;
    rend = i * Nrp;
    
    data_sync = sync(rstart:rend);
    data_sync(data_sync < 0) = 0;
    data_sync(data_sync > 0) = 1;

    sync_matrices{i} = data_sync;
end


plot(sync_matrices{5});
xlabel('Data Sample Number');
ylabel('Amplitude');
title('Position parsed sync data');


figure;

Tp = 0.;
N = Tp * Fs;

% Initialize a cell array to store integrated data for each segment
integrated_data = cell(M, 1);
% 
% % Iterate through the data_matrices and perform integration for each segment
% % for i = 1:M
% %     % Check if the data_matrices{i} is not empty
% %     if ~isempty(data_matrices{i}) 
% %                 
% %         % Sum only the up-chirp data for this segment
% %         integrated_segment = sum(data_matrices{i}(1:(Nrp+1)/2, :));
% %        
% %         % Store the integrated data in the cell array
% %         integrated_data{i} = integrated_segment;
% %     end
% % end
% 

ifft_data = cell(M, 1);

% % Iterate through the data_matrices and perform Hilbert transform and other operations
% for i = 1:M
%     % Step 1: Perform FFT on each row
%     fft_result = fft(integrated_data{i}, 4*N, 2);
% 
%     % Step 2: ifft of the second half of the vector
%    
%     ifft_result = ifft(fft_result(floor(N/2)+1:N));
% 
%     % Step 3: Replace NaN values with 1e-30
%     ifft_result(isnan(ifft_result)) = 1e-30;
% 
%     ifft_data{i} = ifft_result;
% end

% plot(ifft_result);
% figure;
% 
% fstart = 1e9; % Start frequency in Hz
% fstop = 2e9; % Stop frequency in Hz
% c = 3e8; % Speed of light in m/s
% 
% % Calculate Kr matrix
% Kr_start = (4 * pi / c) * fstart;
% Kr_stop = (4 * pi / c) * fstop;
% Kr = linspace(Kr_start, Kr_stop, N/2); 
% 
% % Calculate Xa 
% L = N/2;
% Xa = linspace(-L / 2, L / 2, L);
% 
% %Hann window
% H = 0.5*(1+cos(2*pi*Xa/L));
% plot(H);
% xlabel('Data Sample Number');
% ylabel('Amplitude');
% title('Hann Window (time domain)')
% figure;
% 
% for i= 1:M
%     ifft_data{i} = ifft_data{i} .* H;
% end
% 
% %Plots
% for i = 1:M
%     imagesc(Kr, Xa, angle(ifft_data{i}));
%     clim([-3, 3]);  
%     colorbar;
%     title("SAR");
% end

% Assuming sync_pulse_data contains your sync and pulse data
threshold = 0.1; % Set your threshold value
binary_signal = data_sync > threshold; % Thresholding

% Detect rising edges (starts of upchirp cycles)
rising_edges = diff(binary_signal) == 1;


num_upchirps = sum(rising_edges);

% Assuming rising_edges is the variable containing the rising edge indices

% Find the positions of upchirps
upchirp_positions = find(rising_edges);

% Now, upchirp_positions contains the indices where each upchirp starts in your data






% Assuming upchirp_positions contains the indices where each upchirp starts
% and sync_pulse_data is your sync and pulse data

% Initialize variables
num_upchirps = length(upchirp_positions);
integrated_data = zeros(1, num_upchirps);

% Iterate through upchirps and integrate data within each upchirp
for i = 1:num_upchirps
    % Extract upchirp cycle data
    if i < num_upchirps
        upchirp_data = data_sync(upchirp_positions(i):upchirp_positions(i+1)-1);
    else
        % For the last upchirp, extract until the end of the data
        upchirp_data = data_sync(upchirp_positions(i):end);
    end
    
    % Integrate the upchirp data
    integrated_data(i) = trapz(upchirp_data);
end

% Calculate the average integrated value across all upchirps
average_integrated_value = mean(integrated_data);


% Define Trp (minimum range profile time duration) in seconds and sampling frequency Fs
Trp = 0.020; % 20ms
Fs = 1000; % Replace with your actual sampling frequency in Hz

% Calculate N based on Trp and Fs
N = Trp * Fs;

% Initialize the data matrix with zeros
data_matrix = zeros(N - 1, N / 2);

% Assuming upchirp_positions contains the indices where each upchirp starts
% and sync_pulse_data is your sync and pulse data

% Iterate through positions and integrate up-chirp data
for i = 1:(N - 1)
    % Find upchirp positions for this position i
    upchirp_indices = upchirp_positions(upchirp_positions >= (i - 1) * N / 2 & upchirp_positions < i * N / 2);

    % Extract and integrate up-chirp data for this position
    integrated_value = 0;
    for j = 1:length(upchirp_indices) - 1
        upchirp_data = sync_pulse_data(upchirp_indices(j):upchirp_indices(j + 1) - 1);
        integrated_value = integrated_value + trapz(upchirp_data);
    end

    % Calculate average integrated value for this position
    average_integrated_value = integrated_value / (length(upchirp_indices) - 1);

    % Populate the data matrix
    data_matrix(i, :) = average_integrated_value;
end

% Assuming data_matrix is your processed data matrix

% Step 1: IFFT on each row (position)
ifft_data = ifft(data_matrix, [], 2);

% Step 2: IFFT on the positive frequencies of each row
positive_frequencies_data = ifft_data(:, 1:end/2);

% Step 3: Replace NaN values with 1e-30
positive_frequencies_data(isnan(positive_frequencies_data)) = 1e-30;

