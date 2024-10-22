% Define constants
c = 3e8; % Speed of light in meters/second
startFreq = 2.408e9; % Start frequency in Hz
stopFreq = 2.495e9; % Stop frequency in Hz
N = 1024; % Number of samples

% Read the audio file
filename = 'Range_Test_File.m4a';
[y, fs] = audioread(filename);

% Extract sync data and radar backscatter data
syncThreshold = 0.1; % Adjust this threshold as needed
syncData = y(:, 1); % Assuming sync data is in the first channel
radarData = y(:, 2); % Assuming radar backscatter data is in the second channel

% Data inversion correction
radarData = -radarData;

% Find sync points using threshold
syncIndices = find(syncData > syncThreshold);

% Initialize variables to store results
numChirps = length(syncIndices) - 1;
rangeProfile = zeros(numChirps, N);

% Parse up-chirp data and apply zero-padding
for i = 1:numChirps
    startIndex = syncIndices(i) + 1;
    endIndex = syncIndices(i + 1);
    chirpData = radarData(startIndex:endIndex);
    
    % Apply zero-padding
    chirpData = [chirpData; zeros(3 * N, 1)]; % Zero-padding factor of 4
    
    % Perform inverse FFT
    rangeProfile(i, :) = ifft(chirpData, N);
end

% Calculate range values
deltaR = c / (2 * (stopFreq - startFreq));
Rmax = N * deltaR / 2;
rangeValues = linspace(0, Rmax, N);

% Plot the range profile
figure;
imagesc(rangeValues, 1:numChirps, 20*log10(abs(rangeProfile)));
xlabel('Range (meters)');
ylabel('Chirp Index');
title('Range Profile');
colorbar;

% Display max range
fprintf('Max Range: %.2f meters\n', Rmax);