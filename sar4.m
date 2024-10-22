clear;
clc;
[aud, Fs] = audioread('SAR_Test_File.m4a');
aud = -aud;

sync = aud(:, 2);
back = aud(:, 1);

Trp = 0.250;
Nrp = Trp * Fs;

threshold = 0.5;
sync_clean = sync > threshold;

%getting different indices for different positions:
n = 25;
c = 0;
ind = 1;
pos = zeros(n, 1);
for i = 1:length(sync_clean)
    if sync_clean(i) == 0
        c = c+1;
    else 
        if (c >= 5000)  || ind == 1
            pos(ind) = i;
            ind = ind + 1;
        end
        c = 0;
    end
end

back_psn_parse = zeros(n, Nrp);
sync_psn_parse = zeros(n, Nrp);

for i = 1:25
    back_psn_parse(i, : ) = back(pos(i): pos(i) + Nrp - 1);
    sync_psn_parse(i, : ) = sync_clean(pos(i): pos(i) + Nrp - 1);
end

% Plots for Position Parsed data
% plot(back_psn_parse(1, :))
% title("Position Parsed down-converted data")
% xlabel("Data Sample Number (Nrp)")
% ylabel("Amplitude")
% 
% figure;
% 
% plot(sync_psn_parse(1,:))
% title("Position parsed sync data")
% xlabel("Data Sample Number (Nrp)")
% ylabel("Amplitude")

%Integration

Tp = 0.02;
N = Tp*Fs;
back_uc_parsed = zeros(n,N/2);
uc_ind = [];
a = 1;
for i = 1:n
    for j = 1:Nrp-1
        if sync_psn_parse(i,j+1) - sync_psn_parse(i,j) > 0 && a<=5
            uc_ind(a) = j;
            a = a+1;
        end
    end
    n_uc = length(uc_ind);
    %Performing integration over 5 up-chirps per position
    row_i = back_psn_parse(i, : );
    for k = 1:N
        back_uc_parsed(i,k) = sum(row_i(uc_ind + k -1))/n_uc;
    end
end

Tp = 0.02;
N = Tp * Fs;
ifft_data = zeros(n, N/2);
fft_data = [];

% Iterate through the data_matrices and perform Hilbert transform and other operations
for i = 1:n
    % Step 1: Perform FFT on each row
    fft_result = fft(back_uc_parsed(i,:));
    %fft_result = fftshift(fft_result);
    % Step 2: ifft of the second half of the vector
   
    ifft_result = ifft(fft_result(N/2+1:N));

    % Step 3: Replace NaN values with 1e-30
    ifft_result(isnan(ifft_result)) = 0;
    ifft_data(i, :) = ifft_result;
end

f_start = 2.408e9;
f_stop = 2.495e9;
fc = 2.43e9;
c = 3e8;
lambda = c/fc;

% Calculate Kr matrix
Kr_start = (4 * pi / c) * f_start;
Kr_stop = (4 * pi / c) * f_stop;
Kr = linspace(Kr_start, Kr_stop, N/2); 

%Xa array
Xa = lambda/2*linspace(-(n-1)/2,(n-1)/2,n);
L = N/2;

%range of values for Hann Window
Kx = linspace(-L / 2, L / 2, N/2);

%Hann window
H = 0.5*(1+cos(2*pi*Kx/L));

% plot(H);
% xlabel('Data Sample Number');
% ylabel('Amplitude');
% title('Hann Window (time domain)')

ifft_hann = zeros(n,N/2);
for i = 1:n
    ifft_hann(i, :) = H.*ifft_data(i, :); 
end

%Phase of Hilbert and Hann transformed row-wise inverse fft
phase_before = -angle(ifft_hann); 

%Plot for phase before along-track FFT
imagesc(Kr, Xa, phase_before);
caxis([-3 3]);
colorbar
title('phase before along track FFT');
xlabel('Kr(rad/m)');
ylabel('synthetic apeture position, Xa(m)');

find_NaN = isnan(ifft_hann);
ifft_hann(find_NaN) = 0;

zeropad = zeros(1012, N/2);
iffth_zp = [zeropad; ifft_hann; zeropad];

f_shift = fftshift(fft(iffth_zp), 1);
size_fshift = size(f_shift);
Kx = linspace(-2*pi/lambda, 2*pi/lambda, size_fshift(1));

figure(2)
imagesc(Kr, Kx, 20*log10(abs(f_shift)));
caxis([-25 10]);
colorbar;
title('Magnitude after along track FFT');
xlabel('Kr(rad/m)');
ylabel('Kx(rad/m)');

phase_after = -angle(f_shift);

figure(3)
imagesc(Kr, Kx, phase_after);
caxis([-3 3]);
colorbar;
title('Phase After Along Track FFT');
xlabel('Kr (rad/m)');
ylabel('Kx (rad/m)');

%number of rows in fftshift data matrix for the length of Kx and Ky
P = length(Kx);
Ky = zeros(P, length(Kr));

for i = 1:1:P
    Ky(i, :) = sqrt(Kr(1,:).^2-Kx(i).^2);
end

%Creating Kye vector
kye_start = floor(min(Ky, [], 'all'));
kye_end = floor(max(Ky, [], 'all'));
Kye = linspace(kye_start, kye_end, floor(P/2));

data_interpol = [];
for i = 1:1:P
    data_interpol(i, :) = interp1(Ky(i,:),f_shift(i,:),Kye);
end

find_NaN = isnan(data_interpol);
data_interpol(find_NaN) = 0;  

phase_si = -angle(data_interpol);

figure(4)
imagesc(Kye, Kx, phase_si);
caxis([-3 3]);
colorbar;
title('phase after Stolt interpolation');
xlabel('Kr(rad/m)');
ylabel('Kx(rad/m)');

% 4xP zero padding of interpolated data
pf = 4;
data_ip_padded = ifft2(data_interpol, 4*P, 2*P);
interpol_data_mtrx = rot90(data_ip_padded);
interpol_data_mtrx = fliplr(interpol_data_mtrx);

d_range_1 = 1;
d_range_2 = 100;
c_range_1 = -25;
c_range_2 = 25;
delta_x = 0.1;

del_fy = c*(kye_end - kye_start)/(4*pi);
del_ky = (4*pi*del_fy)/c;
Rmax =  (2*P*c)/(2*del_fy);
Rail_Rmax = P*delta_x;

d_index_1 = round((size(interpol_data_mtrx,1)/Rmax)*d_range_1);
d_index_2 = round((size(interpol_data_mtrx,1)/Rmax)*d_range_2);
c_index_1 = round((size(interpol_data_mtrx,2)/Rail_Rmax)*(c_range_1+(Rail_Rmax/2)));
c_index_2 = round((size(interpol_data_mtrx,2)/Rail_Rmax)*(c_range_2+(Rail_Rmax/2)));

[p,q] = size(interpol_data_mtrx);
sub_matrix = interpol_data_mtrx(1:p/2-p/8, :);
[row,col] = size(sub_matrix);

cross_range = linspace(c_range_1, c_range_2, col);
down_range = linspace(-d_range_1, -d_range_2, row);

down_abs_colv = abs(down_range).^2;
fin_mtrx = zeros(row, col);
for i = 1:col
    fin_mtrx(:, i) = down_abs_colv'.*sub_matrix(:,i);
end
fin_mtrx = 20*log10(abs(fin_mtrx));
figure(5)
imagesc(cross_range, down_range, fin_mtrx);
caxis([-60 -10]);
%xlim = ([-80,80]);
colorbar;
%caxis([-3 3]);
title('Final image');
xlabel('Crossrange (meter)');
ylabel('Downrange (meter)');