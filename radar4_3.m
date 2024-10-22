clear;
[aud, Fs] = audioread('SAR_Test_File.m4a');
aud = -aud;

sync = aud(:, 2);
back = aud(:, 1);

Trp = 2;
Nrp = Trp*Fs;

treshold=0.5;
sync_data= sync > treshold;

xn = zeros(25,1);
n = length(xn);
index = 1;
counter = 0;
for i=1:length(sync_data)
    if sync_data(i) == 0
        counter = counter + 1;
    else
        if counter > 2000 || index == 1
            xn(index) = i;
            index = index + 1;
        end
        counter = 0;
    end
end


A = eye(n,Nrp);
B = eye(n,Nrp);
for i=1:n
    A(i,:) = back(xn(i):xn(i)+Nrp-1);
    B(i,:) = sync_data(xn(i):xn(i)+Nrp-1);
end


Tp = 0.02;
N = Tp*Fs;
C = eye(n,N);
for i=1:n
    %plot(B(i,:), '--');
    %hold on
    rising_edges = find(diff(B(i,:)) > 0) + 1;
    up_chirps = intersect(rising_edges(find(rising_edges>1000)),rising_edges(find(rising_edges<10000)));
    l = length(up_chirps);
    index = 2;
    while l > 5
        if up_chirps(index) - up_chirps(index-1) < 2*N
            up_chirps(index) = [];
            l = l - 1;
        else
            index = index + 1;
        end
    end
    N_up_chirps = length(up_chirps);
    %scatter(up_chirps, ones(1,length(up_chirps)), 'r*');
    % integration
    for j=1:N
        %disp(A(i,up_chirps+(j-1)));
        C(i,j) = sum(A(i,up_chirps+(j-1)))/N_up_chirps;
    end
end



F = fft(C,N,2);
IF = ifft(F(1:n,1:N/2),N/2,2);


f_start = 2.408e9;
f_stop = 2.495e9;
fc = 2.43e9;
c = 3e8;
lamda = c/fc;
Kr = 4*pi/c*linspace(f_start,f_stop,N/2);
Xa = lamda/2*linspace(-12,12,length(xn));
L = N/2;
x = linspace(-L/2,L/2,N/2);
H = 1/2*(1+cos(2*pi*x/L));
IFH = H.*IF;

find_NaN = isnan(IFH);
IFH(find_NaN) = 1*(1e-30);



