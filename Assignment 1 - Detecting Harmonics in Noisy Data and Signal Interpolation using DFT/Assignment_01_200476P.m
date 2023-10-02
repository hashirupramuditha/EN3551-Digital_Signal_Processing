%% Assignment 1 - Detecting Harmonics in Noisy Data and Signal Interpolation using DFT
% Name          : Pramuditha A.A.H.
% Index No.     : 200476P


%% 3.1 Harmonic Detection

clc;
clear;
close all;

Fs = 128;            % Sampling rate
start_time = 0;      % Starting time of the sampling
end_time = 14;       % Ending time of the sampling

% Create vector of time range in sampling
t_range = (start_time:1/Fs:end_time);

% Load the signal file
signal_file = load('signal476.mat');
signal_data = signal_file.xn_test;

% Assiging signal values to subsets
S1 = signal_data(1:128);
S2 = signal_data(1:256);
S3 = signal_data(1:512);
S4 = signal_data(1:1024);
S5 = signal_data(1:1792);

signal_subsets = {S1, S2, S3, S4, S5};
dft_magnitudes = {};

% Loop through ech subset to compute the DFT magnitudes
for i = 1:5
    % Compute the DFT of the subset
    dft_transform = fft(signal_subsets{i});
    
    % Calculate the magnitude of the DFT coefficients
    magnitude = abs(dft_transform);
    
    dft_magnitudes{i} = magnitude;
    
    % Plot the magnitude for the current subset
    subplot(length(signal_subsets), 1, i);
    frequencies = (0:length(signal_subsets{i})-1) * Fs / length(signal_data);
    stem(frequencies, magnitude);
    title(['S', num2str(i)]);

end

% Label the x-axis for the last subplot
xlabel('Frequency (Hz)');





%% 3.14 Apply DFT Averaging method

% Initialize the variables/ parameters
total_samples = length(signal_data);
L = 14;
K = floor(total_samples / L);
i = 1; j = K;
mean_dft = 0;

for n = 1:L
    subset = signal_data(i:j);
    dft_magnitude = abs(fft(subset, K));
    mean_dft = mean_dft + dft_magnitude;
    
    i = j+1;
    j = j+K;
end

mean_dft = mean_dft / L;

% Plot the magnitude for the current subset
subplot(1, 2, 1);
plot(mean_dft);
title('Mean of DFT Sequences (Line Plot)');
xlabel('Frequency (Hz)');
grid on;
subplot(1, 2, 2);
stem(mean_dft);
title('Mean of DFT Sequences (Stem Plot)');
xlabel('Frequency (Hz)');
grid on;

% The least value of L which gives 4 distinguishable harmonics which are
% less than 64Hz can be considered as 12 or 13.

% Implement DFT averaging for different L values for the same subset length
for L = 1:14
    x = 1;
    y = 128;
    dft_mean_new = 0;
    
    for n = 1:L
        subset = signal_data(x:y);
        dft_magnitude = abs(fft(subset, 128));
        dft_mean_new = dft_mean_new + dft_magnitude;
    
        x = y+1;
        y = y+128;
    end
    
    dft_mean_new = dft_mean_new / L;
    
    % Plot the average DFTs for different L values with the same subset
    % length.
    figure;
    subplot(1, 1, 1);
    stem(dft_mean_new);
    title(['Mean of DFT Sequences for L=',num2str(L)]);
    xlabel('Frequency (Hz)');
end


% Implementing DFT averaging method for different K values.
K = [100, 135, 512, 256];

for i = 1:length(K)
    x = 1;
    y = K(i);
    dft_mean = 0;
    subset_count = floor(length(signal_data)/K(i));
    
    for j = 1:subset_count
        subset = signal_data(x:y);
        dft_magnitude = abs(fft(subset, K(i)));
        dft_mean = dft_mean + dft_magnitude;
    
        x = y+1;
        y = y+K(i);
    end
    
    dft_mean = dft_mean / subset_count;
    
    % Plot the average DFTs for different L values with the same subset
    % length.
    figure;
    subplot(1, 1, 1);
    stem(dft_mean);
    title(['Mean of DFT Sequences for K=',num2str(K(i))]);
    xlabel('Frequency (Hz)');
end




%% 3.2 Interpolation
clc;
clear;
close all;

% Load the signal and store it in variable "y"
load handel;

% Get first 20,000 sample values
N = 20000;                      % Number of samples needed
signal_subset = y(1:N);

% Define various signals based on the obtained samples
x = y(1:N);
x2 = x(1:2:N);
x3 = x(1:3:N);
x4 = x(1:4:N);

% Plot original signal
figure;
subplot(1, 1, 1);
stem(x(1:50), 'b', 'Marker', 'o');
title('X Original Signal (First 50 Samples)');

% Signal interpolation for x2 by adding zeros into the middle of the signal sample
% in ferquency domain.
K1 = 1;
X2 = fft(x2);
N1 = length(X2);
if mod(N1, 2) == 0
    N = N1/2;
    zero_padded_X2 = [X2(1:N); X2(N+1)/2; zeros((K1*N1)-1, 1); X2(N+1)/2; X2((N+2):N1)];
else
    N = (N1+1)/2;
    zero_padded_X2 = [X2(1:N); zeros(K1*N1, 1); X2((N+1):N1)];
end

interpolated_signal_x2 = (K1+1)*ifft(zero_padded_X2);
edited_x2 = interpolated_signal_x2(1:((K1+1)*(N1-1))+2);
difference_2 = norm(edited_x2 - x);

% Plot the first 50 samples of both X2 and interpolated version of X2
% using stem plots
figure;
subplot(3, 2, 1);
stem(x2(1:50), 'b', 'Marker', 'o');
title('X2 Original Signal');
subplot(3, 2, 2)
stem(edited_x2(1:50), 'r', 'Marker', 'x');
title('Interpolated Version of X2');


% Signal interpolation for x3 by adding zeros into the middle of the signal sample
% in ferquency domain.
K2 = 2;
X3 = fft(x3);
N2 = length(X3);
if mod(N2, 2) == 0
    N = N2/2;
    zero_padded_X3 = [X3(1:N); X3(N+1)/2; zeros((K2*N2)-1, 1); X3(N+1)/2; X3((N+2):N2)];
else
    N = (N2+1)/2;
    zero_padded_X3 = [X3(1:N); zeros(K2*N2, 1); X3((N+1):N2)];
end

interpolated_signal_x3 = (K2+1)*ifft(zero_padded_X3);
edited_x3 = interpolated_signal_x3(1:((K2+1)*(N2-1))+2);
difference_3 = norm(edited_x3 - x);

% Plot the first 50 samples of both original and interpolated versions
% using stem plots
subplot(3, 2, 3);
stem(x3(1:50), 'b', 'Marker', 'o');
title('X3 Original Signal');
subplot(3, 2, 4)
stem(edited_x3(1:50), 'r', 'Marker', 'x');
title('Interpolated Version of X3');


% Signal interpolation for x4 by adding zeros into the middle of the signal sample
% in ferquency domain.
K3 = 3;
X4 = fft(x4);
N3 = length(X4);
if mod(N3, 2) == 0
    N = N3/2;
    zero_padded_X4 = [X4(1:N); X4(N+1)/2; zeros((K3*N3)-1, 1); X4(N+1)/2; X4((N+2):N3)];
else
    N = (N3+1)/2;
    zero_padded_X4 = [X4(1:N); zeros(K3*N3, 1); X4((N+1):N3)];
end

interpolated_signal_x4 = (K3+1)*ifft(zero_padded_X4);
edited_x4 = interpolated_signal_x4(1:((K3+1)*(N3-1))+4);
difference_4 = norm(edited_x4 - x);

% Plot the first 50 samples of both original and interpolated versions
% using stem plots
subplot(3, 2, 5);
stem(x4(1:50), 'b', 'Marker', 'o');
title('X4 Original Signal');
subplot(3, 2, 6)
stem(edited_x4(1:50), 'r', 'Marker', 'x');
title('Interpolated Version of X4');

% Print the 2-norm difference values between original signal and
% interpolated signal.
fprintf('The difference between X2 and the original signal x in 2-norm: %f\n ', difference_2);
fprintf('The difference between X3 and the original signal x in 2-norm: %f\n ', difference_3);
fprintf('The difference between X4 and the original signal x in 2-norm: %f\n ', difference_4);