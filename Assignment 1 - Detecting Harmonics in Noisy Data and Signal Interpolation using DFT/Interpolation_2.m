%% 3.2 Interpolation

clc;
clear;
close all;

% Load the signal and store it in variable "y"
load handel;

% Get first 20,000 sample values
N = 2000;                      % Number of samples needed
signal_subset = y(1:N);

% Define various signals based on the obtained samples
x = y(1:N);
x2 = x(1:2:N);
x3 = x(1:3:N);
x4 = x(1:4:N);

% 3.1 Signal interpolation by adding zeros between each two samples

K = [1, 2, 3];

% Declaring arrays for original and interpolated signals
interpolated_signals = {[], [], []};
original_signals = {x2, x3, x4};
difference_signals = {[], [], []};

for j = 1:length(original_signals)
    for i = 1:length(original_signals{j})
        x = original_signals{j};
        interpolated_signals{j} = [interpolated_signals{j}, x(i)];

        for m = 1:K(j)
            interpolated_signals{j} = [interpolated_signals{j}, 0];
        end
    end
    
   % Ensure both signals have the same lengths (zero-pad if needed)
   y = original_signals{j};
   y(end+1: length(interpolated_signals{j})) = 0;
   
   % Calculate the difference between original and interpolated version
   difference_signals{j} = interpolated_signals{j} - original_signals{j};
   
   % Plot the first 50 samples of both original and interpolated versions
   % using stem plots
   figure;
   subplot(2, 1, 1);
   stem(original_signals{j}(1:50), 'b', 'Marker', 'o');
   title(['X', num2str(j+1), ' Original Signal']);
    
   hold on;
   subplot(2, 1, 2)
   stem(interpolated_signals{j}(1:50), 'r', 'Marker', 'x');
   title(['Interpolated Version of X',num2str(j+1)]);

   % Add labels and a legend
   xlabel('Time (samples)');
   ylabel('Amplitude');

   % Hold off to stop overlaying plots
   hold off;
end


% % Calculate the difference between signals
% difference_signal_2 = interpolated_signal_2 - x2;
% difference_signal_3 = interpolated_signal_3 - x3;
% difference_signal_4 = interpolated_signal_4 - x4;

