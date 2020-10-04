% Class 4. Clutter, CFAR, AoA - Lesson 4. CA-CFAR
% Carlos Gómez Huélamo

% The following steps here can be used to implement CFAR in the next MATLAB exercise. 
% You can use the code template below to get started as well

% 1. Define the number of training cells and guard cells
% 2. Start sliding the window one cell at a time across the complete FFT 1D array.
%    Total window size should be: 2(T+G) + CUT
% 3. For each step, sum the signal (noise) within all the leading or lagging training cells
% 4. Average the sum to determine the noise threshold
% 5. Using an appropriate offset value scale the threshold
% 6. Now, measure the signal in the CUT, which is T+G+1 from the window starting point
% 7. Compare the signal measured in 5 against the threshold measured in 4
% 8. If the level of signal measured in CUT is smaller than the threshold measured, 
%    then assign 0 value to the signal within CUT

% Implement 1D CFAR using lagging cells on the given noise and target scenario.

% Close and delete all currently open figures
clc, clear, close all;

% Data_points
Ns = 1000;

% Generate random noise
s = abs(randn(Ns,1)); % 0 mean noise. Standard deviation = 1

%Targets location. Assigning bin 100, 200, 300 and 700 as Targets with the amplitudes of 8, 9, 4, 11.
s([100 ,200, 300, 700])=[8 9 4 11]; % At this points, replaces noise values with these amplitudes

%plot the output
plot(s);

% TODO: Apply CFAR to detect the targets by filtering the noise.

% 1. First define a CFAR window and pick the optimal number of Training 
%    and Guard Cells:
% 1a. Training Cells

T = 12;

% 1b. Guard Cells 

G = 4;

% Offset : Adding room above noise threshold for desired SNR 
% Here we are working on linear values, hence we multiply the offset to the
% threshold values
offset=3;

% Vector to hold threshold values 
threshold_cfar = [];

%Vector to hold final signal after thresholding
signal_cfar = [];

% Steo through all the cells from one end to another of the signal vector
% but make sure that you have right spacing from the end
% 2. Slide window across the signal length
for i = 1:(Ns-(G+T+1)) % Note that the last position we evaluate (so, Cell Under Test)
    % is that which associated last leading training cell reaches the end
    % of the matrix
    
    disp(i)

    % 2. - 5. Determine the average noise threshold by measuring it within the training cells
    
    noise_level = sum(s(i:i+T-1)) / T; % Noise is based on lagging training cells
    
    % To determine the threshold take the average of summed noise and
    % multiply it with the offset
    
    threshold = noise_level * offset;
    threshold_cfar = [threshold_cfar, {threshold}]; % Add a new element to the vector

    % 6. Measuring the signal within the CUT
    
    % Now pick the cell under test which is T + G cells away from the first
    % training cell and measure the signal level
    
    signal = s(i+T+G); % Signal value within the CUT

    % 7. 8. If the signal level at CUT is below the threshold, assign it a
    % 0 value
    
    if (signal < threshold)
        signal = 0;
    end

    signal_cfar = [signal_cfar, {signal}];
end

% plot the filtered signal
plot (cell2mat(signal_cfar),'g--');

% plot original sig, threshold and filtered signal within the same figure.
figure,plot(s);
hold on,plot(cell2mat(circshift(threshold_cfar,G)),'r--','LineWidth',2)
hold on, plot (cell2mat(circshift(signal_cfar,(T+G))),'g--','LineWidth',4);
legend('Signal','CFAR Threshold','detection')

