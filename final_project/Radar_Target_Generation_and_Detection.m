% Final Project - Radar Target Generation and Detection
% Carlos Gómez Huélamo

clc, clear all, close all;
%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity remains contant

R_target = 110; % Initial distance of the target (m)
v_target = -20; % Speed of the target (approaching vehicle) (m/s)

% Maximum range can be 200 m
% Speed can range between -70 m/s and 70 m/s

%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

%Operating carrier frequency of Radar 
fc = 77e9; % Frequency of operation (Hz)
max_range = 200; % (m)
range_resolution = 1; % (m)
max_velocity = 70; % (m/s)
c = 3e8; % Speed of light (m/s)

plotting = 1;

% Bandwidth for each chirp for given resolution

bandwidth = c / (2 * range_resolution);

% The sweep time can be computed based on the time needed for the signal to travel the unambiguous
% maximum range. In general, for an FMCW radar system, the sweep time should be at least
% 5 to 6 times the round trip time. This example uses a factor of 5.5.

t_chirp = 5.5 * 2 * max_range / c;

slope = bandwidth / t_chirp;

disp('slope')
disp(slope)
                                                   
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 

Nd = 128; % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp.

Nr = 1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp

t = linspace(0,Nd*t_chirp,Nr*Nd); % total time for samples

disp('Start')
disp(0)
disp('End')
disp(Nd*t_chirp)
disp('Length t')
disp(length(t))
disp('Step')
disp((Nd*t_chirp)/length(t))

% Creating the vectors for Tx, Rx and Mix based on the total samples input.

Tx = zeros(1,length(t)); % transmitted signal
Rx = zeros(1,length(t)); % received signal
Mix = zeros(1,length(t)); % beat signal (Substracted signal)

% Similar vectors for range_covered and time delay.

r_t = zeros(1,length(t));
td = zeros(1,length(t));

%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i = 1:length(t)             
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    
    r_t(i) = R_target + v_target*t(i); % r_t = vector of range covered
    
    % *%TODO* :
    %For each time sample we need update the transmitted and received signal.  
    
    % The receiving signal is the time delayed version of the transmit
    % signal. That is, t is substituted by t-tau, where tau represent the
    % delay time
    
    tau = (r_t(i) * 2) / c;  
    t_aux = t(i) - tau;
    
    Tx(i) = cos(2*pi*(fc*t(i) + (1/2)*slope*(t(i)^2))); % Equation of transmitting signal
    Rx(i)  = cos(2*pi*(fc*t_aux + (1/2)*slope*(t_aux^2)));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    
    Mix(i) = Tx(i).*Rx(i); % Element wise multiplication
end

%% RANGE MEASUREMENT

 % *%TODO* :
% reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
% Range and Doppler FFT respectively.

Mix = reshape(Mix,[Nr,Nd]);

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.

signal_fft = fft(Mix,Nr); % Second argument: Y = fft(X,n) returns the n-point DFT. 
% If no value is specified, Y is the same size as X.

 % *%TODO* :
% Normalize by the maximum y-axis value, and then take the absolute value of 
% FFT output

signal_fft = signal_fft ./ max(signal_fft); % Element wise operation
signal_fft = abs(signal_fft);

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
% We are only interested in the positive half of signal length L

L = length(signal_fft);
signal_fft = signal_fft(1:L/2+1); % from A to B-1, to we sum + 1

%plotting the range
%figure ('Name','Range from First FFT')

 % *%TODO* :
 % plot FFT output 

if (plotting)
    figure
    plot(signal_fft)
    title('Range from First FFT')
    xlabel('range (m)')
    ylabel('Normalized Signal Amplitude')
    axis ([0 200 0 1]);
end


%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map. You will implement CFAR on the generated RDM (Range Doppler Map)

% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix = reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.

sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.

sig_fft2 = sig_fft2(1:Nr/2, 1:Nd); % Middle of rows
sig_fft2 = fftshift (sig_fft2);
Range_Doppler_Map = abs(sig_fft2);
Range_Doppler_Map = 10*log10(Range_Doppler_Map) ;

% Use the surf function to plot the output of 2DFFT and to show axis in both
% dimensions

if (plotting)
    doppler_axis = linspace(-100, 100, Nd);
    range_axis = linspace(-200, 200, Nr/2) * ((Nr/2)/400);
    figure
    surf(doppler_axis, range_axis, Range_Doppler_Map);
    colorbar;
    title('Results after FFT2 to analyze Doppler Effect');
    xlabel('Speed (Doppler axis), m/s')
    ylabel('Range (Range axis), m')
    zlabel('Normalized Signal Amplitude')
end

%% CFAR implementation

% Slide Window through the complete Range Doppler Map

% 1. Determine the number of Training cells for each dimension Tr and Td. 
%    Similarly, pick the number of guard cells Gr and Gd.

% *%TODO* :
%Select the number of Training Cells in both the dimensions.

Tr = 10; % Professor's
Td = 8; % Professor's

% *%TODO* :
% Select the number of Guard Cells in both dimensions around the Cell under 
% test (CUT) for accurate estimation

Gr = 4; % Professor's
Gd = 4; % Professor's

% *%TODO* :
% offset the threshold by SNR value in dB

offset = 6.0; % dB

% 2. Slide the Cell Under Test (CUT) across the complete cell matrix.

% *%TODO* :
% design a loop such that it slides the CUT across range doppler map by
% giving margins at the edges for Training and Guard Cells.

% For every iteration sum the signal level within all the training
% cells. To sum convert the value from logarithmic to linear using db2pow
% function. Average the summed values for all of the training
% cells used. After averaging convert it back to logarithimic using pow2db.

% 3. Select the grid that includes the training, guard and test cells. Grid Size = (2Tr+2Gr+1)·(2Td+2Gd+1).
%    Note that 1 means the size of CUT in each dimension.
% 4. The total number of cells in the guard region and cell under test. (2Gr+1)(2Gd+1).
% 5. This gives the Training Cells : (2Tr+2Gr+1)(2Td+2Gd+1) - (2Gr+1)(2Gd+1).

grid_size = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1);
guard_CUT_size = (2*Gr+1)*(2*Gd+1);
number_training_cells = grid_size - guard_CUT_size;

for p = Tr+Gr+1 : Nr/2 - (Tr+Gr) % rows
    for q = Td+Gd+1 : Nd - (Td+Gd) % columns 
        
        % At this moment, RDM(p,q) represents our CUT (Cell Under Test)
        
        % Iterative through all training cells relative to this CUT
        
        % 6. Measure and average the noise across all the training cells. This gives the threshold.
        
        % *%TODO* :
        %Create a vector to store noise_level for each iteration on training cells

        accumulated_noise_level = zeros(1,1); % Single value (e.g. average of noise)
        
        for a = p - (Tr+Gr) : p + (Tr+Gr)
            for b = q - (Td+Gd) : q + (Td+Gd)
                if (abs(p-a) > Gr || abs(q-b) > Gd) % If training cell
                    accumulated_noise_level = accumulated_noise_level + db2pow(Range_Doppler_Map(a,b));
                end
            end
        end
        
        %disp('cnt')
        %disp(cnt)
        
        % Calculate threshold from average noise 
        
        average_noise = accumulated_noise_level / number_training_cells;
        threshold = average_noise;
        
        % 7. Add the offset (if signal strength is in dB) to the threshold to keep the false alarm to the minimum.
        %    If linear units -> Multiply. If dB (logarithmic) -> sum.
        
        % Further add the offset to it to determine the threshold. Next, compare the
        % signal under CUT with this threshold. If the CUT level > threshold assign
        % it a value of 1, else equate it to 0.
        
        % Retrieve dB
        
        threshold = pow2db(threshold);
        
        % Add offset
        
        threshold = threshold + offset;
        
        CUT = Range_Doppler_Map(p,q); % Value of the signal in dB in the CUT
        
%         disp('CUT')
%         disp(CUT)
%         disp('Threshold')
%         disp(threshold)

        if (CUT < threshold)
            Range_Doppler_Map(p,q) = 0;
        else
            Range_Doppler_Map(p,q) = 1;
        end
    end
end

% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
% than the Range Doppler Map as the CUT cannot be located at the edges of
% matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 

Range_Doppler_Map(Range_Doppler_Map ~= 0 & Range_Doppler_Map ~= 1) = 0; % Evaluate condition over the whole matrix (if the element
% is different from 0 and different from 1 (so, no evalauted) -> Assign 0)
 
% *%TODO* :
% display the CFAR output using the Surf function like we did for Range
% Doppler Response output.

if (plotting)
    figure
    surf(doppler_axis, range_axis, Range_Doppler_Map);
    colorbar;
    title('CA-CFAR Filtered RDM Surface plot');
    xlabel('Speed (Doppler axis), m/s')
    ylabel('Range (Range axis), m')
    zlabel('Normalized Signal Amplitude')
end
 