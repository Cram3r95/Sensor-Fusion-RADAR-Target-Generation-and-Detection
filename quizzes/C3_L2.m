% Class 3. Range-Doppler Estimation - Lesson 2. Range Estimation
% Carlos Gómez Huélamo

% Calculate the range of four targets with respective measured beat
% frequencues (0 MHz, 1.1 MHz, 13 MHz, 24 MHz)

% Radar maximum range (m)

Rmax = 300;

% Range resolution (m)

range_resolution = 1;

% Speed of the light (m/s)

c = 3*10^8;

% TODO: Find the Bandwidth sweep of chirp for 1 m resolution
% Bandwidth = (c*range_resolution (aka delta_r)) / 2

bandwidth = c / (2 * range_resolution);

% TODO: Calculate the chirp time based on the Radar Max Range

% Note : The sweep time can be computed based on the time needed for the signal to travel 
% the maximum range. In general, for an FMCW radar system, the sweep time should be at least 
% 5 to 6 times the round trip time. This example uses a factor of 5.5:

t_chirp = (5.5 * 2 * Rmax) / c; % t_chirp is the time the electromagnetic wave takes
% to cover the whole bandwidth (from fA to fB, ej: from 77 MHz to 77.4 MHz)

% TODO define the frequency shifts

beat_freq = [0 1.1e6 13e6 24e6];
calculated_range = (c * t_chirp * beat_freq) / (2 * bandwidth);

% Display results

disp(calculated_range)





