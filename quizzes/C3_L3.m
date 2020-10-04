% Class 3. Range-Doppler Estimation - Lesson 3. Doppler Estimation
% Carlos Gómez Huélamo

% complete the TODOs to calculate the velocity in m/s of four targets 
% with following doppler frequency shifts: [3 KHz, -4.5 KHz, 11 KHz, -3 KHz]

% Radar operating frequency (Hz)

f = 77e9; 

% Speed of light (m/s)

c = 3e8;

% TODO : Calculate the wavelength

lambda = c/f;

% TODO : Define the doppler shifts in Hz using the information from above 

fd = [3e3, -4.5e3, 11e3, -3e3];

% TODO : Calculate the velocity of the targets  fd = 2*vr/lambda

vr = (fd * lambda) / 2;

% TODO: Display results

disp(vr)

% Question 1/2: From the workspace above, estimate the location of the targets. 
% If all the targets were at 200 m range ahead of the ego (radar) vehicle and the velocity 
% of the ego vehicle is 5 m/s. Then in next 5 seconds which of the targets would be closest 
% to the ego vehicle

