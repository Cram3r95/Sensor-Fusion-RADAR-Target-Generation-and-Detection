% Class 2. Radar Principles - Lesson 7. Radar Range Equation
% Carlos Gómez Huélamo

%Operating frequency (Hz)
fc = 77.0e9;

%Transmitted power (W)
Pt = 3e-3;

%Antenna Gain (linear)
G =  10000;

%Minimum Detectable Power (dBm)
Ps = 1e-10;

%RCS of a car (m^2)
RCS = 100;

%Speed of light (m/s)
c = 3*10^8;

%TODO: Calculate the wavelength

% c = f*lambda (operating frequency * associated_wavelength)

lambda = c/fc;

%TODO : Measure the Maximum Range a Radar can see

% Y = nthroot(X,N) returns the real nth root of the elements of X. Both X and N must be real scalars or arrays of the same size.
% If an element in X is negative, then the corresponding element in N must be an odd integer

aux = (Pt * G^2 * lambda^2 * RCS) / (Ps * (4*pi)^3);
Range = nthroot(aux,4)