% Class 3. Range-Doppler Estimation - Lesson 4. Fast Fourier Transform (FFT)
% Carlos Gómez Huélamo

% Use a Fourier transform to find the frequency components of a signal buried in noise. 
% Specify the parameters of a signal with a sampling frequency of 1 kHz and a signal duration 
% of 1.5 seconds:

% Set parameters

Fs = 1000; % Sampling frequency (Hz)
T = 1/Fs; % Sampling period (s)
L = 1500; % Length of signal (ms)
t = (0:L-1)*T; % Time vector from 0 to 1.5 using T as step (each step is 1 ms), so t is in seconds

A1 = 0.7;
f1 = 77;

A2 = 2;
f2 = 43;

% Define a signal containing a 77 Hz sinusoid of amplitude 0.7 and a 43Hz sinusoid of amplitude 2
% So, this signal contains two sinusoidal waves!!!!

signal = A1*sin(2*pi*f1*t) + A2*sin(2*pi*f2*t);

% Corrupt the signal with noise (white noise, thermal noise)

X = signal + 2*randn(size(t));

% Plot the noisy signal in the time domain. 
% It is difficult to identify the frequency components by looking at the signal X(t)

% plot(1000*t(1:50), X(1:50)) % t is in seconds, if we multiply by 1000, it is miliseconds
% title('Signal Corrupted with Zero-Mean Random Noise')
% xlabel('t (miliseconds)')
% ylabel('X(t)')

% Run the FFT for dimension of samples N

%signal_fft = fft(X,L);
signal_fft = fft(X);

% The output of FFT processing of a signal is a complex number (a+jb). 
% Since, we just care about the magnitude we take the absolute value 
% (sqrt(a^2+b^2)) of the complex number

signal_fft = abs(signal_fft/L);
%signal_fft = abs(signal_fft);
%signal_fft = signal_fft ./ max(signal_fft);

% FFT output generates a mirror image of the signal. 
% But we are only interested in the positive half of signal length L, 
% since it is the replica of negative half and has all the information we need

% In other words: We just compute the single-sided spectrum as we reject
% the mirror image

signal_fft = signal_fft(1:L/2+1);

% Plotting

f = Fs*(0:(L/2))/L;
plot(f,signal_fft)
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1 (f)|')

