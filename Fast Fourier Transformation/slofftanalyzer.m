% this function demands a set of (non double) reflectivity and aequidistant
% delayvalues ('reflectivity' and 'delay'), a lower index 'k' and an upper
% index 'l', cuts the data in that index range, produces a smoothed one with
% a given 'smoothrange', subtracts the unsmoothed with the smoothed one,
% applies a window function and outputs the FFT of this ('fftcut' over
% 'fftcut'). Also for the FFT provide information for the pulsewidth
% ('pulsew') remark: Window-function 'window_1.m' necessary. parameter
% 'optionall' = 1, cutdelay and rdy2windowref is given for the whole spectrum
% (incl the area before zero)
%
%
% created by n1oftheabove, 2013 at Potsdam University. Use as you like.

function [ fftcutfreq, fftcut, cutdelay, rdy2windowref ] = slofftanalyzer(delay, reflectivity, k, l, smoothrange, pulsew, optionall)
cutdelay = delay(k:l);
cutreflectivity =  reflectivity(k:l);

smoothedreflectivity = smooth(cutreflectivity,smoothrange);
rdy2windowref = cutreflectivity - smoothedreflectivity;

win_width = 0.99*abs(max(cutdelay)-min(cutdelay));
win_cent = mean(cutdelay);
win_edge = 0.02;
window = window_1(win_width, win_edge, win_cent,cutdelay);
windwd_corcut = window.*rdy2windowref;

if max(delay).*min(delay) < 0
    delayrange = abs(max(delay))+abs(min(delay));
else
    delayrange = abs(max(delay))-abs(min(delay));
end
samples= length(delay);

T = delayrange./samples;                                        % time intervall between two measured points
Fs = 1/T;                                                       % Sampling Frequency in THz
      
NFFT = length(cutdelay);
Y = fft(windwd_corcut, NFFT)/NFFT;
f = Fs/2 * linspace(0,1,floor(NFFT/2)+1);
Ampdens = 2*abs(Y(1:floor(NFFT/2)+1))';

% first index in frequency vector, from which on frequency is bigger than inverse pulse width

fmax = find(f >(1000*1./pulsew),1,'first');
fftcutfreq= f(1:fmax);
fftcut= Ampdens(1:fmax);

if optionall == 1
    cutrefall = zeros(l,1);
    cutrefall(1:(k-1),1) = reflectivity(1:(k-1)) - smooth(reflectivity(1:(k-1)),5);
    cutrefall(k:l,1) = cutreflectivity - smoothedreflectivity;
    rdy2windowref = cutrefall;
    cutdelay = delay(1:l);
end
