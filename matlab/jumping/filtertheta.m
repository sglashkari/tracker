function [theta, phase, mag]= filtertheta(t, lfp, wcl, wch)
% FILTERTHETA filters out theta from the raw LFP signal
% filtertheta(time, signal, lower cutoff frequency, higher cutoff frequency)

if nargin < 2
    [t, lfp] = read_bin_csc;
end

Ts = mean(diff(t));
SamplingFrequency = 1/Ts;
if nargin < 3
    wcl = 5; % lower cutoff frequency of 5 Hz
end
if nargin < 4
    wch = 12; % higher cutoff frequency of 12 Hz
end

wn = [wcl wch]/(SamplingFrequency/2);
order = 2;

[b,a] = butter(order, wn);
theta = filtfilt(b,a,double(lfp));

if nargout == 0
    close all;
    plot(t, lfp, 'Color', uint8([230 230 230]));
    hold on
    plot(t, theta,'r');
    clear theta;
elseif nargout > 1
    z = hilbert(theta);
    phase = angle(z);
    phase = rad2deg(phase);
    mag = abs(z);
end    
    
end