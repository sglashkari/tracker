function [theta, phase]= filtertheta(t, lfp, wcl, wch)
% FILTERTHETA filters out theta from the raw LFP signal
% filtertheta(time, signal, lower cutoff frequency, higher cutoff frequency)

if nargin < 2
    [t, lfp] = readcsc;
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
theta = filtfilt(b,a,lfp);

if nargout == 0
    close all;
    plot(t, lfp, 'Color', uint8([230 230 230]));
    hold on
    plot(t, theta,'r');
    clear theta;
elseif nargout > 1
    z = hilbert(theta);
    phase = angle(-z); % negative because trough is zero
    phase = rad2deg(phase);
end    
    
end