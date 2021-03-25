function [rhythm, phase, mag]= filterlfp(t, lfp, w, wch)
% FILTERLFP filters out theta from the raw LFP signal
% 
% (I)   filterlfp(time, signal, w)
%
% Two methods:
%
% (1) w = {'delta'; 'theta'; 'beta'; 'gamma'}
% (2) w = [wcl wch] e.g [6 9] for theta
%       wcl : lower cutoff frequency
%       wch : higher cutoff frequency
%
% (II)  filterlfp(time, signal, wcl, wch)
%

if nargin < 2
    [t, lfp] = readcsc;
end

Ts = mean(diff(t));
SamplingFrequency = round(1/Ts);

if nargin < 3
    w = 'theta';
elseif nargin == 4
    w = [w wch];
end

if isstring(w) || ischar(w)
    switch w
        case 'delta'
            w = [0.5 3.5];
        case 'theta'
            w = [6 9];
        case 'beta'
            w = [10 20];
        case 'gamma'
            w = [30 50];
        otherwise
            w = [6 9];
            disp('not recognized, but theta is chosen.')
    end
end

wn = w/(SamplingFrequency/2);
order = 2;

[b,a] = butter(order, wn);
rhythm = filtfilt(b,a,double(lfp));

if nargout == 0
    close all;
    plot(t, lfp, 'Color', uint8([230 230 230]));
    hold on
    plot(t, rhythm,'r');
    clear rhythm;
elseif nargout > 1
    z = hilbert(rhythm);
    phase = angle(z);
    phase = rad2deg(phase);
    mag = abs(z);
end    
    
end