function [rhythm, phase, mag]= thetafilter(t, lfp)
% FILTERLFP filters out theta from the raw LFP signal
% By Yotaro (modifed 2023-02-14)
%
Ts = mean(diff(t));
Fs = round(1/Ts);
Fstop1 = 3;     % First Stopband Frequency
Fpass1 = 6;     % First Passband Frequency
Fpass2 = 12;    % Second Passband Frequency
Fstop2 = 15;    % Second Stopband Frequency
Astop1 = 60;    % First Stopband Attenuation (dB)
Apass  = 1;     % Passband Ripple (dB)
Astop2 = 60;    % Second Stopband Attenuation (dB)

h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1, ...
    Fpass2, Fstop2, Astop1, Apass, Astop2, Fs);

Hd = design(h, 'butter', ...
    'MatchExactly', 'stopband', ...
    'SOSScaleNorm', 'Linf');

% Extract theta phase by taking hilbert transform of filtered lfp
rhythm = filtfilt(Hd.sosMatrix,Hd.ScaleValues,double(lfp));

rhythm2= filterlfp(t, lfp);

if nargout == 0
    close all;
    plot(t, lfp, 'Color', uint8([200 200 200]));
    hold on
    plot(t, rhythm,'r');
    plot(t, rhythm2,'b');
    clear rhythm;
end
if nargout > 1
    z = hilbert(rhythm);
    phase = rad2deg(angle(z)); 
end
if nargout > 2
    mag = abs(z);
end
    
end