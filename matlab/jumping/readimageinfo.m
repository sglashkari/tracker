function [Timestamp, GPIO] = readimageinfo(filename, N)
%READIMAGEINFO takes an image filename and returns its embedded information
%   [Timestamp, GPIO] = readimageinfo(filename, N) returns the timestamp and
%   gpio for N frames of the filename
%   Timestamp:          N x 1
%   GPIO(pin 0 to 3):   N x 4 
% For more information check this out
% <a href="matlab:web('https://www.flir.com/globalassets/support/iis/knowledge-base/flir-machine-vision-camera-register-reference.pdf#page=88','-browser')">Embedded Image Information</a>

if nargin < 1
    filename = 'C:\Users\Shahin\test\matlab\automation\fc2_save_2020-10-22-222431-';
    N = 10;
elseif nargin <2
    N = 1;
end

Timestamp = zeros(N,1);
GPIO = zeros(N,4);
for i = 1:N
    A1 = imread([filename num2str(i-1,'%04.f') '.pgm'],'pgm');
    
    % Timestamp
    B1 = A1(1,1:4);
    hexStr = dec2bin(B1,8);
    Z = hexStr';
    Z = Z(:)';
    Timestamp(i) = bin2dec(Z(1:7))+125e-6*bin2dec(Z(8:20)); % [0,128)
    
    % GPIO
    B2 = A1(1,5);
    hexStr = dec2bin(B2,8);
    GPIO(i,1:4) = hexStr(1:4)-'0000';
end
Timestamp = unwrap((Timestamp-64)/64*pi)/pi*64+64; %unwrap the timestamp