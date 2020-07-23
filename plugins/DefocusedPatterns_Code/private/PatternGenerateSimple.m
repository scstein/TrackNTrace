function [ pattimg ] = PatternGenerateSimple( xcent, ycent, be, al, NA, lamem, mag, focus, pixelsize, nn, field)
%PATTERNGENERATESIMPLE Summary of this function goes here
%   Detailed explanation goes here

% NA          - Objective NA
% lamem       - emission wavelength [µm]
% mag = 100   - magnification
% focus = 0.9 -(de-)focus position
% pixelsize   - camera pixel size [µm];
% nn          - radius of disk pattern is calculated for
% field       - size of generated image is 2*field+1
% be          - in plane angle
% al          - out of plane angle

z = 0; % always 0
pattimg = PatternGeneratePos(xcent,ycent,z, NA, [], [], [], [], [], [], lamem, mag, focus, [], [], pixelsize, nn, field, be, al, []);
pattimg = pattimg/sum(pattimg(:));

end

