function [ pattimg ] = PatternGenerateSimple_SEP( xcent, ycent, be, al, field, SetupParams, SEPdata)
%PATTERNGENERATESIMPLE Summary of this function goes here
%   Detailed explanation goes here

pattimg = PatternGeneratePos_SEP(xcent,ycent,0, SetupParams.NA, [], [], [], [], [], [], SetupParams.lamem, SetupParams.mag, SetupParams.focus, [], [], SetupParams.pixelsize, SetupParams.nn, field, be, al, [], SEPdata);
pattimg = pattimg/sum(pattimg(:));


end

