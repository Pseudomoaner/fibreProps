function [] = plotLocalOrientationHistogram(axes,inDat,cMap)
%PLOTLOCALORIENTATIONHISTOGRAM plots the histogram of the local
%orientations of fibres detected in the input AFM data.
%
%   INPUTS:
%       -fibreProps: A collection of measurements about the fibre network,
%       including the localOrientation field.
%
%   Author: Oliver J. Meacock, (c) 2021

noBins = 18;

angSet = linspace(0,pi,noBins+1);
angList = (angSet(1:end-1) + diff(angSet(1:2))/2)';
angList = [angList;angList + pi;angList(1)];

for i = 1:size(inDat,1)
    N = histcounts(inDat{i}+pi/2,'BinEdges',angSet);

    %Repeat N so the plot goes the full 360 degrees.
    N = [N';N';N(1)];
    
    polarplot(axes,angList,N,'Color',cMap(i,:),'lineWidth',1.5)
end

axes.ThetaLim = [0,180];