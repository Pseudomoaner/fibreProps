function [] = plotLocalOrientationHistogram(axes,cMap,normalise,inDat,labels)
%PLOTLOCALORIENTATIONHISTOGRAM plots the histogram of the local
%orientations of fibres detected in the input AFM data.
%
%   INPUTS:
%       -axes: The target axes you want to plot into
%       -cMap: The colourmap you want to use to distinguish different
%       samples
%       -normalise: Whether to use raw orientation bin counts or normalised
%       (PDF) values.
%       -inDat: The data for plotting. In a cell array, with each cell
%       containing a different vector of angles.
%       -labels: Strings defining the labels by which each dataset will be
%       referred to in the legend.
%
%   Author: Oliver J. Meacock, (c) 2021

noBins = 18;

angSet = linspace(0,pi,noBins+1);
angList = (angSet(1:end-1) + diff(angSet(1:2))/2)';
angList = [angList;angList + pi;angList(1)];

for i = 1:size(inDat,1)
    if normalise
        N = histcounts(inDat{i}+pi/2,'BinEdges',angSet,'Normalization','pdf');
    else
        N = histcounts(inDat{i}+pi/2,'BinEdges',angSet,'Normalization','count');
    end

    %Repeat N so the plot goes the full 360 degrees.
    N = [N';N';N(1)];
    
    polarplot(axes,angList,N,'Color',cMap(i,:),'lineWidth',1.5)
end

legend(axes,labels,'Location','northwest','Interpreter','none')

axes.ThetaLim = [0,180];