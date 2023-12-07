function [] = plotFibreHistogram(axes,cMap,normalise,data,labels)
%PLOTFIBREHISTOGRAM plots the histogram of specified dataset associated
%with the automatically-detected fibres in AFM images.
%
%   INPUTS:
%       -axes: The target axes you want to plot into
%       -cMap: The colourmap you want to use to distinguish different
%       samples
%       -normalise: Whether to use raw orientation bin counts or normalised
%       (PDF) values.
%       -data: The data for plotting. In a cell array, with each cell
%       containing a different vector of angles.
%       -labels: Strings defining the labels by which each dataset will be
%       referred to in the legend.
%
%   Author: Oliver J. Meacock, (c) 2023

for i = 1:size(data,1)
    if normalise
        histogram(axes,data{i},'Normalization','pdf','FaceColor',cMap(i,:),'FaceAlpha',0.5,'EdgeColor',cMap(i,:),'lineWidth',1.5)
    else
        histogram(axes,data{i},'FaceColor',cMap(i,:),'FaceAlpha',0.5,'EdgeColor',cMap(i,:),'lineWidth',1.5)
    end
end

legend(axes,labels,'Location','northwest','Interpreter','none')