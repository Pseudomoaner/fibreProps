function [] = plotFibreBoxplot(axes,data,labels)
%PLOTFIBREBOXPLOT generates a boxplot of the sample-wide averages for each
%condition.
%
%   INPUTS:
%       -axes: The target axes you want to plot into
%       -data: The data for plotting. In a cell array, with each cell
%       containing a different vector of sample-wide averages.
%       -labels: Strings defining the labels by which each dataset will be
%       referred to in the legend.
%
%   Author: Oliver J. Meacock, (c) 2021

if ~isempty(data)
    %Convert data into boxplot-compatible format
    plotDat = horzcat(data{:});
    groups = cell(numel(plotDat),1);
    runSum = 1;
    for i = 1:size(data,1)
        [groups{runSum:runSum+numel(data{i})-1}] = deal(labels{i});
        runSum = runSum + numel(data{i});
    end
    boxchart(axes,categorical(groups),plotDat)
    set(axes,'TickLabelInterpreter','none');
    set(axes,'XTickLabelRotation',30);
end