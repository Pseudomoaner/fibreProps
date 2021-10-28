function [] = plotFibreOrientationMass(fibreOrientations,fibreSizes)
%PLOTFIBREORIENTATIONMASS plots the total amount of 'mass' (fibre length)
%present in the input AFM fibre data in a given orientational bin.
%   INPUTS:
%       -fibreOrientations: The orientations of the fibres in the image,
%       output as a field of fibreProps by measureFibres. An Nx1 vector.
%       -fibreSizes: The total length of each fibre, output as a field of
%       fibreProps by measureFibres. An Nx1 vector.
%
%   Author: Oliver J. Meacock, (c) 2021

noBins = 18;

angSet = linspace(-pi/2,pi/2,noBins);
[~,~,angBin] = histcounts(vertcat(fibreOrientations),'BinEdges',angSet);
binMasses = accumarray(angBin,fibreSizes);

polarplot(angSet(1:end-1) + diff(angSet(1:2))/2,binMasses,'r','lineWidth',1.5)

ax = gca;
ax.ThetaLim = [-90,90];