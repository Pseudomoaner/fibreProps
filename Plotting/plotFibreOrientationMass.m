function [] = plotFibreOrientationMass(fibreOrientations,fibreSizes)
%PLOTFIBREORIENTATIONMASS plots the total amount of 'mass' (fibre length)
%present in the input AFM fibre data in a given orientational bin.
%
%   INPUTS:
%       -fibreOrientations: The orientations of the fibres in the image,
%       output as a field of fibreProps by measureFibres. An Nx1 vector.
%       -fibreSizes: The total length of each fibre, output as a field of
%       fibreProps by measureFibres. An Nx1 vector.
%
%   Author: Oliver J. Meacock, (c) 2021

noBins = 18;

angSet = linspace(0,pi,noBins+1);
[~,~,angBin] = histcounts(vertcat(fibreOrientations)+pi/2,'BinEdges',angSet);
binMasses = accumarray(angBin',fibreSizes');

%Need to pad out binMasses with zeros if the final bins didn't contain
%anything
binMasses = [binMasses;zeros(noBins-size(binMasses,1)-1,1)];

%Repeat binMasses so the plot goes the full 360 degrees.
binMasses = [binMasses;binMasses;binMasses(1)];
angList = (angSet(1:end-1) + diff(angSet(1:2))/2)';
angList = [angList;angList + pi;angList(1)];

polarplot(angList,binMasses,'r','lineWidth',1.5)

ax = gca;
ax.ThetaLim = [0,180];