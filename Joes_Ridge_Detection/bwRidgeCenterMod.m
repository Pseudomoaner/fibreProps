function [bwridge,Width,Nsc,bwridgeCpy] = bwRidgeCenterMod(I,  scales)
%BWRIDGECENTERMOD locates ridges of a given scale in an input grayscale
%image.
%
%   INPUTS:
%       -I: The input image, as a matrix of double values.
%       -scales: The scale (width) of ridges you wish to find in I, in
%       pixels. Can be given as a vector of values, in which case a sort of
%       'maximal' ridge magnitude will be calculated across all scales.
%
%   OUTPUTS:
%       -bwridge: Binary image of ridges in I.
%       -Width: Spatial scale that gave the maximal response at each
%       spatial location.
%
%   Authors: Joeseph Harvey (c) 2014 and Oliver J. Meacock (c) 2019

%Analytical parameters
LpThreshFac = 1.5;
NThreshFac = 0.85; %Originally 1 (22/10/2020)
minRidgeArea1 = 15;
minRidgeArea2 = 30;
maxPoreArea = 20;
waterThresh = 3;

% Extract the stationary points of scale-space valleys
[N,Lp] = im_scalablehess2(-I, scales);

LpThresh = mean(abs(Lp(:)))*LpThreshFac;
NThresh = mean(N(:))*NThreshFac;

nearPeaks = Lp;
nearPeaks(or(nearPeaks>LpThresh,nearPeaks<-LpThresh)) = nan;
A = nearPeaks < 0;
B = nearPeaks > 0;

se = strel('sphere',1);
Adil = imdilate(A,se);
Bdil = imdilate(B,se);

stationaryPts = and(Adil,Bdil);

%Find maximal N-scores for each pixels and associated scales
[Nsc, MaxIdx] = max(N, [], 3);
Width = scales(MaxIdx);

maxPts = zeros(size(stationaryPts));

%Fairly certain there's a faster/more elegant way of doing this...
for i = 1:size(maxPts,1)
    for j = 1:size(maxPts,2)
        if MaxIdx(i,j) ~= 1 && Nsc(i,j) > NThresh
            maxPts(i,j,MaxIdx(i,j)) = 1;
        end
    end
end

bwridge3D = and(maxPts,stationaryPts);
bwridge = sum(bwridge3D,3);
bwridge = bwareaopen(bwridge,minRidgeArea1,4);
bwridgeCpy = bwridge; %Copy that will be exported from function

%Apply a weak watershed transform to the image:
dists = -bwdist(bwridge);
distA = imhmin(dists,waterThresh);
distW = watershed(distA);

bwridge(distW == 0) = 1;

%Clean up ridges
bwridge = ~bwareaopen(~bwridge,maxPoreArea,4);
bwridge = bwareaopen(bwridge,minRidgeArea2,4);
bwridge = bwmorph(bwridge,'thin',Inf);
bwridge = bwskel(imclose(bwridge,strel('disk',1)));