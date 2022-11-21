function [bwridge,Width,Nsc,bwridgeCpy] = bwRidgeCenterMod(I, scales, waterThresh, dx)
%BWRIDGECENTERMOD locates ridges of a given scale in an input grayscale
%image.
%
%   INPUTS:
%       -I: The input image, as a matrix of double values.
%       -scales: The scale (width) of ridges you wish to find in I, in
%       pixels. Can be given as a vector of values, in which case a sort of
%       'maximal' ridge magnitude will be calculated across all scales.
%       -waterThresh: The inclusion threshold for the watershed seeds.
%       -dx: The spacing (in nm) between pixels in the image.
%
%   OUTPUTS:
%       -bwridge: Binary image of ridges in I.
%       -Width: Spatial scale that gave the maximal response at each
%       spatial location. In pixels.
%       -Nsc: Scores associated with each maximal response scale
%       -bwridgeCpy: The original image of the fibres as detected by the
%       ridge-detection algorithm. Useful for reconstructing fibres.
%
%   Authors: Joeseph Harvey (c) 2014 and Oliver J. Meacock (c) 2019

progressbar(0)

%Analytical parameters
LpThreshFac = 1.5; %Ridge strength (Lp) detection threshold scaling factor
NThreshFac = 0.6; %Ridge score (N) detection threshold scaling factor. Originally 1 (22/10/2020)
minRidgeArea1 = round(2/(dx^2)); %Ridges in the original ridge-detected image must be at least this large to be included
minRidgeArea2 = round(4/(dx^2)); %Ridges in the pore-filled image must be at least this large to be included
maxPoreArea = round(1.5/(dx^2)); %Pores must be at least this large not to be filled in and removed

% Extract the stationary points of scale-space valleys
[N,Lp] = im_scalablehess2(-I, scales);

progressbar(0.2)

LpThresh = median(abs(Lp(:)))*LpThreshFac;
NThresh = median(N(:))*NThreshFac;

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

progressbar(0.4)

bwridge3D = and(maxPts,stationaryPts);
bwridge = sum(bwridge3D,3);
bwridge = bwareaopen(bwridge,minRidgeArea1,4);
bwridgeCpy = bwridge; %Copy that will be exported from function to be used in fibre reconstruction

%Apply a weak watershed transform to the image:
dists = -graydist(I,bwridge);
distA = imhmin(dists,waterThresh*(prctile(I(:),99)-prctile(I(:),1)));
distW = watershed(distA);

bwridge(distW == 0) = 1;

progressbar(0.6)

%Clean up ridges
bwridge = ~bwareaopen(~bwridge,maxPoreArea,4);
bwridge = bwareaopen(bwridge,minRidgeArea2,4);
bwridge = bwmorph(bwridge,'thin',Inf);
bwridge = bwskel(imclose(bwridge,strel('disk',1)));

%Remove any ridges on image edges
bwridge([1:3,end-2:end],:) = 0;
bwridge(:,[1:3,end-2:end]) = 0;

progressbar(1)