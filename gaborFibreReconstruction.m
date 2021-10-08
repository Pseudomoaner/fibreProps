function outStruct = gaborFibreReconstruction(flatImage,dx)
%GABORFIBRERECONSTRUCTION applies a Gabor filterbank to the input image to
%locate long, strandlike structures within it. Then measures these
%structures in a separate script.
%
%   INPUTS:
%       -flatImage: Grayscale image containing fibers.
%       -dx: Distance between adjacent pixels (in um)
%
%   OUTPUTS:
%       -outStruct: Structure containing locations, orientations and 2D
%       reconstructions of individual fibers.
%
%   Author: Oliver J. Meacock, (c) 2020

%Analytical parameters
angleRange = 0:1:179; %Range of angles over which fiber orientations should be calculated
spatAspRat = 0.2; %Spatial aspect ratio of gabor filterset
spatExtent = 8; %Spatial extent of gabor filterset
magStdThresh = 3; %Number of standard deviations above mean filter output magnitude should be to be segmented out
phaseThresh = 0.5; %Threshold used to deliniate fibres from non-fibres (phase-based)
minFibre = 100; %Minimum size (in voxels) for a detected 'fibre' to be passed to the next stages

flatImage = imgaussfilt(flatImage,2); %Remove any spatial scale introduced by the rasterization

gaborMagStore = zeros([size(flatImage),size(angleRange,2)]);
gaborPhaseStore = zeros([size(flatImage),size(angleRange,2)]);
for i = 1:size(angleRange,2)
    [gaborMagStore(:,:,i),gaborPhaseStore(:,:,i)] = imgaborfilt(flatImage,spatExtent,angleRange(i),'SpatialAspectRatio',spatAspRat);
end

magThresh = mean(gaborMagStore(:)) + magStdThresh*std(gaborMagStore(:));
gaborStore = zeros(size(flatImage,1),size(flatImage,2),size(angleRange,2));
for i = 1:size(angleRange,2)
    phaseGood = and(gaborPhaseStore(:,:,i)<phaseThresh,gaborPhaseStore(:,:,i)>-phaseThresh); %Ensures you're just getting bright ridges on a dark background
    magGood = gaborMagStore(:,:,i)>magThresh;
    gaborStore(:,:,i) = and(phaseGood,magGood);
end
gaborStore = bwareaopen(gaborStore,minFibre,18);

outStruct = measureGaborFibers(gaborStore,dx);