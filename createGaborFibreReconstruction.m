%Script that locates fibres using a Gabor filter bank. Please ensure
%you have run txtToMat.m first to generate the input matrices.

flattenScale = 15; %Scale of the upper DOG filter
angleRange = 0:1:179; %Range of angles over which fiber orientations should be calculated
spatAspRat = 0.2; %Spatial aspect ratio of gabor filterset
spatExtent = 8; %Spatial extent of gabor filterset
magStdThresh = 3; %Number of standard deviations above mean filter output magnitude should be to be segmented out
phaseThresh = 0.5; %Threshodl used to deliniate fibres from non-fibres (phase-based)
minFibre = 100; %Minimum size (in voxels) for a detected 'fibre' to be passed to the next stages

visualise = true;

%If the output image location is already occupied, delete what is there so you
%don't get frames appended on top of it
if exist(fullfile(root,outputImg),'file')
    delete(fullfile(root,outputImg))
end

load(fullfile(root,inputFile))

origZvals = AFMmat(:,:,3) - min(min(AFMmat(:,:,3))); %Sets the smallest value to zero
flattenedImg = origZvals - imgaussfilt(origZvals,flattenScale);
flattenedImg = flattenedImg - min(flattenedImg(:));

flattenedImg = imgaussfilt(flattenedImg,2); %Remove any spatial scale introduced by the rasterization

gaborMagStore = zeros([size(flattenedImg),size(angleRange,2)]);
gaborPhaseStore = zeros([size(flattenedImg),size(angleRange,2)]);
for i = 1:size(angleRange,2)
    [gaborMagStore(:,:,i),gaborPhaseStore(:,:,i)] = imgaborfilt(flattenedImg,spatExtent,angleRange(i),'SpatialAspectRatio',spatAspRat);
end

magThresh = mean(gaborMagStore(:)) + magStdThresh*std(gaborMagStore(:));
gaborStore = zeros(size(origZvals,1),size(origZvals,2),size(angleRange,2));
for i = 1:size(angleRange,2)
    phaseGood = and(gaborPhaseStore(:,:,i)<phaseThresh,gaborPhaseStore(:,:,i)>-phaseThresh); %Ensures you're just getting bright ridges on a dark background
    magGood = gaborMagStore(:,:,i)>magThresh;
    gaborStore(:,:,i) = and(phaseGood,magGood);
end
gaborStore = bwareaopen(gaborStore,minFibre,18);

outStruct = measureGaborFibers(gaborStore,0.1);