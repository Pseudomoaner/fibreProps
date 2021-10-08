%Script that converts an AFM heightmap into a strand image. Please ensure
%you have run txtToMat.m first to generate the input matrices.

clear all
close all

root = 'C:\Users\olijm\Desktop\Laia analysis\Raw data';
inputFile = '2020.07.22--WT_internal_XYZ_testdata.mat';
outputSurf = '2020.07.22_FibreReconstruction_Surf.tif';
outputRidge = '2020.07.22_FibreReconstruction_Ridge.tif';
output2Dridge = '2020.07.22_FibreReconstruction_2DRidge.tif';

upscFac = 1; %Set greater than 1 to generate smoother (but slower to process) images
flattenScale = 15; %Scale of the upper DOG filter
ridgeScale = 20; 
ridgeThresh = 0.1;
ridgeSizeFilt = [100,inf]; %Area thresholds used to accept/reject ridges
poreSizeFilt = [25,inf]; %Area thresholds used to accept/reject pores
minBranch = 5;
dil3D = 2;
upscZ = 2; %Factor by which Z should be over-egged compared to the X and Y dimensions in the final output

visualise = true;

load(fullfile(root,inputFile))

origZvals = AFMmat(:,:,3) - min(min(AFMmat(:,:,3))); %Sets the smallest value to zero
upscZvals = imresize(origZvals,upscFac,'bilinear');
flattenedImg = upscZvals - imgaussfilt(upscZvals,flattenScale);
flattenedImg = flattenedImg - min(flattenedImg(:));

ridgeImg = bwRidgeCenterMod(-flattenedImg,ridgeScale,ridgeThresh);
ridgeImg = bwareafilt(ridgeImg,ridgeSizeFilt);
ridgeImg = 1-bwareafilt(logical(1-ridgeImg),poreSizeFilt);

fibreImg = bwskel(logical(imresize(ridgeImg,1/upscFac)),'MinBranchLength',minBranch);

if visualise
    h1 = imshow((flattenedImg-min(flattenedImg(:)))./(max(flattenedImg(:))-min(flattenedImg(:))));
    hold on
    h2 = imshow(fibreImg);
    set(h2,'AlphaData',0.2)
end

heightMap = origZvals; %Choose either origZvals or flattenedImg, depending on whether you want to remove large-scale height fluctuations or not

Zrange = range(heightMap(:));
noZlevels = ceil(Zrange*upscZ/dx);

%Inherit fibre z-positions from original image
[fibreY,fibreX] = ind2sub(size(fibreImg),find(fibreImg));
fibreHeights = round(heightMap(fibreImg)*upscZ/dx);
fibreHeights(fibreHeights<1) = 1;
fibreHeights(fibreHeights>noZlevels) = noZlevels;
fibreInds = sub2ind([size(fibreImg),noZlevels],fibreY,fibreX,fibreHeights);

[surfX,surfY] = meshgrid(1:size(origZvals,2),1:size(origZvals,1));
surfX = reshape(surfX,numel(origZvals),1);
surfY = reshape(surfY,numel(origZvals),1);
surfHeights = round(reshape(heightMap,numel(origZvals),1)*upscZ/dx);
surfHeights(surfHeights<1) = 1;
surfHeights(surfHeights>noZlevels) = noZlevels;
surfInds = sub2ind([size(heightMap),noZlevels],surfY,surfX,surfHeights);

se2D = strel('disk',dil3D);
fibreHeights2D = round(heightMap*upscZ/dx);
fibreHeights2D = fibreHeights2D - min(fibreHeights2D(:)) + 1;
fibreHeights2D(~logical(imdilate(fibreImg,se2D))) = 0;

%Create 3D matrix indicating fibre positions in full 3D space
bin3D = zeros(size(origZvals,1),size(origZvals,2),round(noZlevels*upscFac));
bin3D(fibreInds) = 1;
se = strel('sphere',dil3D);
fibre3D = imdilate(bin3D,se);

%Create 3D matrix indication surface positions in equivalent 3D space
surf3D = zeros(size(origZvals,1),size(origZvals,2),round(noZlevels*upscFac));
surf3D(surfInds) = 1;

imwrite(uint8(fibreHeights2D),fullfile(root,output2Dridge))
for i = 1:size(fibre3D,3)
    imwrite(fibre3D(:,:,i),fullfile(root,outputRidge),'WriteMode','append')
    imwrite(surf3D(:,:,i),fullfile(root,outputSurf),'WriteMode','append')
end