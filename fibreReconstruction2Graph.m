%Script that converts an AFM heightmap into a 3D projection of the extracted fibre
%network and the orignal surface data.

clear all
close all

root = 'C:\Users\olijm\Desktop\Laia analysis\MachineLearnTest';
inputFile = 'I3_Ori.txt';
outputFib = 'I3_Ori_Fibre3D.tif';
outputSurf = 'I3_Ori_Surf3D.tif';

upscFac = 1; %Set greater than 1 to generate smoother (but slower to process) images
zZoomFac = 2; %Set greater than one to increase the scaling of the z-axis relative to the x and y axes
flattenScale = 15; %Scale of the upper DOG filter

visualise = false;

[AFMmat,dx] = txtToMat(root,inputFile);

origZvals = AFMmat(:,:,3) - min(min(AFMmat(:,:,3))); %Sets the smallest value to zero
upscZvals = imresize(origZvals,upscFac,'bilinear');
flattenedImg = upscZvals - imgaussfilt(upscZvals,flattenScale);
flattenedImg = flattenedImg - min(flattenedImg(:));

[netNodes,netLinks,fibreGroups] = networkReconstruction(flattenedImg,origZvals,1);

%Reconstruct tidied ridge image using graph of network
linkImg = zeros(size(flattenedImg));
for i = 1:size(netLinks,2)
    linkImg(netLinks(i).point) = 1;
end
fibreImg = logical(linkImg);

heightMap = flattenedImg; %Choose either origZvals or flattenedImg, depending on whether you want to remove large-scale height fluctuations or not
heightMap = heightMap*zZoomFac;

Zrange = range(heightMap(:));
noZlevels = ceil(Zrange/dx);

%Inherit fibre z-positions from original image
[fibreY,fibreX] = ind2sub(size(fibreImg),find(fibreImg));
fibreHeights = round(heightMap(fibreImg)/dx);
fibreHeights(fibreHeights<1) = 1;
fibreHeights(fibreHeights>noZlevels) = noZlevels;
fibreInds = sub2ind([size(fibreImg),noZlevels],fibreY,fibreX,fibreHeights);

%Create 3D matrix indicating fibre positions in full 3D space
fibre3D = zeros(size(origZvals,1),size(origZvals,2),round(noZlevels*upscFac));
fibre3D(fibreInds) = 1;

%Create separate 3D matrix indicating original heightmap image
[surfX,surfY] = meshgrid(1:size(heightMap,2),1:size(heightMap,1));
surfX = surfX(:); surfY = surfY(:);
surfHeights = round(heightMap/dx);
surfHeights(surfHeights<1) = 1;
surfHeights(surfHeights>noZlevels) = noZlevels;
surfInds = sub2ind([size(fibreImg),noZlevels],surfY,surfX,surfHeights(:));

surf3D = zeros(size(origZvals,1),size(origZvals,2),round(noZlevels*upscFac));
surf3D(surfInds) = 1;

seSurf = strel('sphere',2);
seFib = strel('sphere',4);
dilatedSurf = imdilate(surf3D,seSurf);
dilatedFib = imdilate(fibre3D,seFib);

%Because you will be in the -append write mode, make sure to delete any
%versions of the binary images that have already been generated (so they
%don't just get appended to)
if exist([root,filesep,outputFib],'file')
    delete([root,filesep,outputFib])
end
if exist([root,filesep,outputSurf],'file')
    delete([root,filesep,outputSurf])
end

for i = 1:size(fibre3D,3)
    imwrite(dilatedFib(:,:,i),[root,filesep,outputFib],'WriteMode','append')
    imwrite(dilatedSurf(:,:,i),[root,filesep,outputSurf],'WriteMode','append')
end

%And finally visualise results (if requested)
if visualise    
    %The original heightmap
    p1 = patch(isosurface(surf3D));
    p1.FaceColor = 'blue';
    p1.EdgeColor = 'none';
    
    %The fibre network
    p2 = patch(isosurface(dilatedFib));
    p2.FaceColor = 'red';
    p2.EdgeColor = 'none';
    axis equal
end