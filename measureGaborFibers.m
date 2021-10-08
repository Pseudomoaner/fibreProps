function outStruct = measureGaborFibers(inputMat,dx)
%MEASUREGABORFIBERS creates a structure (outStruct) containing the
%locations and measurements of the fibres detected by the Gabor filterbank.
%
%   INPUTS:
%       -inputMat: Binary array indicating the spatial coordinates in the
%       first two dimensions and the Gabor filter's orientation in the
%       third. Should wrap around (so that inputMat(:,:,end) leads directly
%       onto inputMat(:,:,1)).
%       -dx: Distance between sampling points in the array. Assumed that
%       the image is isotropic.
%
%   Author: Oliver J. Meacock, (c) 2020

%To take care of the 'seam', we will need to concatinate the start of the
%input array to the end of the input array (in dimension 3, the angular
%direction) and then eliminate any detected objects with identical
%centroids.
noAngs = size(inputMat,3);
buffAngs = round(noAngs/5); %Using 20% of input angles as the buffer
inputMat = cat(3,inputMat,inputMat(:,:,1:buffAngs));
inputMat = imclearborder(inputMat,18); %clear object touching borders to suppress edge effects and partially detected fibers

CCs = bwconncomp(logical(inputMat),18);
props3D = regionprops3(CCs,'Centroid','VoxelList');

%Go through all the objects, and remove the second instance of any with
%repeated centroids
badIDs = [];
cents = double.empty(0,2);
for i = 1:size(props3D,1)
    repInds = find(sum(cents == repmat(props3D.Centroid(i,1:2),size(cents,1),1),2) == 2);
    if ~isempty(repInds)
        badIDs = [badIDs;i];
    else
        cents = [cents;props3D.Centroid(i,1:2)];
    end
end

props3D(badIDs,:) = [];

%Find properties of remaining objects, including average orientation and 2D
%projection of pixel profile
outStruct = struct('Orientation',[],'PixelList',[],'Centroid',[]);
dAng = 180/noAngs;
for i = 1:size(props3D,1)
    outStruct(i).Centroid = props3D.Centroid(i,1:2)*dx;
    outStruct(i).Orientation = wrapTo180(props3D.Centroid(i,3)*dAng);
    if outStruct(i).Orientation < 0
        outStruct(i).Orientation = outStruct(i).Orientation + 180;
    end
    outStruct(i).Orientation = mod(outStruct(i).Orientation+90,180)-90;
    
    outStruct(i).PixelList = unique(props3D.VoxelList{i}(:,1:2),'rows');
    outStruct(i).PixelIndList = sub2ind([size(inputMat,1),size(inputMat,2)],outStruct(i).PixelList(:,1),outStruct(i).PixelList(:,2));
end