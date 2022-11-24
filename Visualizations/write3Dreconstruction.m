function [] = write3Dreconstruction(fibreProps,heightmap,upscZ,dx,destination)
%WRITE3DRECONSTRUCTION creates a 3D reconstruction of the input set of
%fibres based on the selected heightmap.
%
%   INPUTS:
%       -fibreProps: A field of the fibreProps structure output by the
%       fibreFinder GUI.
%       -heightmap: The source heigtmap you want to read the z-coordinates
%       for each position on the fibre from. Typically either taken from
%       the origImgs or flatImgs outputs of fibreFinder.
%       -upscZ: How much you want to upscale the z-dimension relative to
%       the x and y-dimensions to amplify differences in ridge height.
%       -dx: The pixel size of the image (0.5 for fibreFinder)
%       -destination: The target output filename (without file extension).
%
%   Author: Oliver J. Meacock, 2022

%Determine the number of layers this reconstruction will need
Zrange = range(heightmap(:));
noZlevels = ceil(Zrange*upscZ/dx);

%Reconstruct a 2d fibre image
fibreImg = false(size(heightmap));
for i = 1:size(fibreProps,2)
    fibreImg(sub2ind(size(heightmap),fibreProps(i).backList(:,2),fibreProps(i).backList(:,1))) = true;
end

%Inherit fibre z-positions from original image
[fibreY,fibreX] = ind2sub(size(fibreImg),find(fibreImg));
fibreHeights = round(heightmap(fibreImg)*upscZ/dx);
fibreHeights(fibreHeights<1) = 1;
fibreHeights(fibreHeights>noZlevels) = noZlevels;
fibreInds = sub2ind([size(fibreImg),noZlevels],fibreY,fibreX,fibreHeights);

[surfX,surfY] = meshgrid(1:size(heightmap,2),1:size(heightmap,1));
surfX = reshape(surfX,numel(heightmap),1);
surfY = reshape(surfY,numel(heightmap),1);
surfHeights = round(reshape(heightmap,numel(heightmap),1)*upscZ/dx);
surfHeights(surfHeights<1) = 1;
surfHeights(surfHeights>noZlevels) = noZlevels;
surfInds = sub2ind([size(heightmap),noZlevels],surfY,surfX,surfHeights);

%Create 3D matrix indicating fibre positions in full 3D space
bin3D = zeros(size(heightmap,1),size(heightmap,2),noZlevels);
bin3D(fibreInds) = 1;
se = strel('sphere',5);
fibre3D = imdilate(bin3D,se);

%Create 3D matrix indication surface positions in equivalent 3D space
surf3D = zeros(size(heightmap,1),size(heightmap,2),noZlevels);
surf3D(surfInds) = 1;

%Write 3D images
for i = 1:size(fibre3D,3)
    if i == 1
        imwrite(fibre3D(:,:,i),[destination,'_Fibres.tif'])
        imwrite(surf3D(:,:,i),[destination,'_Surf.tif'])        
    else
    imwrite(fibre3D(:,:,i),[destination,'_Fibres.tif'],'WriteMode','append')
    imwrite(surf3D(:,:,i),[destination,'_Surf.tif'],'WriteMode','append')
    end
end