function [AFMmat,dx] = txtToMat(root,inputFile)
%TXTTOMAT converts text format AFM data into a matlab matrix.
%
%   INPUTS:
%       -root: The root directory analysis is being performed in.
%       -inputFile: The name of the specific file you want to analyse
%
%   OUTPUTS:
%       -AFMmat: The AFM data in Matlab matrix format.
%       -dx: The pixel size, in physical units. The physical units can vary
%       depending on how the file was saved, and the units may not be saved
%       ine the file header, so be cautious about interpreting this output.
%
%   Author: Oliver J. Meacock, (c) 2021

fileID = fopen(fullfile(root,inputFile));
C = textscan(fileID,'%f %f %f','Delimiter',' ','CommentStyle','#');
fclose(fileID);

xList = C{1};
yList = C{2};
zList = C{3};

%Find the first reoccurance of an x-coordinate, and mark the distance
%between this position and the the first position as the number of
%x-coordinates. Number of y-coordinates can be calculated from this
xDists = find(xList == xList(1),2);
pxX = xDists(2) - 1;
pxY = round(size(yList,1)/pxX);

%Reshape matrices...
xMat = reshape(xList,pxX,pxY);
yMat = reshape(yList,pxX,pxY);
zMat = reshape(zList,pxX,pxY);

dx = (xMat(end,1)-xList(1,1))/size(xMat,1);

%...and store as a single 3D matrix
AFMmat = cat(3,xMat,yMat,zMat);