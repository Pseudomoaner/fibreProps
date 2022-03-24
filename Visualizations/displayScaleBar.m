function [] = displayScaleBar(ax,dx,origImg)
%DISPLAYSCALEBAR creates a scalebar in the bottom-right-hand corner of the
%specified axes, based on the current image size.
%
%   INPUTS:
%       -ax: The target axes
%       -dx: The pixel size, in nm
%       -origImg: The original image with pixel size dx
%
%   Author: Oliver J. Meacock, 2022

possBarLens = [0.5,1,2,5,10,20,50,100,200,500]; %In nm

imgWid = size(origImg,2)*dx;

[~,closeInd] = min(abs(possBarLens - imgWid/3));

actualBarLen = possBarLens(closeInd)/dx;

barHeight = size(origImg,1)/50;

rectangle(ax,'Position',[14*size(origImg,2)/15-actualBarLen,19*size(origImg,1)/20,actualBarLen,barHeight],'FaceColor','k')
text(ax,14*size(origImg,2)/15-actualBarLen,14*size(origImg,1)/15,[num2str(possBarLens(closeInd)),' nm'],'FontSize',15,'FontWeight','bold')