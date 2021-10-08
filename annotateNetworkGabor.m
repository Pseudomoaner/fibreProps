function [noteNode,noteLink] = annotateNetworkGabor(inNode,inLink,fibreStruc,zImg)
%ANNOTATENETWORKGABOR provides annotations of the input graph based on the input
%data, including fibre assignments.
%
%   INPUTS:
%       -inNode: input set of nodes
%       -inLink: input set of links between nodes
%       -fibreStruc: Structure containing fibre data. Output by
%       gaborFibreReconstruction.m.
%       -zImg: The original image containing the z-heights of the fibre
%       network
%
%   OUTPUTS:
%       -noteNode: The version of the input nodes with annotations added
%       -noteLink: The version of the input links with annotations added
%
%   Author: Oliver J. Meacock, (c) 2020

%Analytical parameters
fibreRad = 5; %Amount fibre profile will be dilated by
angularTol = 45; %Range of angles in degrees (plus and minus) that a link should be wrt fibre to be counted as belonging to that fiber 

%Initialise fibre index storage
for i = 1:size(inLink,2)
    inLink(i).Fibres = [];
end

%Begin by drawing out the network of links, using a separate numerical
%label for each one
linkImg = zeros(size(zImg));
for l = 1:size(inLink,2)
    linkImg(inLink(l).point) = l;
end

%Now dilate each fibre profile, and see which links fit within
se = strel('disk',fibreRad);
for f = 1:size(fibreStruc,2)
    fibImg = zeros(size(zImg));
    fibImg(fibreStruc(f).PixelIndList) = 1;
    fibImg = imdilate(fibImg,se)';
    
    currLinks = linkImg;
    currLinks(~logical(fibImg)) = 0;
    
    %For each link included now in currLinks, check that it fulfills
    %certain criteria to be annotated as part of this fibre
    currLinkList = unique(currLinks(:));
    currLinkList = currLinkList(2:end); %First element is always 0, which we can remove
    if ~isempty(currLinkList)
        for l = currLinkList'
            goodLink = true;
            
            %Check 1: Is at least half of link inside (dilated) fiber
            %perimeter?
            inFib = sum(currLinks(:) == l);
            if inFib < numel(inLink(l).point)/2
                goodLink = false;
            end
            
            %Check 2: Is the fibre (roughly) oriented in the same direction as the link?
            nd1Pos = [inNode(inLink(l).n1).comx,inNode(inLink(l).n1).comy];
            nd2Pos = [inNode(inLink(l).n2).comx,inNode(inLink(l).n2).comy];
            lkAng = atand((nd1Pos(2) - nd2Pos(2))/(nd1Pos(1) - nd2Pos(1)));
            angDiff = mod(lkAng - fibreStruc(f).Orientation+270,180)-90;
            if abs(angDiff) > angularTol
                goodLink = false;
            end
            
            %If this link passes checks, annotate it with this fibre index
            if goodLink
                inLink(l).Fibres = [inLink(1).Fibres;f];
            end
        end
    end
end

noteNode = inNode;
noteLink = inLink;