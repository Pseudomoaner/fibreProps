function [outNodes,outLinks,fibreGroups] = networkReconstruction(flatImage,zImage,dx,widFac)
%NETWORKRECONSTRUCTION generates a graph of the input AFM fibre network. In
%the first stage, a scale-selecting ridge-detection algorithm is used to detect the
%locations and widths of fibres. The resulting network is then skeletonized
%and converted to a graph using the Skel2Graph3D function. Edits are now
%made to this graph to tidy it up.
%
%   INPUTS:
%       -flatImage: The flattened version of-
%       -zImage: The original AFM image, converted from a .txt file
%       using the txtToMat.m function.
%       -dx: The spacing (in nm) between pixels in the image.
%       -widFac: Correction factor used to convert detected fibre widths 
%       based on the scale of a Gaussian filter) to true widths. Requires
%       manual calibration.
%
%   OUTPUTS:
%       -outNodes: Structure conforming (mostly) to the formatting of the
%       node output of the Skel2Graph3D function. Contains the information 
%       (positions, heights, associated links etc.) about the nodes in the
%       reconstruted network graph.
%       -outLinks: Structure conforming (mostly) to the formatting of the
%       link oubput of the Skel2Graph3D function. Contains the information 
%       (paths, heights, associated nodes etc.) about the links in the
%       reconstruted network graph. 
%       -fibreGroup: Image of the same size as the flatImage/zImage inputs,
%       containing labelled segmentations of the original fibres as
%       detected by the ridge detection process. Used later for assigning
%       links and nodes to different fibres.
%
%   Author: Oliver J. Meacock, (c) 2021

%Analysis parameters
minScale = round(max(15,10/dx)); %Ensures the minimum scale is at least 12 pixels - ridge detection is pretty unreliable below this value
maxScale = round(max(45,30/dx));
ridgeScale = linspace(minScale,maxScale,20); %Scales are in nm; converted to pixels  
minBranch = round(2/dx); %Defines the minimal length of an isolated branch - branches shorter than this are excluded

%Create image of ridges in image, and sekeletonize
[ridgeImg,Width,Nscores,fibreGroups] = bwRidgeCenterMod(flatImage,ridgeScale,dx,widFac);
fibreImg = bwskel(logical(ridgeImg),'MinBranchLength',minBranch);

%Convert to graph using Skel2Graph3D
[~,node,link] = Skel2Graph3D(fibreImg,0);

%Reformat links and nodes
node = rmfield(node,{'comz','ep','ptComz'});
for l = 1:size(link,2)
    [link(l).comx,link(l).comy] = ind2sub(size(flatImage),link(l).point);
end

%This function 'heals' the graph, removing individual nodes that have only
%two links and joining the corresponding links together. Seems to be a
%fault in the Skel2Graph3D code that it spits these out.
[healNode,healLink] = healNetwork(node,link,false);

%This function 'zips together' pairs of nodes that are too close together
%for the link between them to be reliably analysed
[zipNode,zipLink] = zipNetwork(healNode,healLink,zImage);

%We also want to associate the original z-values to links and nodes.
%Inherit from original image.

%Nodes first
for n = 1:size(zipNode,2)
    zipNode(n).comz= zImage(zipNode(n).idx);
    zipNode(n).widths = Width(zipNode(n).idx);
    zipNode(n).Nscore = Nscores(zipNode(n).idx);
end

%And then links
for l = 1:size(zipLink,2)
    zipLink(l).comz = zImage(zipLink(l).point);
    zipLink(l).widths = Width(zipLink(l).point);
    zipLink(l).Nscore = Nscores(zipLink(l).point);
end

outNodes = zipNode;
outLinks = zipLink;