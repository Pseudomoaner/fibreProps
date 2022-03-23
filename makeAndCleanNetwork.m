function [outNodes,outLinks] = makeAndCleanNetwork(fibreImg,zImage,Width,Nscores)
%MAKEANDCLEANNETWORK generates a graph of the input AFM fibre network. A 
%ridge-based input is converted to a graph the Skel2Graph3D function. Edits
%are now made to this graph to tidy it up.
%
%   INPUTS:
%       -fibreImg: The ridge-based fibre network estimate, calculated using
%       bwRidgeCenterMod.
%       -zImage: The original AFM image, converted from a .txt file
%       using the txtToMat.m function.
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
%
%   Author: Oliver J. Meacock, (c) 2021

%Convert to graph using Skel2Graph3D
[~,node,link] = Skel2Graph3D(fibreImg,0);

%Reformat links and nodes
node = rmfield(node,{'comz','ep','ptComz'});
for l = 1:size(link,2)
    [link(l).comx,link(l).comy] = ind2sub(size(zImage),link(l).point);
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