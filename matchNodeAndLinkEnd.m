function [outXlist,outYlist,outPtList,flipFlag] = matchNodeAndLinkEnd(inNode,inLink,nInd,lInd)
%MATCHNODEANDLINKEND finds which of the two ends of a given link is closest
%to a given node, and returns the link coordinates with the node end first
%
%   INPUTS:
%       -inNode: Node structure
%       -inLink: Link structure
%       -nInd: Index of the node of interest
%       -lInd: Index of the link of interest
%
%   OUTPUTS:
%       -outXlist: List of link x-coordinates, node-closest element first
%       -outYlist: List of link y-coordinates, node-closest element first
%       -outPtList: List of link pixel indices, node-closest element first
%       -flipFlag: Whether you've applied flipping to the above data or not
%
%   Author: Oliver J. Meacock, (c) 2020

distE1 = sqrt((inLink(lInd).comx(1) - inNode(nInd).ptComx)^2 + (inLink(lInd).comy(1) - inNode(nInd).ptComy)^2);
distE2 = sqrt((inLink(lInd).comx(end) - inNode(nInd).ptComx)^2 + (inLink(lInd).comy(end) - inNode(nInd).ptComy)^2);

if distE1 > distE2
    outXlist = flip(inLink(lInd).comx);
    outYlist = flip(inLink(lInd).comy);
    outPtList = flip(inLink(lInd).point);
    flipFlag = true;
else
    outXlist = inLink(lInd).comx;
    outYlist = inLink(lInd).comy;
    outPtList = inLink(lInd).point;
    flipFlag = false;
end