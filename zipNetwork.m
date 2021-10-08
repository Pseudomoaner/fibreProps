function [zipNode,zipLink] = zipNetwork(inNode,inLink,inImg)
%ZIPNETWORK removes small links between triply connected nodes (i.e. links 
%in the following configuration: >o-o<). Associated nodes are fused
%together.
%
%   INPUTS:
%       -inNode: Original version of the node structure
%       -inLink: Original version of the link structure
%       -inImg: Original image all metrics are calculated from
%
%   OUPUTS:
%       -zipNode: Version of the node structure with small links zipped out
%       -zipLink: Version of the link structure with small links zipped out
%
%   Author: Oliver J. Meacock, (c) 2020

%Analytical parameters:
minLinkLen = 8; %Links shorter than this (in number of elements) will be targeted

Isz = size(inImg);
upstream = ceil(minLinkLen/2);

%Begin by compiling a list of target links
badLinks = [];
for l = 1:size(inLink,2)
    %Not particularly worried about links joining to end-nodes. Ensure we
    %ignore such links.
    n1 = inLink(l).n1;
    n2 = inLink(l).n2;
    if numel(inNode(n1).links) == 1 || numel(inNode(n2).links) == 1
        skip = true;
    else
        skip = false;
    end
    
    %All the links associated with the ends of this node must also be
    %longer than the minimum length
    testLkList = [inNode(n1).links(inNode(n1).links ~= l),inNode(n2).links(inNode(n2).links ~= l)];
    for tL = testLkList
        if numel(inLink(tL).point) <= ceil(minLinkLen/2)+1
            skip = true;
        end
    end
    
    %And the network of links and nodes at this location must be the
    %archetypal '>-<' structure - ensure no wierld loopy structures etc.
    nodeList1 = inNode(n1).conn;
    nodeList2 = inNode(n2).conn;
    if sum(nodeList1 == n2) > 1 || sum(nodeList2 == n1) > 1
        skip = true;
    else
        nodeList1(nodeList1 == n2) = [];
        nodeList2(nodeList2 == n1) = [];
        allNodes = [nodeList1,nodeList2];
        if numel(unique(allNodes)) ~= numel(allNodes)
            skip = true;
        end
    end
    
    if numel(inLink(l).point) < minLinkLen && ~skip
        badLinks = [badLinks;l];
    end
end

%Now loop through these links, keeping node 1 and moving it to the midpoint
%between the two nodes and interpolating between this position and current
%link endpoints
badNodes = [];
for bL = badLinks'
    gN = inLink(bL).n1;
    bN = inLink(bL).n2;
    
    badNodes = [badNodes;bN];
    
    %Interpolate node position
    inNode(gN).ptComx = mean([inNode(gN).ptComx,inNode(bN).ptComx]);
    inNode(gN).ptComy = mean([inNode(gN).ptComy,inNode(bN).ptComy]);
    inNode(gN).comx = round(inNode(gN).ptComx);
    inNode(gN).comy = round(inNode(gN).ptComy);
    inNode(gN).idx = sub2ind(Isz,inNode(gN).comx,inNode(gN).comy);

    %Reindex node
    inNode(gN).links = unique([inNode(gN).links(inNode(gN).links ~= bL),inNode(bN).links(inNode(bN).links ~= bL)]); %Concatenate the link indices of the 'good' and 'bad' nodes and store in the good node structure
    for i = 1:size(inNode(gN).links,2) %Loop over these links
        aL = inNode(gN).links(i); %Link you want to find the other node associated with (for the 'conn' field)
        testNds = [inLink(aL).n1,inLink(aL).n2]; %One of these will be the good/bad node, the other will be the node we want to add to the conn field
        newTgt = testNds(and(testNds ~= gN,testNds ~= bN)); %Isolate the desired node index (not the good/bad node index)
        inNode(gN).conn(i) = newTgt;
        
        %Reindex new target's own conn field so it no longer points to the
        %bad node (if needed)
        inNode(newTgt).conn(inNode(newTgt).conn == bN) = gN;
    end
    
    %Interpolate associated links
    for gL = inNode(gN).links
        [xList,yList,~] = matchNodeAndLinkEnd(inNode,inLink,gN,gL);
        
        src = [xList(upstream),yList(upstream)];
        tgt = [inNode(gN).ptComx,inNode(gN).ptComy];
        
        interpx = linspace(tgt(1),src(1),upstream*10);
        interpy = linspace(tgt(2),src(2),upstream*10);
        
        roundPts = unique(round([interpx',interpy']),'stable','rows');
        
        inLink(gL).comx = [roundPts(:,1)',xList(upstream+1:end)];
        inLink(gL).comy = [roundPts(:,2)',yList(upstream+1:end)];
        
        inLink(gL).point = sub2ind(Isz,inLink(gL).comx,inLink(gL).comy);
    end
    
    %Reindex associated links
    for gL = inNode(gN).links
        if inLink(gL).n1 == gN || inLink(gL).n1 == bN
            inLink(gL).n1 = gN;
        else
            inLink(gL).n2 = gN;
        end
    end
end

inNode(badNodes) = [];
inLink(badLinks) = [];

%Reindex links and nodes
for i = 1:size(inNode,2)
    for j = 1:size(inNode(i).links,2)
        inNode(i).links(j) = inNode(i).links(j) - sum(badLinks <= inNode(i).links(j));
        inNode(i).conn(j) = inNode(i).conn(j) - sum(badNodes <= inNode(i).conn(j));
    end
end

for i = 1:size(inLink,2)
    inLink(i).n1 = inLink(i).n1 - sum(badNodes <= inLink(i).n1);
    inLink(i).n2 = inLink(i).n2 - sum(badNodes <= inLink(i).n2);
end

zipNode = inNode;
zipLink = inLink;