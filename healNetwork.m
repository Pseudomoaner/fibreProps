function [healNode,healLink] = healNetwork(inNode,inLink,visualise)
%HEALNETWORK makes the following adjustments to a network graph:
%   - Removes isolated nodes
%   - Removes nodes connected to two links, and fuses the two links
%   together
%   - Removes nodes connected to a single, very short link
%   - Ensures each pair of nodes is connected by a single link, removing
%   any excess links.
%
%   INPUTS:
%       -inNode: The original node structure
%       -inLink: The original link structure
%
%   OUTPUTS:
%       -healNode: The node structure with double-nodes removed
%       -healLink: The link structure with links associated with
%       double-nodes fused together
%
%   Author: Oliver J. Meacock, (c) 2020

minSingNodeEdgeLen = 15; %How long the link associated with a singly connected node should be for both to be retained

%% Part 0: Recast various fields to avoid use of awkward uint16 type
for l = 1:size(inLink,2)
    inLink(l).n1 = double(inLink(l).n1);
    inLink(l).n2 = double(inLink(l).n2);
end

%% Part 1: Remove isolated nodes
badNodes = [];
for i = 1:size(inNode,2)
    if isempty(inNode(i).links)
        badNodes = [badNodes;i];
    end
end
inNode(badNodes) = [];

%Need to reindex nodes associated with each link
for i = 1:size(inLink,2)
    inLink(i).n1 = inLink(i).n1 - sum(badNodes <= inLink(i).n1);
    inLink(i).n2 = inLink(i).n2 - sum(badNodes <= inLink(i).n2);
end

%And also the nodes associated with the .conn field of the node structure
for i = 1:size(inNode,2)
    for j = 1:size(inNode(i).conn,2)
        inNode(i).conn(j) = inNode(i).conn(j) - sum(badNodes <= inNode(i).conn(j));
    end
end

%% Part 2: Remove nodes connecting two links
badNodes = [];
for i = 1:size(inNode,2)
    if numel(inNode(i).links) == 2
        badNodes = [badNodes;i];
    end
end

%For each bad node, one link will also need to be deleted. Associated data
%will need to be inherited by the remaining link.
badLinks = [];
for i = 1:size(badNodes,1)
    bN = badNodes(i);
    bL = inNode(bN).links(2);
    gL = inNode(bN).links(1);
    
    [glComx,glComy,glPt] = matchNodeAndLinkEnd(inNode,inLink,bN,gL);
    [blComx,blComy,blPt] = matchNodeAndLinkEnd(inNode,inLink,bN,bL);
    
    %Storage order will be [good link, bad node, bad link]
    glPt = flip(glPt);
    glComx = flip(glComx);
    glComy = flip(glComy);
    
    inLink(gL).point = [glPt,inNode(bN).idx',blPt];
    inLink(gL).comx = [glComx,inNode(bN).comx',blComx];
    inLink(gL).comy = [glComy,inNode(bN).comy',blComy];
    
    %Reindex connections between kept link and good nodes
    if inLink(bL).n1 == bN
        newTgt = inLink(bL).n2;
    else
        newTgt = inLink(bL).n1;
    end
    if inLink(gL).n1 == bN
        inLink(gL).n1 = newTgt;
        src = inLink(gL).n2;
    else
        inLink(gL).n2 = newTgt;
        src = inLink(gL).n1;
    end
    tgtNdChangeInd = inNode(newTgt).conn == bN;
    srcNdChangeInd = inNode(src).conn == bN;
    inNode(newTgt).conn(tgtNdChangeInd) = src;
    inNode(src).conn(srcNdChangeInd) = newTgt;
    inNode(newTgt).links(tgtNdChangeInd) = gL;
    inNode(src).links(srcNdChangeInd) = gL;
    
    badLinks = [badLinks;bL];
end

%Remove bad links and nodes
inNode(badNodes) = [];
inLink(badLinks) = [];

%Reindex node links and nodes
for i = 1:size(inNode,2)
    for j = 1:size(inNode(i).links,2)
        inNode(i).links(j) = inNode(i).links(j) - sum(badLinks <= inNode(i).links(j));
        inNode(i).conn(j) = inNode(i).conn(j) - sum(badNodes <= inNode(i).conn(j));
    end
end

%And reindex link nodes
for i = 1:size(inLink,2)
    inLink(i).n1 = inLink(i).n1 - sum(badNodes <= inLink(i).n1);
    inLink(i).n2 = inLink(i).n2 - sum(badNodes <= inLink(i).n2);
end

%% Part 3: Remove nodes connecting to a single, very short edge, along with associated links
badNodes = [];
for i = 1:size(inNode,2)
    if numel(inNode(i).links) == 1
        if numel(inLink(inNode(i).links).point) < minSingNodeEdgeLen
            badNodes = [badNodes;i];
        end
    end
end

badLinks = [];
for i = 1:size(badNodes,1)
    bN = badNodes(i);
    bL = inNode(bN).links(1);
    gN = inNode(bN).conn(1);
    
    if sum(badNodes == inNode(bN).conn) == 0 %Special case if you have two singly connected nodes joined by a single short edge. Both will end up deleted anyway, so don't worry about changing associated references        
        %Delete references to bad node and link
        inNode(gN).conn(inNode(gN).conn == bN) = [];
        inNode(gN).links(inNode(gN).links == bL) = [];
    end
    
    badLinks = [badLinks;bL];
end

badLinks = unique(badLinks);

%Remove bad links and nodes
inNode(badNodes) = [];
inLink(badLinks) = [];

%Reindex node links and nodes
for i = 1:size(inNode,2)
    for j = 1:size(inNode(i).links,2)
        inNode(i).links(j) = inNode(i).links(j) - sum(badLinks <= inNode(i).links(j));
        inNode(i).conn(j) = inNode(i).conn(j) - sum(badNodes <= inNode(i).conn(j));
    end
end

%And reindex link nodes
for i = 1:size(inLink,2)
    inLink(i).n1 = inLink(i).n1 - sum(badNodes <= inLink(i).n1);
    inLink(i).n2 = inLink(i).n2 - sum(badNodes <= inLink(i).n2);
end

%% Part 4: Repeat double link node removal, now singly linked nodes have been removed
badNodes = [];
for i = 1:size(inNode,2)
    if numel(inNode(i).links) == 2
        badNodes = [badNodes;i];
    end
end

%For each bad node, one link will also need to be deleted. Associated data
%will need to be inherited by the remaining link.
badLinks = [];
for i = 1:size(badNodes,1)
    bN = badNodes(i);
    bL = inNode(bN).links(2);
    gL = inNode(bN).links(1);
    
    [glComx,glComy,glPt] = matchNodeAndLinkEnd(inNode,inLink,bN,gL);
    [blComx,blComy,blPt] = matchNodeAndLinkEnd(inNode,inLink,bN,bL);
    
    %Storage order will be [good link, bad node, bad link]
    glPt = flip(glPt);
    glComx = flip(glComx);
    glComy = flip(glComy);
    
    inLink(gL).point = [glPt,inNode(bN).idx',blPt];
    inLink(gL).comx = [glComx,inNode(bN).comx',blComx];
    inLink(gL).comy = [glComy,inNode(bN).comy',blComy];
    
    %Reindex connections between kept link and good nodes
    if inLink(bL).n1 == bN
        newTgt = inLink(bL).n2;
    else
        newTgt = inLink(bL).n1;
    end
    if inLink(gL).n1 == bN
        inLink(gL).n1 = newTgt;
        src = inLink(gL).n2;
    else
        inLink(gL).n2 = newTgt;
        src = inLink(gL).n1;
    end
    tgtNdChangeInd = inNode(newTgt).conn == bN;
    srcNdChangeInd = inNode(src).conn == bN;
    inNode(newTgt).conn(tgtNdChangeInd) = src;
    inNode(src).conn(srcNdChangeInd) = newTgt;
    inNode(newTgt).links(tgtNdChangeInd) = gL;
    inNode(src).links(srcNdChangeInd) = gL;
    
    badLinks = [badLinks;bL];
end

%Remove bad links and nodes
inNode(badNodes) = [];
inLink(badLinks) = [];

%Reindex node links and nodes
for i = 1:size(inNode,2)
    for j = 1:size(inNode(i).links,2)
        inNode(i).links(j) = inNode(i).links(j) - sum(badLinks <= inNode(i).links(j));
        inNode(i).conn(j) = inNode(i).conn(j) - sum(badNodes <= inNode(i).conn(j));
    end
end

%And reindex link nodes
for i = 1:size(inLink,2)
    inLink(i).n1 = inLink(i).n1 - sum(badNodes <= inLink(i).n1);
    inLink(i).n2 = inLink(i).n2 - sum(badNodes <= inLink(i).n2);
end

%% Part 5: Remove links that connect back to a single node
%Locate the indices of guilty links
nList = zeros(size(inLink,2),2);
for l = 1:size(inLink,2)
    nList(l,1) = inLink(l).n1;
    nList(l,2) = inLink(l).n2;
end

badLinks = find(nList(:,1) == nList(:,2));

%Remove bad links
inLink(badLinks) = [];

%Reindex node links, and remove any necessary links
for i = 1:size(inNode,2)
    remLinks = [];
    for j = 1:size(inNode(i).links,2)
        if sum(badLinks == inNode(i).links(j)) == 1
            remLinks = [remLinks;j];
        else
            inNode(i).links(j) = inNode(i).links(j) - sum(badLinks < inNode(i).links(j));
            possNds = [inLink(inNode(i).links(j)).n1,inLink(inNode(i).links(j)).n2];
            inNode(i).conn(j) = possNds(possNds ~= i);
        end
    end
    inNode(i).links(remLinks) = [];
    inNode(i).conn(remLinks) = [];
end

%% Part 6: Remove all but one link between identical nodes (e.g. cases where both o^o and o-o are true)
%Strictly speaking, these instances could correspond to real situations. But
%greatly simplifies downstream analysis if you chuck out the rare instances
%where this happens.

%Locate the indices of guilty links
nList = zeros(size(inLink,2),2);
for l = 1:size(inLink,2)
    nList(l,1) = inLink(l).n1;
    nList(l,2) = inLink(l).n2;
end

[~,ia,~] = unique(nList,'rows');
badLinks = setdiff(1:size(inLink,2),ia);

%Remove bad links
inLink(badLinks) = [];

%Reindex node links, and remove any necessary links
for i = 1:size(inNode,2)
    remLinks = [];
    for j = 1:size(inNode(i).links,2)
        if sum(badLinks == inNode(i).links(j)) == 1
            remLinks = [remLinks;j];
        else
            inNode(i).links(j) = inNode(i).links(j) - sum(badLinks < inNode(i).links(j));
            possNds = [inLink(inNode(i).links(j)).n1,inLink(inNode(i).links(j)).n2];
            inNode(i).conn(j) = possNds(possNds ~= i);
        end
    end
    inNode(i).links(remLinks) = [];
    inNode(i).conn(remLinks) = [];
end

%% Part 7: Repeat double link node removal, now 'eyes' have been removed
badNodes = [];
for i = 1:size(inNode,2)
    if numel(inNode(i).links) == 2
        badNodes = [badNodes;i];
    end
end

%For each bad node, one link will also need to be deleted. Associated data
%will need to be inherited by the remaining link.
badLinks = [];
for i = 1:size(badNodes,1)
    bN = badNodes(i);
    bL = inNode(bN).links(2);
    gL = inNode(bN).links(1);
    
    [glComx,glComy,glPt] = matchNodeAndLinkEnd(inNode,inLink,bN,gL);
    [blComx,blComy,blPt] = matchNodeAndLinkEnd(inNode,inLink,bN,bL);
    
    %Storage order will be [good link, bad node, bad link]
    glPt = flip(glPt);
    glComx = flip(glComx);
    glComy = flip(glComy);
    
    inLink(gL).point = [glPt,inNode(bN).idx',blPt];
    inLink(gL).comx = [glComx,inNode(bN).comx',blComx];
    inLink(gL).comy = [glComy,inNode(bN).comy',blComy];
    
    %Reindex connections between kept link and good nodes
    if inLink(bL).n1 == bN
        newTgt = inLink(bL).n2;
    else
        newTgt = inLink(bL).n1;
    end
    if inLink(gL).n1 == bN
        inLink(gL).n1 = newTgt;
        src = inLink(gL).n2;
    else
        inLink(gL).n2 = newTgt;
        src = inLink(gL).n1;
    end
    tgtNdChangeInd = inNode(newTgt).conn == bN;
    srcNdChangeInd = inNode(src).conn == bN;
    inNode(newTgt).conn(tgtNdChangeInd) = src;
    inNode(src).conn(srcNdChangeInd) = newTgt;
    inNode(newTgt).links(tgtNdChangeInd) = gL;
    inNode(src).links(srcNdChangeInd) = gL;
    
    badLinks = [badLinks;bL];
end

%Remove bad links and nodes
inNode(badNodes) = [];
inLink(badLinks) = [];

%Reindex node links and nodes
for i = 1:size(inNode,2)
    for j = 1:size(inNode(i).links,2)
        inNode(i).links(j) = inNode(i).links(j) - sum(badLinks <= inNode(i).links(j));
        inNode(i).conn(j) = inNode(i).conn(j) - sum(badNodes <= inNode(i).conn(j));
    end
end

%And reindex link nodes
for i = 1:size(inLink,2)
    inLink(i).n1 = inLink(i).n1 - sum(badNodes <= inLink(i).n1);
    inLink(i).n2 = inLink(i).n2 - sum(badNodes <= inLink(i).n2);
end

%% Part 7: Optional: Validation of output
if visualise
    figure(1)
    hold on
    for n = 1:size(inNode,2)
        plot(inNode(n).comx(1),inNode(n).comy(1),'r.')
        text(inNode(n).comx(1),inNode(n).comy(1),num2str(n),'Color','r')
    end
    for l = 1:size(inLink,2)
        plot(inLink(l).comx,inLink(l).comy)
        text(inLink(l).comx(round(size(inLink(l).comx,2)/2)),inLink(l).comy(round(size(inLink(l).comx,2)/2)),num2str(l))
    end
    
    disp('Breakpoint') %Can trigger the debugger here
end

healNode = inNode;
healLink = inLink;