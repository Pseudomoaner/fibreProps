function [noteNode,noteLink] = annotateNetwork_old(inNode,inLink)
%ANNOTATENETWORK finds long fibres within the input network (input should
%be in graph format, output by skel2graph3d) based purely on the morphology
%of the network.
%
%   INPUTS:
%       -inNode: Input version of the nodes forming the network
%       -inLink: Input version of the links forming the network
%
%   OUTPUTS:
%       -noteNode: Output version of nodes, with nodes forming part of a
%       fibre annotated.
%       -noteLink: Output version of links, with links forming part of a
%       fibre annotated.
%
%   Author: Oliver J. Meacock, (c) 2020

%Analytical parameters
scoreThresh = 0.5;
curvW = 1;
angW = 0.4;
widW = 0.02;
angSampPts = 7; %Number of points away from node you should sample to determine each link's orientation (local to node)
curveSampPts = 10; %Number of points away from node you should sample to check local curvature
widthSampPts = angSampPts;

%Initialise empty fibre index storage
for n = 1:size(inNode,2)
    inNode(n).Fibres = [];
end
for l = 1:size(inLink,2)
    inLink(n).Fibres = [];
end

linkedLinks = double.empty(0,3);
for n = 1:size(inNode,2)
    %Calculate the local angle for each link emanating from this node
    currAngs = zeros(size(inNode(n).conn));
    for i = 1:size(inNode(n).conn,2)
        lkID = inNode(n).links(i);
        end1dist = sqrt((inLink(lkID).comx(1) - inNode(n).ptComx)^2 + (inLink(lkID).comy(1) - inNode(n).ptComy)^2);
        end2dist = sqrt((inLink(lkID).comx(end) - inNode(n).ptComx)^2 + (inLink(lkID).comy(end) - inNode(n).ptComy)^2);
        
        if end1dist < end2dist
            if numel(inLink(lkID).comx) > angSampPts
                xList = inLink(lkID).comx(1:angSampPts);
                yList = inLink(lkID).comy(1:angSampPts);
            else
                xList = inLink(lkID).comx;
                yList = inLink(lkID).comy;
            end
        else
            if numel(inLink(lkID).comx) > angSampPts
                xList = flip(inLink(lkID).comx(end-angSampPts+1:end));
                yList = flip(inLink(lkID).comy(end-angSampPts+1:end));
            else
                xList = flip(inLink(lkID).comx);
                yList = flip(inLink(lkID).comy);
            end
        end
        
        %Need to deal with 0 and Inf gradient cases
        if sum(yList == yList(1)) == numel(yList)
            if xList(end) > xList(1)
                currAngs(i) = 0;
            elseif xList(end) < xList(1)
                currAngs(i) = -pi;
            end
        elseif sum(xList == xList(1)) == numel(xList)
            if yList(end) > yList(1)
                currAngs(i) = pi/2;
            elseif yList(end) < yList(1)
                currAngs(i) = -pi/2;
            end
        else
            grad = (xList' - xList(1))\(yList' - yList(1));
            if mean(xList - xList(1)) > 0
                currAngs(i) = atan(grad);
            else
                currAngs(i) = wrapToPi(atan(grad)+pi);
            end
        end
    end
    inNode(n).angs = currAngs;
    
    %Now, for each pair of links eminating from this node, score that pair
    %of links. If score is below threshold, pair together as part of a
    %single fibre.
    IDs1 = [];
    IDs2 = [];
    scList = [];
    for l1Ind = 1:size(inNode(n).links,2)
        for l2Ind = (l1Ind+1):size(inNode(n).links,2)
            %Curvature of segment-segment join
            l1 = inNode(n).links(l1Ind);
            l2 = inNode(n).links(l2Ind);
            [xL1,yL1,~,flipFlag1] = matchNodeAndLinkEnd(inNode,inLink,n,l1);
            [xL2,yL2,~,flipFlag2] = matchNodeAndLinkEnd(inNode,inLink,n,l2);
            
            l1Max = min(numel(xL1),curveSampPts);
            l2Max = min(numel(xL2),curveSampPts);
            
            xList = [xL1(1:l1Max),xL2(1:l2Max)];
            yList = [yL1(1:l1Max),yL2(1:l2Max)];
            
            circMeasures = CircleFitByPratt([xList',yList']);
            curve = 1/circMeasures(3);
            
            %Angular distance between segment emission angles, minus pi radians (expectation)
            angDist = abs(wrapToPi(currAngs(l1Ind) - currAngs(l2Ind) - pi));
            
            %Difference in widths
            l1Max = min(numel(xL1),widthSampPts);
            if flipFlag1
                l1Wid = mean(inLink(l1).widths(end-l1Max+1:end));
            else
                l1Wid = mean(inLink(l1).widths(1:l1Max));
            end
            
            l2Max = min(numel(xL2),widthSampPts);
            if flipFlag2
                l2Wid = mean(inLink(l2).widths(end-l2Max+1:end));
            else
                l2Wid = mean(inLink(l2).widths(1:l2Max));
            end
            widDiff = abs(l1Wid - l2Wid);
            
            currSc = (curve*curvW + angDist*angW + widDiff*widW);
            if currSc < scoreThresh
                IDs1 = [IDs1;l1];
                IDs2 = [IDs2;l2];
                scList = [scList;currSc];
            end
        end
    end
    
    %As a special case, if you have detected three or four links belonging to the
    %same fibre eminating from this node, delete all but the two
    %highest-scoring pairings (this ensures crossing points aren't marked
    %as belonging to a single fibre)
    if numel(IDs1) > 2 && numel(inNode(n).conn) == 4
        [~,scOrd] = sort(scList);
        IDs1(scOrd(3:end)) = [];
        IDs2(scOrd(3:end)) = [];
    end
    
    linkedLinks = [linkedLinks;IDs1,ones(size(IDs1))*n,IDs2];
end

%Progress through the 'linkedLinks' array, marking off fibres as you go
fibInd = 1;
nodeStore = {};
linkStore = {};
while ~isempty(linkedLinks)
    currEnds = [linkedLinks(1,1);linkedLinks(1,3)];
    nodeStore{fibInd} = linkedLinks(1,2);
    linkStore{fibInd} = currEnds;
    
    linkedLinks(1,:) = [];
    
    while ~isempty(currEnds)
        repEnds = [];
        for i = 1:size(currEnds,1)
            endList1 = find(linkedLinks(:,1) == currEnds(i));
            endList2 = find(linkedLinks(:,3) == currEnds(i));
            
            newEnds = [linkedLinks(endList1,3);linkedLinks(endList2,1)];
            newNodes = [linkedLinks(endList1,2);linkedLinks(endList2,2)];
            
            repEnds = [repEnds;newEnds];
            nodeStore{fibInd} = [nodeStore{fibInd};newNodes];
            linkStore{fibInd} = [linkStore{fibInd};newEnds];
            
            linkedLinks(endList1,:) = [];
            linkedLinks(endList2,:) = [];
        end
        currEnds = repEnds;
    end
    fibInd = fibInd + 1;
end

for i = 1:size(nodeStore,2)
    for n = 1:size(nodeStore{i},1)
        inNode(nodeStore{i}(n)).Fibres = [inNode(nodeStore{i}(n)).Fibres;i];
    end
    
    for l = 1:size(linkStore{i},1)
        inLink(linkStore{i}(l)).Fibres = [inLink(linkStore{i}(l)).Fibres;i];
    end
end

noteNode = inNode;
noteLink = inLink;