function [noteNode,noteLink] = annotateNetwork(inNode,inLink,fibreGroups,dx)
%ANNOTATETWORK finds long fibres within the input network (input should
%be in graph format, output by skel2graph3d) based on the morphology
%of the network and the original labelling of the fibres by the ridge
%detection algorithm.
%
%   INPUTS:
%       -inNode: Input version of the nodes forming the network
%       -inLink: Input version of the links forming the network
%       -fibreGroups: Image containing labelled binary segmentations of the
%       fibres, as output by the ridge detection algorithm.
%       -dx: Spacing between pixels, in nm.
%
%   OUTPUTS:
%       -noteNode: Output version of nodes, with nodes forming part of a
%       fibre annotated.
%       -noteLink: Output version of links, with links forming part of a
%       fibre annotated.
%
%   Author: Oliver J. Meacock, (c) 2020

%Analytical parameters
angSampPts = round(2/dx); %Number of points away from node you should sample to determine each link's orientation (local to node)
angThreshHi = deg2rad(60); %A pair of links coming out of a single node must be within pi +- angThreshHi radians of each other not to be cut
angThreshLo = deg2rad(40); %A pair of links coming out of a single node must be within pi +- angThreshLo radians of each other to be fused into a single fibre

%% Part 1: assign angles to each link eminating from each node
for n = 1:size(inNode,2)
    %Calculate the local angle for each link emanating from this node
    currAngs = zeros(size(inNode(n).conn));
    for i = 1:size(inNode(n).conn,2)
        lkID = inNode(n).links(i);
        end1dist = sqrt((inLink(lkID).comx(1) - inNode(n).ptComx)^2 + (inLink(lkID).comy(1) - inNode(n).ptComy)^2);
        end2dist = sqrt((inLink(lkID).comx(end) - inNode(n).ptComx)^2 + (inLink(lkID).comy(end) - inNode(n).ptComy)^2);
        
        %New version of angle detection code just uses a single point the specified distance away from the target node 
        if end1dist < end2dist
            if numel(inLink(lkID).comx) > angSampPts
                xList = inLink(lkID).comx([1,angSampPts]);
                yList = inLink(lkID).comy([1,angSampPts]);
            else
                xList = inLink(lkID).comx([1,end]);
                yList = inLink(lkID).comy([1,end]);
            end
        else
            if numel(inLink(lkID).comx) > angSampPts
                xList = inLink(lkID).comx([end,end-angSampPts+1]);
                yList = inLink(lkID).comy([end,end-angSampPts+1]);
            else
                xList = inLink(lkID).comx([end,1]);
                yList = inLink(lkID).comy([end,1]);
            end
        end
%         if end1dist < end2dist
%             if numel(inLink(lkID).comx) > angSampPts
%                 xList = inLink(lkID).comx(1:angSampPts);
%                 yList = inLink(lkID).comy(1:angSampPts);
%             else
%                 xList = inLink(lkID).comx;
%                 yList = inLink(lkID).comy;
%             end
%         else
%             if numel(inLink(lkID).comx) > angSampPts
%                 xList = flip(inLink(lkID).comx(end-angSampPts+1:end));
%                 yList = flip(inLink(lkID).comy(end-angSampPts+1:end));
%             else
%                 xList = flip(inLink(lkID).comx);
%                 yList = flip(inLink(lkID).comy);
%             end
%         end
        
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
            grad = (yList(2) - yList(1))/(xList(2) - xList(1));
            if xList(2) - xList(1) > 0
                currAngs(i) = atan(grad);
            else
                currAngs(i) = wrapToPi(atan(grad)+pi);
            end
        end
    end
    inNode(n).angs = currAngs;
end

%% Part 2: Assign links and nodes to fibres according to the fibreGroups image
%Initialise empty fibre index storage
for n = 1:size(inNode,2)
    inNode(n).Fibres = [];
end

%Each link can be associated with only a single fibre...
fibreInds = bwlabel(fibreGroups,4)';
miniFibInd = max(fibreInds(:)) + 1; %Ensures mini-fibre indices are all greater than the main fibre indices
for l = 1:size(inLink,2)
    indLocs = sub2ind(size(fibreGroups),inLink(l).comx,inLink(l).comy);
    if numel(indLocs) > 5 %If link is at least 5 pixels long
        fibIndList = fibreInds(indLocs); %The indicies associated with this 
        fibIndList(fibIndList == 0) = [];
        
        if sum(fibIndList == mode(fibIndList)) > numel(fibIndList)/2 %If at least half of the points in this link are associated with a fibre in the original ridge-detected image...
            inLink(l).Fibre = mode(fibIndList);

        else %Associate link with a new 'mini fibre', unconnected to any other links or nodes.
            inLink(l).Fibre = miniFibInd;
            miniFibInd = miniFibInd+1; 
        end
    else %Associate link with a new 'mini fibre', unconnected to any other links or nodes.
        inLink(l).Fibre = miniFibInd;
        miniFibInd = miniFibInd+1;
    end
end

%Nodes can be associated with multiple fibres
for n = 1:size(inNode,2)
    fibList = double.empty(0,1);
    for el = 1:size(inNode(n).links,2)
        fibList = [fibList,inLink(inNode(n).links(el)).Fibre];
    end
    inNode(n).Fibres = fibList;
end

%% Part 3: Break fibre assignments where more than two links at a single node belong to a single fibre
%Find maximal fibre index
maxFib = max([inNode.Fibres]);

for n = 1:size(inNode,2)
    fibCnts = histcounts(inNode(n).Fibres,0.5:maxFib+0.5);
    badInds = find(fibCnts > 2); %Find the indicies of any fibres coming out of this node that are dodgy (more than two links associated with them)
    %Go through each dodgy index...
    for bi = badInds
        tgtLks = inNode(n).Fibres == bi;
        tgtLksInds = find(tgtLks);
        angList = inNode(n).angs(tgtLks);
        
        smallestDiff = 2*pi;
        for a1 = 1:size(angList,2)
            for a2 = 1:size(angList,2)
                if wrapToPi(circ_dist(angList(a1),angList(a2))-pi) < smallestDiff && a1~=a2
                    smallestDiff = wrapToPi(circ_dist(angList(a1),angList(a2))-pi);
                    keepi = a1;
                    keepj = a2;
                end
            end
        end
        
        ridInds = true(size(angList));
        ridInds(keepi) = false;
        ridInds(keepj) = false;
        ridInds = tgtLksInds(ridInds);
        ridLinks = inNode(n).links(ridInds); %Target links to disconnect
        ridNodes = inNode(n).conn(ridInds);
        
        %Now reindex everything connected to the links you've found you want to disconnect...
        for ri = 1:size(ridLinks,2)
            [chngNds,chngLks] = findSubfibreIndices(inNode,inLink,ridNodes(ri),ridLinks(ri),20);
            
            maxFib = maxFib + 1;
            for cn = 1:size(chngNds,2)
                inLink(chngLks(cn)).Fibre = maxFib;
                inNode(chngNds(cn)).Fibres(inNode(chngNds(cn)).Fibres == bi) = maxFib;
            end
            inNode(n).Fibres(ridInds(ri)) = maxFib; %Finally update the origin node
        end
    end
end

%% Part 4: Break fibre assignments where link pairs are more than angThreshHi away from each other at a given node
for n = 1:size(inNode,2)
    fibCnts = histcounts(inNode(n).Fibres,0.5:maxFib+0.5);
    queryInds = find(fibCnts == 2); %Find the indicies of any fibres coming out of this node that are of interest (exactly two links associated with them)
    %Go through each query index...
    for qi = queryInds
        tgtLks = inNode(n).Fibres == qi;
        angList = inNode(n).angs(tgtLks);
        angDiff = wrapToPi(circ_dist(angList(1),angList(2))+pi);
        
        if abs(angDiff) > angThreshHi %If the difference between the link angles is too great, break the assignment to the same fibre
            tgtInds = find(tgtLks);
            ridLink = inNode(n).links(tgtInds(1));
            ridNode = inNode(n).conn(tgtInds(1));
            
            [chngNds,chngLks] = findSubfibreIndices(inNode,inLink,ridNode,ridLink,20);
            
            maxFib = maxFib + 1;
            for cn = 1:size(chngNds,2)
                inLink(chngLks(cn)).Fibre = maxFib;
                inNode(chngNds(cn)).Fibres(inNode(chngNds(cn)).Fibres == qi) = maxFib;
            end
            inNode(n).Fibres(tgtInds(1)) = maxFib; %Finally update the origin node
        end
    end
end

%% Part 5: Break apart any loop structures in the fibres
for f = 1:maxFib
    linkList = [];
    nodeListA = [];
    nodeListB = [];
    for l = 1:size(inLink,2)
        if inLink(l).Fibre == f
            linkList = [linkList;l];
            nodeListA = [nodeListA;inLink(l).n1];
            nodeListB = [nodeListB;inLink(l).n2];
        end
    end
    
    %If fibre is part of one (or several) loops, the number of unique nodes
    %should be less than or equal to the number of links
    noLinks = numel(linkList);
    noNodes = numel(unique([nodeListA;nodeListB]));
    
    if noLinks > 0 && noNodes > 0 && noLinks >= noNodes %If you've found a loop, break it apart into separate link
        col = rand(1,3);
        for l = 1:size(linkList,1)
            maxFib = maxFib + 1;
            lk = linkList(l);
            nd1 = nodeListA(l);
            nd2 = nodeListB(l);
            
            inNode(nd1).Fibres(inNode(nd1).links == lk) = maxFib;
            inNode(nd2).Fibres(inNode(nd2).links == lk) = maxFib;
            inLink(lk).Fibre = maxFib;
        end
    end
end

%% Part 6: Fuse fibres that terminate at a common node and links are less than angThreshLo away from each other at that node
for n = 1:size(inNode,2)
    angList = inNode(n).angs';
    angDists = abs(triu(squareform(wrapToPi(pdist(angList,@circ_dist)+pi)),1));
    angDists(angDists == 0) = nan;
    
    [smallestDiff,minLoc] = min(angDists(:));
    while smallestDiff < angThreshLo
        [linkInd1,linkInd2] = ind2sub(size(angDists),minLoc);
        angDists([linkInd1,linkInd2],:) = nan;
        angDists(:,[linkInd1,linkInd2]) = nan;
        
        [smallestDiff,minLoc] = min(angDists(:));
        
        %Combine detected links
        joinLink = inNode(n).links(linkInd2);
        joinNode = inNode(n).conn(linkInd2);
        extendFib = inNode(n).Fibres(linkInd1);
        ridFib = inNode(n).Fibres(linkInd2);
        
        if extendFib ~= ridFib && sum(inNode(n).Fibres == extendFib) == 1 %Second check ensures you don't proceed if fibre already enters and exits node
            [chngNds,chngLks] = findSubfibreIndices(inNode,inLink,joinNode,joinLink,20);
            
            for cn = 1:size(chngNds,2)
                inLink(chngLks(cn)).Fibre = extendFib;
                inNode(chngNds(cn)).Fibres(inNode(chngNds(cn)).Fibres == ridFib) = extendFib;
            end
            inNode(n).Fibres(linkInd2) = extendFib; %Finally update the origin node
        end
    end
end

noteNode = inNode;
noteLink = inLink;