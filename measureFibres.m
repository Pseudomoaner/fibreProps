function fibreProps = measureFibres(inNodes,inLinks,origImg,dx)
%MEASUREFIBRES measures the properties of fibres detected in the given
%input image.
%
%   INPUTS:
%       -inNodes: The fully processed and labelled version of the
%       information about the network graph's nodes.
%       -inLinks: The fully processed and labelled version of the
%       information about the network graph's nodes.
%       -origImg: The original AFM image, as output by txtToMat.
%       -dx: The spacing between adjacent pixels, in nm.
%
%   OUTPUTS:
%       -fibreProps: Measured properties of the input fibres. Structure
%       consists of the following fields:
%           -orientation: Angle made by fibre relative to image reference 
%           frame (in radians, -pi/2 to pi/2)
%           =midPoint: vector specifying the midpoint of the fibre
%           -midWidth: Measured width of the fibre at the midpoint position
%           (in pixels)
%           -size: Length of the fibre (in nm)
%           -backbone: Binary image, specifying the single pixel-wide
%           'backbone' of the fibre.
%           -score: Average score assigned to the entire collection of
%           pixels associated with the fibre by the ridge detection
%           algorithm.
%           -width: Average width of the fibre (in nm)
%           -branchNo: Number of branches in the fibre. In current version
%           of code, should always be 0.
%           -rawInd: Index of the fibre as stored within the graph
%           structures.
%
%   Author: Oliver J. Meacock

%Find the list of fibre indices that are actually represented in this
%dataset
FInds = [];
for l = 1:size(inLinks,2)
    FInds = [FInds;inLinks(l).Fibre];
end
FInds = unique(FInds);

fibreProps = struct();

measInd = 1; %Index used to index into the fibreProps structure
for F = FInds'
    %Reconstruct fibre
    fibrePx = zeros(size(origImg));
    scoreSet = [];
    widthSet = [];
    branchNo = 0;
    for l = 1:size(inLinks,2)
        if sum(inLinks(l).Fibre == F) > 0
            fibrePx(inLinks(l).point) = 1;
            scoreSet = [scoreSet,inLinks(l).Nscore];
            widthSet = [widthSet,inLinks(l).widths];
        end
    end
    for n = 1:size(inNodes,2)
        if sum(inNodes(n).Fibres == F) > 0
            fibrePx(inNodes(n).idx) = 1;
            scoreSet = [scoreSet,inNodes(n).Nscore'];
            widthSet = [widthSet,inNodes(n).widths'];
            
            %Check the number of links coming out of this node that belong
            %to this fibre. If more than 2, implies this is a branch point
            fibCnt = 0;
            for l = 1:size(inNodes(n).links,2)
                if sum(inLinks(inNodes(n).links(l)).Fibre == F) > 0
                    fibCnt = fibCnt + 1;
                end
            end
            branchNo = branchNo + max(fibCnt - 2,0);
        end
    end
    
    fibrePx = bwmorph(fibrePx,'thin',inf);
    fibrePx = bwmorph(fibrePx,'spur');
    
    %Find the mid-point of the fibre (defined as the half-way point between
    %its ends, so only properly defined for branchless fibres)
    ends = bwmorph(fibrePx,'endpoints');
    startpoint = find(ends(:) == 1,1);
    if branchNo == 0 && ~isempty(startpoint) %Fibre cannot be branched or loop-shaped
        [startY,startX] = ind2sub(size(origImg),startpoint);
        bound = bwtraceboundary(fibrePx,[startY,startX],'S');
        if size(bound,1) < 10
            bound = bwtraceboundary(fibrePx,[startY,startX],'N');
        end
        fibreProps(measInd).midpoint = bound(round(size(bound,1)/4),:);
        fibreProps(measInd).backList = bound(1:floor(size(bound,1)/2),:);

        midInd = sub2ind(size(origImg),fibreProps(measInd).midpoint(1),fibreProps(measInd).midpoint(2));

        %Search for the width of the fibre at its midpoint in the links
        %structure
        foundFlag = false;
        for l = 1:size(inLinks,2)
            if sum(inLinks(l).point == midInd) == 1
                fibreProps(measInd).midwidth = inLinks(l).widths(logical(inLinks(l).point == midInd));
                foundFlag = true;
                break
            end
        end
        
        if ~foundFlag
            for n = 1:size(inNodes,2)
                if sum(inNodes(n).idx == midInd) == 1
                    fibreProps(measInd).midwidth = inNodes(n).widths(inNodes(n).idx == midInd);
                    foundFlag = true;
                    break
                end
            end
        end
        
        if ~foundFlag
            fibreProps(measInd).midwidth = nan(1); %Couldn't find the midpoint for some reason
        end
    else
        fibreProps(measInd).backList = nan(1);
        fibreProps(measInd).midpoint = nan(1,2);
        fibreProps(measInd).midwidth = nan(1);
    end
    
    %Find the fibre orientation
    [fibrePxY,fibrePxX] = ind2sub(size(origImg),find(fibrePx));
    if abs(fibrePxX(end) - fibrePxX(1)) < abs(fibrePxY(end) - fibrePxY(1)) %Least-squares fitting relies on having lots of independant x-values - if these are not available, switch x and y coordinates and add pi/2 to result
        fibreProps(measInd).meanOrientation = atan((fibrePxY-mean(fibrePxY))\(fibrePxX-mean(fibrePxX)));
    else
        unwrapped = atan((fibrePxX-mean(fibrePxX))\(fibrePxY-mean(fibrePxY)));
        fibreProps(measInd).meanOrientation = -wrapToPi(unwrapped*2 + pi)/2;
    end
%     fitLine = polyfit(fibrePxX,fibrePxY,1);
%     fibreProps(measInd).meanOrientation = atan(fitLine(1));

    %Half-width of moving window used to calculate local fibre orientatino
    oriWindHalfWidth = 10;
    if ~isnan(fibreProps(measInd).backList)
        fibreProps(measInd).localOrientation = zeros(size(fibreProps(measInd).backList,1),1);
        for i = 1:size(fibreProps(measInd).backList,1)
            ind1 = max(1,i-oriWindHalfWidth);
            ind2 = min(size(fibreProps(measInd).backList,1),i+oriWindHalfWidth);

            pxList = fibreProps(measInd).backList(ind1:ind2,:);
            if abs(pxList(end,1) - pxList(1,1)) < abs(pxList(end,2) - pxList(1,2)) %Least-squares fitting relies on having lots of independant x-values - if these are not available, switch x and y coordinates and add pi/2 to result
                fibreProps(measInd).localOrientation(i) = atan((pxList(:,2)-mean(pxList(:,2)))\(pxList(:,1)-mean(pxList(:,1))));
            else
                unwrapped = atan((pxList(:,1)-mean(pxList(:,1)))\(pxList(:,2)-mean(pxList(:,2))));
                fibreProps(measInd).localOrientation(i) = -wrapToPi(unwrapped*2 + pi)/2;
            end
%             fitLine = polyfit(pxList(:,1),pxList(:,2),1);
%             fibreProps(measInd).localOrientation(i) = atan(1/fitLine(1));
        end
    else
        fibreProps.localOrientation = NaN;
    end
    
    fibreProps(measInd).size = sum(fibrePx(:))*dx; %Length of the fibre
    fibreProps(measInd).backbone = fibrePx;
    fibreProps(measInd).score = mean(scoreSet);
    fibreProps(measInd).width = mean(widthSet)*dx;
    fibreProps(measInd).branchNo = branchNo;
    fibreProps(measInd).rawInd = F;
    
    measInd = measInd + 1;
end