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
        
        midInd = sub2ind(size(origImg),fibreProps(measInd).midpoint(1),fibreProps(measInd).midpoint(2));
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
        fibreProps(measInd).midpoint = nan(1,2);
        fibreProps(measInd).midwidth = nan(1);
    end
    
    %Find the fibre orientation
    [fibrePxY,fibrePxX] = ind2sub(size(origImg),find(fibrePx));
    fibreProps(measInd).orientation = atan((fibrePxX-mean(fibrePxX))\(fibrePxY-mean(fibrePxY)));
    
    fibreProps(measInd).size = sum(fibrePx(:))*dx; %Length of the fibre
    fibreProps(measInd).backbone = fibrePx;
    fibreProps(measInd).score = mean(scoreSet);
    fibreProps(measInd).width = mean(widthSet)*dx;
    fibreProps(measInd).branchNo = branchNo;
    fibreProps(measInd).rawInd = F;
    
    measInd = measInd + 1;
end