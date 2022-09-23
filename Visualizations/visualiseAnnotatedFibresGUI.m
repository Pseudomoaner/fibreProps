function [] = visualiseAnnotatedFibresGUI(fibreProps,noteNodes,noteLinks,origZ,colType,dx,ax)
%VISUALISEANNOTATEDFIBRES draws and displays a reconstruction of
%automatically detected and measured fibres in an AFM image.
%
%   INPUTS:
%       -fibreProps: Structure containing the detected fibres and
%       associated properties, as output by measureFibres.m.
%       -noteNods: Annotated node structure, as output by
%       annotateNetwork.m.
%       -origZ: Either the original AFM heightmap, or the flattened version
%       of the same.
%       -colType: String specifying the type of colourscheme you want for
%       the reconstruction. Either 'localOrientation' (shows the local
%       orientation at each position in each fibre), 'globalOrientation'
%       (shows the global average orientation of each fibre) or 'Length
%       (shows the length of each fibre).
%
%   Author: Oliver J. Meacock, 2022

widReconFac = dx*4; %Will display fibres at half measured width.

hold(ax,'on')

%% Draw fibres

imgLow = prctile(origZ(:),1);
imgHi = prctile(origZ(:),99);
sharpImg = (origZ-imgLow)/(imgHi-imgLow);
sharpImg(sharpImg > 1) = 1; sharpImg(sharpImg < 0) = 0;

%Paste backprojection into each channel
rCh = sharpImg;
gCh = sharpImg;
bCh = sharpImg;

switch colType
    case {'localOrientation','globalOrientation'}
        cmap = colormap(ax,'hsv'); %For angular variables
    case 'Length'
        cmap = colormap(ax,'turbo'); %For linear variables
end

switch colType
    case 'Length'
        showFracLo = 0; %What fraction of the largest fibres should be shown
        showFracHi = 1;
        %Sort fibres by a value associated with them (length, width etc.)
        fibScores = zeros(size(fibreProps));
        for F = 1:size(fibreProps,2)
            fibScores(F) = size(fibreProps(F).backList,1);
        end
        [sortFibScores,order] = sort(fibScores);

        for F = round(size(sortFibScores,2)*showFracLo)+1:round(size(sortFibScores,2)*showFracHi) %Loop through fibres
            if sortFibScores(F) > 1 %1 for any fibres that weren't associated with links
                currCInd = ceil(((sortFibScores(F)-min(sortFibScores))/(max(sortFibScores)-min(sortFibScores)))*size(cmap,1));
                currCInd = min(currCInd,size(cmap,1));
                currCInd = max(currCInd,1);

                cVals = cmap(currCInd,:);
                Find = order(F);
                
                backboneImg = zeros(size(origZ));
                backboneImg(sub2ind(size(origZ),fibreProps(Find).backList(:,1),fibreProps(Find).backList(:,2))) = 1;

                se = strel('disk',round(fibreProps(Find).width/widReconFac));
                currInds = logical(imdilate(backboneImg',se));

                rCh(currInds) = rCh(currInds)/2 + cVals(1)/2;
                gCh(currInds) = gCh(currInds)/2 + cVals(2)/2;
                bCh(currInds) = bCh(currInds)/2 + cVals(3)/2;
            end
        end
    case 'meanOrientation'
        %Concateanate fibre scores
        fibScores = [fibreProps(F).meanOrientation];

        for F = 1:size(fibreProps,2)
            if fibScores(F) > 1
                currCInd = ceil(((fibreProps(F).meanOrientation(i)+pi/2)/pi)*size(cmap,1));
                currCInd = min(currCInd,size(cmap,1));
                currCInd = max(currCInd,1);

                cVals = cmap(currCInd,:);

                backboneImg = zeros(size(origZ));
                backboneImg(sub2ind(size(origZ),fibreProps(F).backList(:,1),fibreProps(F).backList(:,2))) = 1;

                se = strel('disk',round(fibreProps(F).width/widReconFac));
                currInds = logical(imdilate(backboneImg',se));

                rCh(currInds) = rCh(currInds)/2 + cVals(1)/2;
                gCh(currInds) = gCh(currInds)/2 + cVals(2)/2;
                bCh(currInds) = bCh(currInds)/2 + cVals(3)/2;
            end
        end
    case 'localOrientation' %Rather expansive for what it does... may be able to improve
        for F = 1:size(fibreProps,2)
            if sum(isnan(fibreProps(F).localOrientation)) == 0
                se = strel('disk',round(fibreProps(F).width/widReconFac));
                for i = 1:size(fibreProps(F).localOrientation,1)
                    currCInd = ceil(((fibreProps(F).localOrientation(i)+pi/2)/pi)*size(cmap,1));
                    currCInd = min(currCInd,size(cmap,1));
                    currCInd = max(currCInd,1);

                    cVals = cmap(currCInd,:);
                    
                    backboneImg = zeros(size(origZ));
                    backboneImg(sub2ind(size(origZ),fibreProps(F).backList(:,1),fibreProps(F).backList(:,2))) = 1;
                    currInds = logical(imdilate(backboneImg,se));

                    rCh(currInds) = rCh(currInds)/2 + cVals(1)/2;
                    gCh(currInds) = gCh(currInds)/2 + cVals(2)/2;
                    bCh(currInds) = bCh(currInds)/2 + cVals(3)/2;
                end
            end
        end
end

showIm = cat(3,rCh,gCh,bCh);
imshow(showIm,'Parent',ax)

%% Draw fibre-fibre crossing points

for n = 1:size(noteNodes,2)
    if numel(noteNodes(n).links) == 4
        skip = false;
        %Ensure that this node isn't connected to any (nearby) terminal nodes
        for l = 1:4
            if numel(noteNodes(noteNodes(n).conn(l)).conn) == 1 && numel(noteLinks(noteNodes(n).links(l)).point) < 15
                skip = true;
            end
        end
        
        fibList = zeros(4,1);
        for l = 1:4
            %Also ensure that these are two annotated fibres that are crossing
            if numel(noteLinks(noteNodes(n).links(l)).Fibre) ~= 1
                skip = true;
            else
                fibList(l) = noteLinks(noteNodes(n).links(l)).Fibre;
            end            
        end

        if numel(unique(fibList)) ~= 2
            skip = true;
        end
        
        if ~skip
            plot(ax,noteNodes(n).ptComx,noteNodes(n).ptComy,'o','MarkerFaceColor','k','MarkerEdgeColor','w')
        end
    end
end

%% Display fibre indices

% for F = 1:size(fibreProps,2)
%     text(ax,fibreProps(F).midpoint(1),fibreProps(F).midpoint(2),['F',num2str(F)])
% end
% 
% title('Back-projected fibres')