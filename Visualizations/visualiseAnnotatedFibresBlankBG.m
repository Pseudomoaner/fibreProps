function showIm = visualiseAnnotatedFibresBlankBG(fibreProps,origZ,colType,dx)
%VISUALISEANNOTATEDFIBRESBLANKBG draws and displays a reconstruction of
%automatically detected and measured fibres in an AFM image. In contrast to
%visualiseAnnotatedFibres, the original image is not used as a background.
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

%% Draw fibres

%Set background as 
rCh = ones(size(origZ));
gCh = ones(size(origZ));
bCh = ones(size(origZ));

switch colType
    case {'localOrientation','meanOrientation'}
        cmap = colormap('hsv'); %For angular variables
    case 'Length'
        cmap = colormap('jet'); %For linear variables
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
            if ~isnan(sortFibScores(F)) %NaN for any fibres that weren't associated with links
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
        fibScores = [fibreProps(:).meanOrientation];

        for F = 1:size(fibreProps,2)
            if ~isnan(fibScores(F))
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
        [xGrid,yGrid] = meshgrid(1:size(origZ,1),1:size(origZ,2));
        for F = 1:size(fibreProps,2)
            if sum(isnan(fibreProps(F).localOrientation)) == 0
                for i = 1:size(fibreProps(F).localOrientation,1)
                    currCInd = ceil(((fibreProps(F).localOrientation(i)+pi/2)/pi)*size(cmap,1));
                    currCInd = min(currCInd,size(cmap,1));
                    currCInd = max(currCInd,1);

                    cVals = cmap(currCInd,:);
                    
                    xLoc = fibreProps(F).backList(i,1);
                    yLoc = fibreProps(F).backList(i,2);

                    localWidth = round(fibreProps(F).width/widReconFac);
                    currInds = sqrt((xGrid-xLoc).^2 + (yGrid-yLoc).^2) < localWidth; 

                    rCh(currInds) = rCh(currInds)/2 + cVals(1)/2;
                    gCh(currInds) = gCh(currInds)/2 + cVals(2)/2;
                    bCh(currInds) = bCh(currInds)/2 + cVals(3)/2;
                end
            end
            disp(['F is ',num2str(F),' of ',num2str(size(fibreProps,2)),'.'])
        end
end

showIm = cat(3,rCh,gCh,bCh);
