function centCoord = findCircleCentre(links,maxSize)
%FINDCIRCLECENTRE detects the centre of a series of concentric circles in
%the given fibre map. Will attempt to do so automatically before passing
%over control to the user if this fails.
%
%   INPUTS:
%       links: The fully processed and labelled version of the
%       information about the network graph's edges.
%       maxSize: The maximum size of the image in both dimensions (in pixels)
%
%   OUTPUTS:
%       centCoord: The coordinates of the centre of the circle (in pixels).

%Initial attempt below based on entropy minimisation of orientation
%distribution

sampPts = 50;
xList = linspace(maxSize(2)/(sampPts*2),maxSize(2) - maxSize(2)/(sampPts*2),50);
yList = linspace(maxSize(1)/(sampPts*2),maxSize(1) - maxSize(1)/(sampPts*2),50);
scoreGrid = zeros(sampPts);

%Create lists of coordinates for fibre segments
orientList = zeros(size(links,2),1);
centList = zeros(size(links,2),2);
for k = 1:size(links,2)
    comx = links(k).comx';
    comy = links(k).comy';
    if abs(comx(end) - comx(1)) < abs(comy(end) - comy(1)) %Least-squares fitting relies on having lots of independant x-values - if these are not available, switch x and y coordinates and add pi/2 to result
        meanOrient = -wrapToPi(atan((comy-mean(comy))\(comx-mean(comx)))*2 + pi)/2;
    else
        unwrapped = atan((comx-mean(comx))\(comy-mean(comy)));
        meanOrient = -wrapToPi((-wrapToPi(unwrapped*2)/2)*2)/2;
    end
    orientList(k) = meanOrient;
    centList(k,1) = mean(comx);
    centList(k,2) = mean(comy);
end

%For each lattice site, calculate orientation scores
for i = 1:sampPts
    for j = 1:sampPts
        ptX = xList(j);
        ptY = yList(i);
        reorientList = zeros(size(links,2),1);
        for k = 1:size(links,2)
            corrAng = atan2(centList(k,2)-ptY,centList(k,1)-ptX);
            reorientList(k) = orientList(k)-corrAng;
        end
        reorientList = wrapToPi(reorientList*2)/2;
        %Calculate the entropy of the resulting orientation distribution
        [counts,binCenters] = hist(reorientList,20);
        binWidth = diff(binCenters);
        binWidth = [binWidth(end),binWidth]; % Replicate last bin width for first, which is indeterminate.
        nz = counts>0; % Index to non-zero bins
        frequency = counts(nz)/sum(counts(nz));
        scoreGrid(i,j) = -sum(frequency.*log(frequency./binWidth(nz)));
    end
end

[minEnt,centInd] = min(scoreGrid,[],'all');

%If minimised entropy is low enough, perform automated centre detection.
%Otherwise, ask for it from user.
if minEnt < 1.11
    [centY,centX] = ind2sub([sampPts,sampPts],centInd);
    centCoord = [xList(centX),yList(centY)];
else
    warn = warndlg('Could not find centre of circular structure automatically. Please select manually.','Manual input required');
    uiwait(warn)
    centCoord = NaN;
end

%Approach below based on the tensor method for orientation analysis
 
%Keep stepping up filter size until only one pair of +1/2 defects remains
% [gX,gY] = imgradientxy(flatImg);
% Ixx = gX.*gX;
% Iyy = gY.*gY;
% Ixy = gX.*gY;
% 
% tIxx = imgaussfilt(Ixx,tensorScale,'FilterSize',2*ceil(tensorScale*8)+1); %Extra filter size should suffice to ensure smoothness of image gradients (i.e. no zero values)
% tIyy = imgaussfilt(Iyy,tensorScale,'FilterSize',2*ceil(tensorScale*8)+1);
% tIxy = imgaussfilt(Ixy,tensorScale,'FilterSize',2*ceil(tensorScale*8)+1);
% 
% orients = 0.5 * atan2(2*tIxy,tIyy-tIxx);
