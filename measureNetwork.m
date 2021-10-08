function networkMeasures = measureNetwork(noteNodes,noteLinks,fibreProps,origZvals,dx,widFac)
%MEASURENETWORK measures summary statistics about the given AFM network
%reconstruction.
%
%   INPUTS:
%       -noteNodes: Structure containing information about the graph's
%       nodes, with fibres annotated by the annotateNetwork function.
%       -noteLinks: Structure containing information about the graph's
%       links, with fibres annotated by the annotateNetwork function.
%       -fibreProps: Measured properties of each fibre, as output by the
%       measureFibres function.
%       -origZ: Original heightmap AFM image.
%       -dx: Distance (in nm) between adjacent pixels in image.
%       -widFac: Callibration factor that compensates for the difference
%       between detected and actual fibre widths. If detected width is X,
%       real width is measured as X/widFac.
%
%   OUTPUTS:
%       -networkMeasures: Structure containing the following fields:
%           -order: Scalar order parameter of the network
%           -fibWidMu: Average width of fibres in the network
%           -fibWidCOV: Coefficient of variation of fibre width
%           -fibLenMu: Average length of fibres in the network
%           -fibLenCOV: Coefficient of variation of fibre length
%           -netDensity: Fraction of total area of image that is filled
%           with links.
%           -nodeDensity: Average number of nodes per unit area.
%
%   Author: Oliver J. Meacock, (c) 2021

sizePrctile = 80;
fibSizes = [fibreProps.size];
threshSize = prctile(fibSizes,sizePrctile);
badFibs = fibSizes < threshSize;
fibreProps(badFibs) = [];

%Will need to scale by dimensional units eventually, but there seems to be
%an issue with the way they were saved for some of these data...

%Kappa of fibre orientation distribution
meanOri = circ_mean([fibreProps.orientation]');
networkMeasures.order = mean(2*cos([fibreProps.orientation]-meanOri).^2 - 1);

%Mean and COV of fibre width distribution
networkMeasures.fibWidMu = mean([fibreProps.width])*dx/widFac;
networkMeasures.fibWidCOV = std([fibreProps.width])/mean([fibreProps.width]); %Non-dimensional

%Mean and COV of fibre length distribution
networkMeasures.fibLenMu = mean([fibreProps.size])*dx;
networkMeasures.fibLenCOV = std([fibreProps.size])/mean([fibreProps.size]); %Non-dimensional

%Network density (total fraction of image covered by binary links)
linkBin = zeros(size(origZvals));
for l = 1:size(noteLinks,2)
    linkBin(noteLinks(l).point) = 1;
end

networkMeasures.netDensity = sum(linkBin(:))/(dx*dx*numel(linkBin));

%Node density (number of nodes divided by image area
networkMeasures.nodeDensity = size(noteNodes,2)/(dx*dx*numel(linkBin));