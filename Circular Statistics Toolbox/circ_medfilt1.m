function outData = circ_medfilt1(inData,span)
%CIRC_MEDFILT1 Performs similar function to medfilt1, but acts on circular data instead of linear data.
%
%   inData = angular data (specified in radians) for median filtering
%
%   span = span of moving window (centered on current location)
%
%   outData = returned moving median value of data.
if rem(span,2) ~= 1
    error('Span should be an odd number!')
end

spanFor = (span - 1)/2;
spanBack = (span - 1)/2;

outData = zeros(size(inData));

for i = 1:length(inData)
    if i <= spanBack
        startInd = 1;
        endInd = span;
    elseif i > length(inData) - spanFor
        startInd = length(inData) - span;
        endInd = length(inData);
    else
        startInd = i - spanBack;
        endInd = i + spanFor;
    end


    outData(i) = circ_median(inData(startInd:endInd));
    
    if outData(i) > pi
        inData(startInd:endInd)
    end
end