function [chngNds,chngLks] = findSubfibreIndices(inNode,inLink,ridNode,ridLink,cnt)
%FINDSUBFIBREINDICES uses a recursive method to find all the links and
%nodes connected as part of the same fiber as the target link and node.
%Note the term cnt is provided to break recursion loops.
%
%   INPUTS:
%       -inNode: The structure containing information about the nodes in
%       the AFM network reconstruction graph.
%       -inLink: The structure containing informatino about the links in
%       the AFM network reconstruction graph.
%       -ridNode: The current target node.
%       -ridLink: The current target link. This link will be ignored when
%       looking for links exiting ridNode belonging to the same fibre.
%       -cnt: Countdown variable that causes code to terminate if it
%       reaches zero. Generally indicates the presence of a loop in your
%       assigned fibre.
%
%   OUTPUTS:
%       -chngNds: Indices of all nodes associated with this fibre that are
%       'downstream' of ridNode and ridLink.
%       -chngLks: Indices of all links associated with this fibre that are
%       'downstream' of ridNode and ridLink.
%
%   Author: Oliver J. Meacock, (c) 2021

%Targetted link contains the target fibre index
tgtFib = inLink(ridLink).Fibre;

%Find all the links out of the target node belonging to that fibre
%(discouting ridLink)
tgtInds = inNode(ridNode).Fibres == tgtFib;
tgtInds(inNode(ridNode).links == ridLink) = false;
tgtInds = find(tgtInds);

%If tgtInds is empty, must be at fibre end. Return just these indices.
chngNds = ridNode;
chngLks = ridLink;

%Otherwise, compile the list of links and nodes emanating from the
%newly targetted links and nodes
if ~isempty(tgtInds) && cnt > 0
    for i = 1:size(tgtInds,2)
        [chngNdsSub,chngLksSub] = findSubfibreIndices(inNode,inLink,inNode(ridNode).conn(tgtInds(i)),inNode(ridNode).links(tgtInds(i)),cnt-1);
        chngNds = [chngNds,chngNdsSub];
        chngLks = [chngLks,chngLksSub];
    end
end

%Remove any duplicates that may have emerged in recursion loops
chngNds = unique(chngNds);
chngLks = unique(chngLks);