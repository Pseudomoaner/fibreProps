function backProjImg = visualiseAnnotatedNetwork(noteNodes,noteLinks,fibreProps,origZ)
%VISUALISEANNOTAEDNETWORK allows you to plot the output of the network
%reconstruction and fibre assignment algorithms to visually ensure
%accuracy. Plots are the following:
%   Fig1, left: Plot of the graph-based reconstruction of the fibre
%   network, overlaid on top of the original image. Contains black lines
%   indicating the links, and cyan circles indicating the nodes. Nodes are
%   labelled with their indices.
%   Fig1, right: Plot of the automatically detected fibre-fibre crossing
%   points, plotted as black circles on top of the original image.
%   Fig2, left: Back projection of the detected network, overlaid on top of
%   the original image. The back projection allows you to visualise the
%   detected thickness of each fibre.
%   Fig2, right: Overlay of detected fibres on top of original image,
%   coloured according to a user defined property. Edit the code in this
%   section to define which property (field of the fibreProps structure)
%   should be shown.
%
%   INPUTS:
%       -noteNodes: Structure containing information about the graph's
%       nodes, with fibres annotated by the annotateNetwork function.
%       -noteLinks: Structure containing information about the graph's
%       links, with fibres annotated by the annotateNetwork function.
%       -fibreProps: Measured properties of each fibre, as output by the
%       measureFibres function.
%       -origZ: Original heightmap AFM image.
%
%   Author: Oliver J. Meacock, (c) 2021

%The first figure will be an overlay of the annotated graph on top of the
%original image
subplot(1,2,1)
ax = gca;
visualiseAnnotatedNetwork(noteLinks,noteNodes,origZ,ax)

%The second figure will plot each node proposed to be a crossing point
subplot(1,2,2)
cla
ax=gca;
hold on
imagesc(ax,origZ)
colormap('jet')
axis equal
axis tight
ax.XTick = [];
ax.YTick = [];
ax.Box = 'on';
ax.YDir = 'reverse';

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
            plot(noteNodes(n).ptComx,noteNodes(n).ptComy,'o','MarkerFaceColor','k','MarkerEdgeColor','w')
        end
    end
end

title('Fibre crossing points')

%The third figure will show each detected fibre 'back projected' to show
%its projected width
widReconFac = 10;

%Reconstruct image showing width of fibers at each detected location
fibWidImg = zeros(size(origZ));
for l = 1:size(noteLinks,2)
    fibWidImg(noteLinks(l).point) = noteLinks(l).widths;
end
fibWidImg = fibWidImg';

backProjImg = zeros(size(origZ));
widSet = unique(fibWidImg(:));
widSet = widSet(2:end); %Get rid of leading 0

for w = widSet'
    currLines = fibWidImg == w;
    se = strel('disk',round(w/widReconFac));
    currProj = imdilate(currLines,se);
    backProjImg = or(backProjImg,currProj);
end

dispImg = cat(3,backProjImg*0.75,origZ/max(origZ(:)),backProjImg*0.75);

figure(2)
subplot(1,2,1)
imshow(dispImg)

title('Back-projected network')

%The fourth figure will overlay a backprojection of each fibre
subplot(1,2,2)
ax = gca;

visualiseAnnotatedFibres(fibreProps,noteNodes,flattenedImg,'localOrientation',ax)


