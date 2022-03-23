function [] = visualiseAnnotatedNetworkGUI(noteLinks,noteNodes,origZ,ax)
%RECONSTRUCTFIBRENETWORK creates a visual reconstruction of an annotated
%AFM fibre graph. Mostly used for debugging purposes.
%
%   INPUTS:
%       -noteLinks: annotated network link structure, as output by
%       annotateNetwork.m
%       -noteNodes: annotated node structure, as output by
%       annotateNetwork.m
%       -origZ: The flattened version of the original AFM heightmap.
%       -ax: Handle to axes into which you want to plot.
%
%   Author: Oliver J. Meacock, 2022

hold(ax,'on')
imagesc(ax,origZ)
colormap(ax,'turbo')

caxis(ax,[prctile(origZ(:),1),prctile(origZ(:),99)])

%Begin with links. Colour each fibre a separate random colour.
fibCols = rand(size(noteLinks,2),3); %Much larger than needed, but avoids needing to run through data first

for i = 1:size(noteLinks,2)
    plot(ax,noteLinks(i).comx,noteLinks(i).comy,'w','LineWidth',2)
    plot(ax,noteLinks(i).comx,noteLinks(i).comy,'Color',[0.35,0.35,0.35],'LineWidth',1.5)
    midInd = round(size(noteLinks(i).comx,2)/2);
%     text(noteLinks(i).comx(midInd),noteLinks(i).comy(midInd),['L',num2str(i)],'Color','k')
end

%And then nodes
for i = 1:size(noteNodes,2)
    plot(ax,noteNodes(i).ptComx,noteNodes(i).ptComy,'w.','MarkerSize',16)
    plot(ax,noteNodes(i).ptComx,noteNodes(i).ptComy,'c.','MarkerSize',13)
%     text(noteNodes(i).ptComx,noteNodes(i).ptComy,['N',num2str(i)],'Color',[0,0.65,0.65])
end

% title('Network reconstruction')