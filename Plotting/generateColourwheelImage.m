function outImg = generateColourwheelImage(cmap)
noCs = size(cmap,1);
[xGrid,yGrid] = meshgrid(-500:500,0:500);
[thet,r] = cart2pol(xGrid,yGrid);
thetList = -pi/((noCs-1)*2):pi/(noCs-1):(pi+pi/((noCs-1)*2));
cwheel = zeros([size(thet),3]);
cmap(cmap == 0) = 0.001;

%Apply map to wedges
for i = 1:noCs
    cwheelInds = and(and(thet>thetList(i),thet<=thetList(i+1)),r<500);
    cwheel(:,:,1) = cwheel(:,:,1) + cwheelInds*cmap(i,1);
    cwheel(:,:,2) = cwheel(:,:,2) + cwheelInds*cmap(i,2);
    cwheel(:,:,3) = cwheel(:,:,3) + cwheelInds*cmap(i,3);
end

%Make black bits white
cwheelInds = cwheel(:,:,1) == 0;
cwheel(:,:,1) = cwheel(:,:,1) + cwheelInds;
cwheel(:,:,2) = cwheel(:,:,2) + cwheelInds;
cwheel(:,:,3) = cwheel(:,:,3) + cwheelInds;

outImg = imrotate(cwheel,180);