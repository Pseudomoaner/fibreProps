clear all
close all

root = 'C:\Users\olijm\Desktop\Laia analysis\MachineLearnTest';
inputFiles = {'I6_Ori'};%{'WT','I1_An','I2_An','I3_An','I4_An','I5_An','I1_Ori','I2_Ori','I3_Ori','I4_Ori','I6_Ori'};
inputExtension = '.txt';
outputExtension = '.tif';

visualise = true;

widFac = 5; %Factor by which the measured 'width' of fibres (based on the scale of a Gaussian filter) should be divided to get the real width (in pixels). Needs manual callibration.
flattenScale = 15; %Scale of the upper DOG filter

for i = 1:size(inputFiles,2)
    [AFMmat,dx] = txtToMat(root,[inputFiles{i},inputExtension]);
    dx = 1; %Need to set to real units once they have been properly recorded
    
    %First flatten the image
    origZvals = AFMmat(:,:,3) - min(min(AFMmat(:,:,3))); %Sets the smallest value to zero
    flattenedImg = origZvals - imgaussfilt(origZvals,flattenScale); %Remove low spatial frequency information. Could use e.g. polynomial fitting instead of this DoG approach.
    flattenedImg = flattenedImg - min(flattenedImg(:)); %Sets the smallest value to zero (again)
    
    %Then extract the fine-scale network
    [netNodes,netLinks,fibreGroups] = networkReconstruction(flattenedImg,origZvals,dx);
    
    %Find long fibres directly from network structure
    [noteNodes,noteLinks] = annotateNetwork(netNodes,netLinks,fibreGroups,dx);
    
    %Measure properties of long fibres
    fibreProps.(inputFiles{i}) = measureFibres(noteNodes,noteLinks,origZvals,dx);
    
    if visualise
        %Reconstruct network visualisation to check it looks OK
        backProjImg = visualiseAnnotatedNetwork(noteNodes,noteLinks,fibreProps.(inputFiles{i}),flattenedImg');
    end
    
    %Collate measures and save
    networkMeasures.(inputFiles{i}) = measureNetwork(noteNodes,noteLinks,fibreProps.(inputFiles{i}),origZvals,dx,widFac);
    
    %Save the flattened image, for ease of visual comparison
    imwrite(uint8(round((flattenedImg-max(flattenedImg(:)))*255/(min(flattenedImg(:))-max(flattenedImg(:))))),fullfile(root,[inputFiles{i},outputExtension]))
end

save(fullfile(root,'networkMeasures.mat'),'networkMeasures','fibreProps')