function anaPred = prepareFibresForManualMeasurement(flatImg,fibreMeasures)

maxVal = max(flatImg(:));

anaPred = double.empty(0,2);

for F = 1:size(fibreMeasures,2)
    if ~isnan(fibreMeasures(F).midwidth) && fibreMeasures(F).size > 40
        currImg = flatImg;
        currImg(logical(fibreMeasures(F).backbone)) = maxVal;
        
        cla
        imshow(currImg,[])
        hold on
        plot(fibreMeasures(F).midpoint(2),fibreMeasures(F).midpoint(1),'r.','MarkerSize',12)
        title(['Fibre index: ', num2str(F)])
        
        export_fig(['C:\Users\olijm\Desktop\Laia analysis\WidthAnalysis_Ex2\Fibre_',num2str(F),'.jpg'],'-m0.9989')
        
        anaPred = [anaPred;F,fibreMeasures(F).midwidth];
    end
end