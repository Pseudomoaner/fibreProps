clear all
close all

input = 'C:\Users\olijm\Desktop\Laia analysis\MachineLearnTest\networkMeasures.mat';

load(input)

datFields = fieldnames(networkMeasures);
measFields = fieldnames(networkMeasures.(datFields{1}));

PCAmatrix = zeros(numel(datFields),numel(measFields));

for i = 1:size(PCAmatrix,1)
    for j = 1:size(PCAmatrix,2)
        PCAmatrix(i,j) = networkMeasures.(datFields{i}).(measFields{j});
    end
end

scPCAmat = PCAmatrix./repmat(std(PCAmatrix),size(PCAmatrix,1),1); %Nondimensionalise by dividing by the std of each variable
[coeff,score,latent] = pca(scPCAmat);

figure(1)
plot(score(:,1),score(:,2),'r.')
hold on
for i = 1:size(score,1)
    text(score(i,1) + 0.1,score(i,2),datFields{i},'Interpreter','none')
end

title('PCA')
xlabel('Component 1')
ylabel('Component 2')

ax1 = gca;
ax1.LineWidth = 1.5;
axis equal

Y = tsne(scPCAmat,'Algorithm','barneshut');

figure(2)
plot(Y(:,1),Y(:,2),'r.')
hold on
for i = 1:size(score,1)
    text(Y(i,1) + 10,Y(i,2),datFields{i},'Interpreter','none')
end

title('t-SNE')

ax2 = gca;
ax2.LineWidth = 1.5;
axis equal