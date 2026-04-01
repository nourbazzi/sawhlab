% This script calculates the sample drift parameters for each sample/FOV in
% sequential FISH by fitting a Gaussian curve to the TetraSpeck/FluoSphere bead signal in each
% round 
% ------------------------------------------------------------------------------------------------------------------------------------------
% Input: multidimensional matrices (3D image stacks as mat files) of beads in each round
% of sequential FISH
%
% SampleNum is the FOV, NumHybs is the number of rounds of secondary hybs
% ------------------------------------------------------------------------------------------------------------------------------------------
% Output: DriftParams file for each sample/FOV containing drift in x,y,z
% from Hyb0 (primary probe round)
% ------------------------------------------------------------------------------------------------------------------------------------------
% Ahilya N. Sawh, PhD
% 2020
% 
% using Zhuang lab fitFoci and selectROI functions
%% ------------------------------------------------------------------------------------------------------------------------------------------
clear all
close all

SampleNum = 6; %try a different sample here
NumHybs = 15; %try a lower number of hybs to start then increase up to 14 if time allows

FileName = ['sequential/488_0_' num2str(SampleNum) '.mat']; 

load(FileName)

ImageMax = max(ImageStack488(:,:,1:20),[],3);
figure(1)
imagesc(ImageMax)
colormap gray
caxis([100 1500])
axis equal
title('bead in hyb0')

selectROI

figure(2)
imagesc(ImageMax)
colormap gray
axis equal
title('bead in hyb0 overlay with all hyb positions')
hold on

for i = 0:NumHybs
    FileName = ['sequential/488_' num2str(i) '_' num2str(SampleNum) '.mat']; 
    load(FileName)
    [Xfit(i+1), Yfit(i+1), Zfit(i+1)] = fitFoci(ImageStack488(:,:,1:20), roiList, i+1,1);
    figure(2)
    j = jet(NumHybs+1);
    plot(Xfit(i+1), Yfit(i+1), 'x', 'MarkerEdgeColor', j(i+1,:));
end

hold off

for i = 1:NumHybs
    Xdrift(i) = Xfit(i+1) - Xfit(1);
    Ydrift(i) = Yfit(i+1) - Yfit(1);
    Zdrift(i) = Zfit(i+1) - Zfit(1);
end

save(['DriftParams' num2str(SampleNum) '.mat'], 'Xdrift', 'Ydrift', 'Zdrift');
    
