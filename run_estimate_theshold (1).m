% This script tests threshold values (LocalMaxThresh) for each probe channel, to be used in
% next program to fit foci in all Hyb rounds
% ------------------------------------------------------------------------------------------------------------------------------------------
% Input: multidimensional matrices (3D image stacks) of foci in a channel
% ------------------------------------------------------------------------------------------------------------------------------------------
% Output: none
% ------------------------------------------------------------------------------------------------------------------------------------------
% Ahilya N. Sawh, PhD
% 06.03.2019
% Version 1.0
% 
% using the Zhuang lab fitMultipleFoci function
%% ------------------------------------------------------------------------------------------------------------------------------------------
clear all
close all

%set xy size and z-stack size
xysz= 1608;
zsz= 101;

SampleNum = 2;
Hyb = 1;

LocalMaxThresh560 = 10;
LocalMaxThresh640 = 5;
%try otsu or adaptive thresh?

Channel = 560;
FileName = ['sequential/' num2str(Channel) '_' num2str(Hyb) '_' num2str(SampleNum) '.mat'];
load(FileName)

ImProbe = max(ImageStack560,[],3);
figure

imagesc(ImProbe)
colormap gray
axis equal
title('Z stack max image 560')%flattened image of nuclei full xyz stack

embryo = drawellipse('StripeColor','y');


BW = createMask(embryo);
Embryo = zeros(xysz,xysz,zsz);
for i = 1:zsz
test = ImageStack560(:,:,i);
test(BW == 0) = 0;
Embryo(:,:,i) = test;
end

figure
set(gcf,'Visible','on')
imagesc(max(Embryo,[],3))
colormap gray
axis equal
title('selected embryo')

roirect=drawrectangle('Color', 'r')
x1 = ceil(roirect.Position(1,1));
x2 = ceil(roirect.Position(1,1))+ceil(roirect.Position(1,3));
y1 = ceil(roirect.Position(1,2));
y3 = ceil(roirect.Position(1,2))+ceil(roirect.Position(1,4));

I = Embryo(y1:y3, x1:x2, :);
figure
ImProbe = max(I,[],3);
imagesc(ImProbe)
colormap gray
axis equal
%caxis([0 3500])
title('selected FOV 560')

noisyV = I;

filteredV = imgaussfilt3(I, 1); %gaussian filter

figure
subplot(1,2,1)
imagesc(max(noisyV,[],3))
colormap gray
caxis([0 300])
title('noisy 560')
subplot(1,2,2)
imagesc(max(filteredV,[],3))
colormap gray
caxis([0 300])
title('filtered 560')

ImageStack560 = filteredV;

ImageMax = max(ImageStack560,[],3);
figure(100)
imagesc(ImageMax)
colormap gray
axis equal
title('560 max') 

[Xfit, Yfit, Zfit, Xgof, Ygof, Zgof, Intensity] = fitMultipleFoci(ImageStack560,LocalMaxThresh560);
figure
hold on
plot(Xfit, Yfit, 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
for j = 1:length(Xfit)
    text(Xfit(j), Yfit(j), num2str(j));
end
hold off
title('560 max foci detected')

figure(100)
hold on
plot(Xfit, Yfit, 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
for j = 1:length(Xfit)
    text(Xfit(j), Yfit(j), num2str(j));
end
hold off
axis equal
set(gca, 'YDir', 'reverse')
title('560 max foci detected')
saveas(gcf,[num2str(Channel) '_' num2str(Hyb) '_' num2str(SampleNum) '_thresh' num2str(LocalMaxThresh560) '_imagemaxfoci'])


%Labels for QC
LabelVol = zeros(size(ImageStack560), 'logical');

radius = 1;

for i = 1:length(Xfit)
    x = round(Xfit(i));
    y = round(Yfit(i));
    z = round(Zfit(i));

    x = max(1, min(size(LabelVol,2), x));
    y = max(1, min(size(LabelVol,1), y));
    z = max(1, min(size(LabelVol,3), z));

    xs = max(1, x-radius):min(size(LabelVol,2), x+radius);
    ys = max(1, y-radius):min(size(LabelVol,1), y+radius);
    zs = max(1, z-radius):min(size(LabelVol,3), z+radius);

    LabelVol(ys, xs, zs) = true;
end


%Pixels Size of Prime95B at 100X (XYZ) (µm):0.11 x 0.11 x 0.20, 1608x1608
%pixels

Zfit = Zfit*0.2; % convert to µm
Xfit = Xfit*110/1000; %µm
Yfit = Yfit*110/1000; %µm
Yfit = 1608*110/1000-Yfit;
    
figure
scatter3(Xfit, Yfit, Zfit, 'ok', 'MarkerFaceColor', 'r');
hold off
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
ylim([0 176.88])
xlim([0 176.88])

%%
Channel = 640;
FileName = ['sequential/' num2str(Channel) '_' num2str(Hyb) '_' num2str(SampleNum) '.mat'];
load(FileName)

ImProbe2 = max(ImageStack640,[],3);
figure

imagesc(ImProbe2)
colormap gray
axis equal
title('Z stack max image 640')%flattened image of nuclei full xyz stack

embryo = drawellipse('StripeColor','y');


BW = createMask(embryo);
Embryo = zeros(xysz,xysz,zsz);
for i = 1:zsz
test = ImageStack640(:,:,i);
test(BW == 0) = 0;
Embryo(:,:,i) = test;
end

figure
set(gcf,'Visible','on')
imagesc(max(Embryo,[],3))
colormap gray
axis equal
title('selected embryo')

roirect=drawrectangle('Color', 'r')
x1 = ceil(roirect.Position(1,1));
x2 = ceil(roirect.Position(1,1))+ceil(roirect.Position(1,3));
y1 = ceil(roirect.Position(1,2));
y3 = ceil(roirect.Position(1,2))+ceil(roirect.Position(1,4));

I = Embryo(y1:y3, x1:x2, :);
figure
ImProbe2 = max(I,[],3);
imagesc(ImProbe2)
colormap gray
axis equal
%caxis([0 3500])
title('selected FOV 640')

noisyV = I;

filteredV = imgaussfilt3(I, 1); %gaussian filter
figure
subplot(1,2,1)
imagesc(max(noisyV,[],3))
colormap gray
caxis([0 300])
title('noisy 640')
subplot(1,2,2)
imagesc(max(filteredV,[],3))
colormap gray
caxis([0 300])
title('filtered 640')

ImageStack640 = filteredV;

ImageMax = max(ImageStack640,[],3);
figure
imagesc(ImageMax)
colormap gray
axis equal
title('640 max')
figure(200)
imagesc(ImageMax)
colormap gray
axis equal

[Xfit, Yfit, Zfit, Xgof, Ygof, Zgof, Intensity] = fitMultipleFoci(ImageStack640,LocalMaxThresh640);
figure
hold on
plot(Xfit, Yfit, 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
hold off

figure(200)
hold on
plot(Xfit, Yfit, 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
for j = 1:length(Xfit)
    text(Xfit(j), Yfit(j), num2str(j));
end
hold off
axis equal
set(gca, 'YDir', 'reverse')
title('640 max foci detected')
saveas(gcf,[num2str(Channel) '_' num2str(Hyb) '_' num2str(SampleNum) '_thresh' num2str(LocalMaxThresh640) '_imagemaxfoci'])

clear LabelVol
%Labels for QC
LabelVol = zeros(size(ImageStack640), 'logical');


radius = 1;

for i = 1:length(Xfit)
    x = round(Xfit(i));
    y = round(Yfit(i));
    z = round(Zfit(i));

    x = max(1, min(size(LabelVol,2), x));
    y = max(1, min(size(LabelVol,1), y));
    z = max(1, min(size(LabelVol,3), z));

    xs = max(1, x-radius):min(size(LabelVol,2), x+radius);
    ys = max(1, y-radius):min(size(LabelVol,1), y+radius);
    zs = max(1, z-radius):min(size(LabelVol,3), z+radius);

    LabelVol(ys, xs, zs) = true;
end

%Pixels Size of Prime95B at 100X (XYZ) (µm):0.11 x 0.11 x 0.20, 1608x1608
%pixels

Zfit = Zfit*0.2; % convert to µm
Xfit = Xfit*110/1000; %µm
Yfit = Yfit*110/1000; %µm
Yfit = 1608*110/1000-Yfit;
    
figure
scatter3(Xfit, Yfit, Zfit, 'ok', 'MarkerFaceColor', 'r');
hold off
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
ylim([0 176.88])
xlim([0 176.88])
