function  [Xfit, Yfit, Zfit, Xgof, Ygof, Zgof, Intensity] = fitMultipleFoci(ImageStack,LocalMaxThresh)
% Assuming the foci in the image stack are well isolated, fit the 3D
% positions of all the foci.

% clear all
% close all
% FileName = ['Sequential/STORMCy3_01_00'];
% NumFrame = 206;
% [ImageStack, InfoFile] = ReadZStack(FileName,NumFrame,5);
% ImageMean = mean(ImageStack,3);
% figure(100)
% imagesc(ImageMean)
% colormap gray
% axis equal

%threshold to find local mixum in the image stack
NoRows = size(ImageStack, 1);
NoColumns = size(ImageStack, 2);
NoZ = size(ImageStack,3);
BW = zeros(NoRows, NoColumns, NoZ);
for i = 4: NoRows-3
    for j = 4: NoColumns-3
        for k = 4: NoZ-3
            if ImageStack(i,j,k)>LocalMaxThresh+... 
                (sum(sum(ImageStack(i-2:i+2, j-3, k)))+...
                sum(sum(ImageStack(i-2:i+2, j+3, k)))+...
                sum(sum(ImageStack(i-3, j-3:j+3, k)))+...
                sum(sum(ImageStack(i+3, j-3:j+3, k))))/24
%             sum(sum(sum(ImageStack(i-3:i+3, j-3:j+3, k-3:k+3),3),2),1)/7^3

                BW(i,j,k) = 1;
            end
        end
    end
end

% MaxBW = max(BW, [], 3);
% figure(200)
% imagesc(MaxBW)
% colormap gray
% axis equal

%identify connected ares as the individual foci
CC = bwconncomp(BW, 26);
S = regionprops(CC, 'Centroid');

%fit the 3D positions of each focus
j = 0;

for i = 1:length(S)
    X0 = round(S(i).Centroid(1));
    Y0 = round(S(i).Centroid(2));
    I = round(S(i).Centroid(3));
    CenterIntensity = ImageStack(Y0, X0, I);
    
    Data1 = mean(ImageStack(Y0,X0-3:X0+3,I),1);
    Data2 = mean(ImageStack(Y0-3:Y0+3,X0,I),2);
    
    GaussEqu = 'a*exp(-(x-b)^2/2/c^2)+d';
    StartPoint1 = [CenterIntensity X0 1 0];
    StartPoint2 = [CenterIntensity Y0 1 0];
    
    [f1, gof1] = fit((X0-3:X0+3)', Data1', GaussEqu, 'Start', StartPoint1);
    [f2, gof2] = fit((Y0-3:Y0+3)', Data2, GaussEqu, 'Start', StartPoint2);
    
    Zprofile = permute(ImageStack(Y0, X0, :), [3 1 2]);
    Irange = (I-3:I+3);
    Data3 = Zprofile(Irange);
    StartPoint3 = [CenterIntensity I 1 0];
    [f3, gof3] = fit(Irange', Data3, GaussEqu, 'Start', StartPoint3);
    
    if f1.b>15 && f1.b<NoColumns-14 && f2.b>15 && f2.b<NoRows-14 && f3.b>1 && f3.b<NoZ
        j = j+1;
        Xfit(j) = f1.b;
        Yfit(j) = f2.b;
        Zfit(j) = f3.b;
        Xgof(j) = gof1;
        Ygof(j) = gof2;
        Zgof(j) = gof3;
        Intensity(j) = CenterIntensity;
    end
end
if j==0
    Xfit = [];
    Yfit = [];
    Zfit = [];
    Xgof = [];
    Ygof = [];
    Zgof = [];
    Intensity = [];
end

