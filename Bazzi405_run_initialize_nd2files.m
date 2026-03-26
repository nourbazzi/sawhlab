% This script opens and converts files from nd2 format to multidemensional matrices for downstream
% processing
% ------------------------------------------------------------------------------------------------------------------------------------------
% Requirements: bfmatlab installed and in the path
% ------------------------------------------------------------------------------------------------------------------------------------------
% Input: path to .nd2 file with multiple series(FOV/samples) and illumination channels
% 
% series: FOV/position
% 
% 
% ------------------------------------------------------------------------------------------------------------------------------------------
% Output: Multidemensional matrices for each channel, series, and hyb#
% ------------------------------------------------------------------------------------------------------------------------------------------
% Ahilya N. Sawh, PhD
% 
%% ------------------------------------------------------------------------------------------------------------------------------------------
clear all 
close all
tic

mkdir sequential_denoised

NumSeries = 6;
NumChannels = 1;
sz = 1608;
zst = 101;
Hyb = 0
%% primary probe file



data = bfopen('ChannelWF 640,WF 561,WF 488,WF 405_Seq0000 - WF405Denoised.nd2');

for series = 1:NumSeries

    ImageStack = zeros(sz,sz,zst*NumChannels);
    for i = 1:zst*NumChannels %planes are organized by z position first then channel
        ImageStack(:,:,i) = data{series, 1}{i,1};
    end

    % %channel 1 = 640
    % ImageStack640 = zeros(sz,sz,zst);
    % ImageStack640 = ImageStack(:,:,1:4:zst*NumChannels);
    % ImageMax = max(ImageStack640,[],3);
    % figure
    % imagesc(ImageMax)
    % axis equal
    % colormap gray
    % title(['640_Hyb' num2str(Hyb) '_FOV' num2str(series)])
    % saveas(gcf,['sequential/640_Hyb' num2str(Hyb) '_FOV' num2str(series)])
    % save( ['sequential/640_' num2str(Hyb) '_' num2str(series) '.mat'],'ImageStack640');
    % clear ImageMax ImageStack640
    % 
    % %channel 2 = 560
    % ImageStack560 = zeros(sz,sz,zst);
    % ImageStack560 = ImageStack(:,:,2:4:zst*NumChannels);
    % ImageMax = max(ImageStack560,[],3);
    % figure
    % imagesc(ImageMax)
    % axis equal
    % colormap gray
    % title(['560_Hyb' num2str(Hyb) '_FOV' num2str(series)])
    % saveas(gcf,['sequential/560_Hyb' num2str(Hyb) '_FOV' num2str(series)])
    % save( ['sequential/560_' num2str(Hyb) '_' num2str(series) '.mat'],'ImageStack560');
    % clear ImageMax ImageStack560
    % 
    % %channel 3 = 488
    % ImageStack488 = zeros(sz,sz,zst);
    % ImageStack488 = ImageStack(:,:,3:4:zst*NumChannels);
    % ImageMax = max(ImageStack488,[],3);
    % figure
    % imagesc(ImageMax)
    % axis equal
    % colormap gray
    % title(['488_Hyb' num2str(Hyb) '_FOV' num2str(series)])
    % saveas(gcf,['sequential/488_Hyb' num2str(Hyb) '_FOV' num2str(series)])
    % save( ['sequential/488_' num2str(Hyb) '_' num2str(series) '.mat'],'ImageStack488');
    % clear ImageMax ImageStack488

    %channel 4 = 405
    ImageStack405 = zeros(sz,sz,zst);
    ImageStack405 = ImageStack(:,:,1:1:zst*NumChannels);
    ImageMax = max(ImageStack405,[],3);
    figure
    imagesc(ImageMax)
    axis equal
    colormap gray
    title(['405_Hyb' num2str(Hyb) '_FOV' num2str(series)])
    saveas(gcf,['sequential_denoised/405_Hyb' num2str(Hyb) '_FOV' num2str(series)])
    save( ['sequential_denoised/405_' num2str(Hyb) '_' num2str(series) '.mat'],'ImageStack405');
    clear ImageStack ImageMax ImageStack405
    close all

end
clear data




toc