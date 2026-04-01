function [Xfit, Yfit, Zfit, Xgof, Ygof, Zgof, Intensity] = V2_ANSrun_fitfoci_nofig_denoised_DNAMERFISH(SampleNum)

tic
disp(['Processing SampleNum = ' num2str(SampleNum)])

%% 
% Load calibration parameters

load(['DriftParams' num2str(SampleNum) '.mat']);  % Xdrift(i), Ydrift(i), Zdrift(i) %fiducial beads
load('DeltaZ.mat');                               % scalar z-offset between 560 and 640 %Tetraspeck/calib beads
load('tform.mat');                                % 2D geometric transform (640 → 560 frame) %Tetraspeck/calib beads

% thresholds
LocalMaxThresh_560 = 20;
LocalMaxThresh_640 = 10;

% Preallocate all outputs for 24 probes
Xfit      = cell(1,24); %creates a 1-by-24 cell array
Yfit      = cell(1,24);
Zfit      = cell(1,24);
Xgof      = cell(1,24);
Ygof      = cell(1,24);
Zgof      = cell(1,24);
Intensity = cell(1,24);

%% 
% FIT 560 CHANNEL (Probes 1–12)
disp('Fitting 560 channel...')

for hyb = 1:12
    probeID = hyb;   % 560 -> probes 1–12

    FileName = ['sequential/560_' num2str(hyb) '_' num2str(SampleNum) '.mat'];
    load(FileName)  % loads ImageStack560
    disp(['Loaded 560 file: ' FileName])

    % Denoise
    ImageStack560 = imgaussfilt3(ImageStack560, 1);

    % Fit foci
    disp(['Fitting 560 loci, hyb ' num2str(hyb)])
    [Xfit{probeID}, Yfit{probeID}, Zfit{probeID}, ...
        Xgof{probeID}, Ygof{probeID}, Zgof{probeID}, Intensity{probeID}] = ...
        fitMultipleFoci(ImageStack560, LocalMaxThresh_560);

    % Drift correction (indexed by hyb)
    Xfit{probeID} = Xfit{probeID} - Xdrift(hyb);
    Yfit{probeID} = Yfit{probeID} - Ydrift(hyb);
    Zfit{probeID} = Zfit{probeID} - Zdrift(hyb);
end

disp('Finished fitting 560 channel')


%% 
% FIT 640 CHANNEL (Probes 13–24)
disp('Fitting 640 channel...')

for hyb = 1:12
    probeID = hyb + 12;   % 640 → probes 13–24

    FileName = ['sequential/640_' num2str(hyb) '_' num2str(SampleNum) '.mat'];
    load(FileName)  % loads ImageStack640
    disp(['Loaded 640 file: ' FileName])

    % Warp 640 images to 560 coordinate system in XY
    for m = 1:size(ImageStack640,3)
        ImageStack640(:,:,m) = imtransform(ImageStack640(:,:,m), tform, ...
            'XData',[1 1608], 'YData',[1 1608]);
    end

    % Denoise
    ImageStack640 = imgaussfilt3(ImageStack640, 1);

    % Fit foci
    disp(['Fitting 640 loci, hyb ' num2str(hyb)])
    [Xfit{probeID}, Yfit{probeID}, Zfit{probeID}, ...
        Xgof{probeID}, Ygof{probeID}, Zgof{probeID}, Intensity{probeID}] = ...
        fitMultipleFoci(ImageStack640, LocalMaxThresh_640);

    % Drift correction (indexed by hyb)
    Xfit{probeID} = Xfit{probeID} - Xdrift(hyb);
    Yfit{probeID} = Yfit{probeID} - Ydrift(hyb);
    Zfit{probeID} = Zfit{probeID} - Zdrift(hyb);

    % Z-channel mismatch: warp 640 → 560 in Z
    Zfit{probeID} = Zfit{probeID} - DeltaZ;
end

disp('Finished fitting 640 channel')


%% 
% PIXEL → MICRON CONVERSION
px_xy = 0.11;   % µm
px_z  = 0.20;   % µm

disp('Converting coordinates to microns...')

for probeID = 1:24
    if isempty(Zfit{probeID})
        continue
    end

    Zfit{probeID} = Zfit{probeID} * px_z;
    Xfit{probeID} = Xfit{probeID} * px_xy;
    Yfit{probeID} = Yfit{probeID} * px_xy;

    % Flip Y axis into microscope coordinate system
    Yfit{probeID} = (1608 * px_xy) - Yfit{probeID};
end


%% 
% SAVE RESULTS
outname = ['result_' num2str(SampleNum) '.mat'];
disp(['Saving: ' outname])

save(outname, 'Xfit','Yfit','Zfit','Xgof','Ygof','Zgof','Intensity')

toc
end
