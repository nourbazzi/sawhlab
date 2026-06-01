% Batch-produce per-FOV 640 TIFFs for the napari PG folder.
% Iterates over every Sample{N}FOV{F}coords.mat under Ncoords_1/parsed/.
clear all
mkdir('tiffs_P_granules');
coordFiles = dir('Ncoords_1/parsed/Sample*FOV*coords.mat');
if isempty(coordFiles)
    error('No coords files found in Ncoords_1/parsed/ — did you run the cropping step?');
end

for k = 1:length(coordFiles)
    name = coordFiles(k).name;     % e.g. 'Sample1FOV2coords.mat'
    tok = regexp(name, 'Sample(\d+)FOV(\d+)coords\.mat', 'tokens', 'once');
    if isempty(tok), continue; end
    sampleNum = str2double(tok{1});
    FOV       = str2double(tok{2});

    % load THIS FOV's crop box
    load(fullfile('Ncoords_1/parsed', name), 'x1', 'x2', 'y1', 'y3');

    % load the matching 640 stack
    matPath = sprintf('sequential/640_0_%d.mat', sampleNum);
    if ~isfile(matPath)
        fprintf('Skip %s: %s not found\n', name, matPath);
        continue
    end
    load(matPath, 'ImageStack640');

    % crop to this FOV and cast to uint16 ONCE
    I = uint16(ImageStack640(ceil(y1):ceil(y3), ceil(x1):ceil(x2), :));

    outputFileName = sprintf('tiffs_P_granules/Raw_Sample%d_FOV%d_640.tif', sampleNum, FOV);
    if isfile(outputFileName)
        delete(outputFileName);   % avoid appending onto a stale file
    end

    for z = 1:size(I, 3)
        imwrite(I(:,:,z), outputFileName, ...
            'WriteMode', 'append', 'Compression', 'none');
    end

    fprintf('Saved %s  (%dx%dx%d)\n', outputFileName, size(I,1), size(I,2), size(I,3));
end
disp('Done.');