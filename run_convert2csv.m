% Specify the folder where the files live.
myFolder = 'Ncoords_1';
mkdir Ncoords_1_txt
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolder);
    uiwait(warndlg(errorMessage));
    myFolder = uigetdir(); % Ask for a new one.
    if myFolder == 0
         % User clicked Cancel
         return;
    end
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.mat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
fprintf(1, 'Now reading folder %s\n', myFolder);
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    % Now do whatever you want with this file name,
    % such as reading it in as an image array with imread()
%     print(fullFileName)
%     imageArray = imread(fullFileName);
%     imshow(imageArray);  % Display image.
%     drawnow; % Force display to update immediately.
    % Load a mat file
    fprintf(1, 'Basename %s\n', baseFileName);
    load(fullFileName)
%     % Write a txt file
%     newname = ''
%     writetable(nstats)
    % writetable(nstats,'nstats_tab.txt','Delimiter','tab')
    % Clear workspace
    ix = strfind(baseFileName,'.')
    name_pre = baseFileName(1:ix(1)-1);
    name_txt = strcat('Ncoords_1_txt/nstats', name_pre, '.txt')
    writetable(nstats,name_txt)%,'Delimiter','tab')

    name_txt = strcat('Ncoords_1_txt/pcellnucID', name_pre, '.txt')
    writematrix(pcellnucID,name_txt)%,'Delimiter','tab')

    name_txt = strcat('Ncoords_1_txt/somcellnucID', name_pre, '.txt')
    writematrix(somcellnucID,name_txt)%,'Delimiter','tab')

    name_txt = strcat('Ncoords_1_txt/pcellncenter', name_pre, '.txt')
    writematrix(pcellncenter,name_txt)%,'Delimiter','tab')

    clear nstats pcellnucID somcellnucID pcellncenter;
end
