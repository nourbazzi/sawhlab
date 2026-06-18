system('curl -o readNPY.m https://raw.githubusercontent.com/kwikteam/npy-matlab/master/npy-matlab/readNPY.m');
system('curl -o readNPYheader.m https://raw.githubusercontent.com/kwikteam/npy-matlab/master/npy-matlab/readNPYheader.m');
addpath(pwd);

HYB           = 0;
NUM_SERIES    = 15;
CHANNEL_NAMES = {'SD561', 'SD488', 'SD405'};

for series = 1:15
%for series = 20
    fprintf('Processing series %d / %d...\n', series, NUM_SERIES);
    for ch = 1:3
        chName  = CHANNEL_NAMES{ch};
        npyFile = sprintf('sequential/%s_%d_%d.npy', chName, HYB, series);

        chStack  = readNPY(npyFile);
        chStack  = permute(chStack, [2 3 1]);

        ImageMax = max(chStack, [], 3);

        fig = figure('Visible', 'off');
        imagesc(ImageMax); axis equal; colormap gray;
        title(sprintf('%s_Hyb%d_FOV%d', chName, HYB, series));
        saveas(fig, sprintf('sequential/%s_Hyb%d_FOV%d', chName, HYB, series));
        close(fig);

        saveData.(sprintf('ImageStack%s', chName)) = chStack;
        save(sprintf('sequential/%s_%d_%d.mat', chName, HYB, series), ...
            '-struct', 'saveData');
        clear saveData chStack ImageMax;
    end
end
fprintf('All done.\n');