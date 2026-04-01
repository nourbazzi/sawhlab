function [FinalFociMatrix, CentroidMatrix] = trackFoci3D(SampleNum, codebookname)

book = readcell(codebookname);

numRounds = 24;
WidthThreshMax = 4; % Remove "foci" that are too wide
WidthThreshMin = 1; % Remove "foci" that are too narrow
AdjrsquareThreshold = 0.7;
focusThresh = 2; % 3D distance threshold for linking foci

% Convert SampleNum to string for file operations
sampleStr = num2str(SampleNum);

% Load data
dataFile = ['results/result' sampleStr '.mat'];

load(dataFile, 'Xfit', 'Yfit', 'Zfit', 'Intensity', ...
               'Xgof', 'Ygof', 'Zgof');

% Initialize linked foci storage
LinkedFoci = {}; 

% Process each round
for i = 1:numRounds
    if isempty(Xfit{i}) || isempty(Yfit{i}) || isempty(Zfit{i})
        continue; % Skip empty rounds
    end

    % Filtering step
    Ind = find( [Xgof{i}.adjrsquare] > AdjrsquareThreshold & ...
                [Ygof{i}.adjrsquare] > AdjrsquareThreshold & ...
                [Zgof{i}.adjrsquare] > AdjrsquareThreshold);

    % Apply filtering
    Xfit{i} = Xfit{i}(Ind);
    Yfit{i} = Yfit{i}(Ind);
    Zfit{i} = Zfit{i}(Ind);
    Intensity{i} = Intensity{i}(Ind);
    Xgof{i} = Xgof{i}(Ind);
    Ygof{i} = Ygof{i}(Ind);
    Zgof{i} = Zgof{i}(Ind);

    % Load directly into P1 (no conversion?)
    P1 = [Xfit{i}', Yfit{i}', Zfit{i}'];

    % If first round, initialize tracking
    if i == 1  
        for j = 1:size(P1,1)
            LinkedFoci{j} = nan(numRounds, 3); % Store (x, y, z) across rounds
            LinkedFoci{j}(i,:) = P1(j,1:3); % Save first appearance
        end
    else
        % Compare current round (i) foci to previous rounds
        for j = 1:size(P1,1)
            foundMatch = false;
            for k = 1:length(LinkedFoci)
                dist = sum((LinkedFoci{k}(1,1:3) - P1(j,1:3)).^2); 
                if dist < focusThresh^2
                    LinkedFoci{k}(i,:) = P1(j,1:3); % Save focus if match
                    foundMatch = true;
                    break;
                end
            end
            if ~foundMatch % No match = create new entry
                LinkedFoci{end+1} = nan(numRounds, 3);
                LinkedFoci{end}(i,:) = P1(j,1:3); 
            end
        end
    end
end

% Create matrix with all foci
FinalFociMatrix = nan(length(LinkedFoci), numRounds * 3); 
for j = 1:length(LinkedFoci)
    FinalFociMatrix(j, :) = reshape(LinkedFoci{j}', 1, []); 
end

% find centroids
CentroidMatrix = cell(length(LinkedFoci), 4); % Store (X_centroid, Y_centroid, Z_centroid, BinaryCode)
for j = 1:length(LinkedFoci)
    CentroidMatrix{j,1} = mean(LinkedFoci{j}(:,1), 'omitnan'); % X centroid
    CentroidMatrix{j,2} = mean(LinkedFoci{j}(:,2), 'omitnan'); % Y centroid
    CentroidMatrix{j,3} = mean(LinkedFoci{j}(:,3), 'omitnan'); % Z centroid
    
    % Track binary code as a character string
    binaryVector = ~isnan(LinkedFoci{j}(:,1)); % 1 if detected, 0 if missing
    binaryString = char('0' + binaryVector'); % Convert logical array to '0'/'1' characters
    CentroidMatrix{j,4} = binaryString'; % Store as character string
end

% Save results 
disp(['Saving results for Sample ' sampleStr]);

save(['TrackedFoci3D_' sampleStr '.mat'], 'FinalFociMatrix');
save(['CentroidMatrix_' sampleStr '.mat'], 'CentroidMatrix');

for j = 1:length(CentroidMatrix)
    binaryString = CentroidMatrix{j,4}; % Get the binary string for this entry
    
    % Find the row in the codebook that matches the binary string
    %rowIndex = find(ismember(book(:,1), binaryString), 1); % Assuming binary string is in the first column of the codebook
    rowIndex = find(ismember(book(:,1), {binaryString}), 1);
    
    if ~isempty(rowIndex)
        % Get the corresponding column number from the codebook
        columnNumber = rowIndex; % Assuming the codebook has the column numbers in the first column
        % Append the column number to the CentroidMatrix entry
        CentroidMatrix{j,5} = columnNumber; 
    else
        CentroidMatrix{j,5} = NaN;
    end
end
disp(['Saving decoded centroid matrix for Sample ' sampleStr]);
save(['centroid_matrix_decoded_' sampleStr '.mat'], 'CentroidMatrix');










toc
end