%% BI 260 Image Processing and Analysis Mammography Project
% Jacob Ellison, Ricardo Flores, Walter Zhao

clear
close all
clc

%% Read images

fprintf("IMAGE ACQUISITION\n")
cd '/Users/rickflores/Documents/BI 260/mdb200data'
pgmFiles = dir('/Users/rickflores/Documents/BI 260/mdb200data/*.pgm');
imgSize = 1024;
num = 200;
imgList = zeros(imgSize,imgSize,num);
for i = 1:num
    fprintf("reading image " + i + "\n")
    imgList(:,:,i) = imread(pgmFiles(i).name);
end
imgList = double(imgList) ./ max(imgList,[],'all');

left = zeros(imgSize,imgSize,num./2);
right = zeros(imgSize,imgSize,num./2);
j = 1;
k = 1;
for i = 1:size(imgList,3)
    if mod(i,2) == 1
        left(:,:,j) = imgList(:,:,i);
        j = j + 1;
    else
        right(:,:,k) = imgList(:,:,i);
        k = k + 1;
    end
end
cd '/Users/rickflores/Documents/BI 260'
fprintf("\n")

%% Extract ground truth data

fprintf("GROUND TRUTH EXTRACTION\n")
info = fopen('/Users/rickflores/Documents/BI 260/mdb200data/Info200.txt');
infoRead = textscan(info, '%s', 'Delimiter', '\n');
fclose(info);
infoRead = vertcat(infoRead{:});
mdb = find(contains(infoRead, 'mdb'));
spec = zeros(length(mdb),1);
for i = mdb(1):mdb(end)
    specimen = char(infoRead(i));
    spec(i + 1 - mdb(1)) = str2num(specimen(4:6));
end
numMass = [];
for i = spec(1):spec(end)
    numMass = [numMass; sum(spec(:) == i)];
end
diag = zeros(length(mdb),1);
density = zeros(length(mdb),1);
tumor = zeros(length(mdb),1);
for i = mdb(1):mdb(end)
    if contains(infoRead(i), 'NORM')
    else
        diag(i + 1 - mdb(1)) = 1;
    end
    if contains(infoRead(i), 'F')
        density(i + 1 - mdb(1)) = 1;
    elseif contains(infoRead(i), 'G')
        density(i + 1 - mdb(1)) = 2;
    else
        density(i + 1 - mdb(1)) = 3;
    end
    if contains(infoRead(i), ' M ')
        tumor(i + 1 - mdb(1)) = 2;
    elseif contains(infoRead(i), ' B ')
        tumor(i + 1 - mdb(1)) = 1;
    else
        tumor(i + 1 - mdb(1)) = 0;
    end
end
for i = mdb(2):mdb(end)
    if spec(i) == spec(i - 1)
        density(i) = [];
        diag(i) = [];
        tumor(i) = [];
    end
end
fprintf("\n")

%% Process images

fprintf("IMAGE PROCESSING\n")
num = 200; % change this to read more images
resize = 256;
processedList = zeros(resize,resize,num);
imgThreshList = zeros(resize,resize,num);
for i = 1:num
    
    fprintf("processing image " + i + "\n")
    
    img = imgList(:,:,i);
    
    % flip right mammograms
    if mod(i,2) == 0
        img = flip(img,2);
    end
    
    % extract image sizes
    [xSize,ySize] = size(img);
    totSize = xSize.*ySize;
    
    % normalize intensities and sharpen using unsharp masking
    img = img ./ max(img,[],'all');
    img = imsharpen(img,'Amount',9);
    
    % binarize image
    imgBW = imbinarize(img);
    
    % reform image using binary mask
    for j = 1:totSize
        if imgBW(j) == 0
            img(j) = 0;
        end
    end
    
    % remove artifacts to isolate breast
    CC = bwconncomp(img);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    for k = 1:CC.NumObjects
        if k ~= idx
            img(CC.PixelIdxList{k}) = 0;
        end
    end
    
    % resize image and multiple threshold
    img = imresize(img,[resize resize]);
    thresh = multithresh(img,2);
    imgThresh = imquantize(img,thresh);
    imgThresh = imgThresh./max(imgThresh,[],'all');

    imgBW = imbinarize(imgThresh,0.34); % binarize image
    
    % reform image using binary mask
    [xSize,ySize] = size(img);
    totSize = xSize.*ySize;
    for j = 1:totSize
        if imgBW(j) == 0
            imgThresh(j) = 0;
        end
    end
    
    % find first point on line
    [m,n] = size(img);
    start = round(150.*(n./1024));
    locY1 = 0;
    for a = n:-1:1
        if imgThresh(start,a) == 1
            locY1 = a;
            break
        end
    end
    
    if locY1 ~= 0
        locX1 = 0;
        for a = start:m
            if imgThresh(a,locY1) ~= 1
                locX1 = a;
                break
            end
        end
        if locX1 == 0
            locX1 = m;
        end
        pos1 = [locX1 locY1];
        
        % find second point on line
        locX2 = 0;
        for a = 1:m
            for b = 1:n
                if imgThresh(a,b) == 1
                    locX2 = a;
                    break
                end
                if locX2 ~= 0
                    break
                end
            end
        end
        locY2 = 0;
        for a = 1:locY1
            if imgThresh(locX2,a) == 1
                locY2 = a;
                break
            end
        end
        if locX1 == 0
            locX1 = 1;
        end
        pos2 = [locX2 locY2];
        
        % define line
        m = (locX1 - locX2)./(locY1 - locY2);
        b = locX1 - m.*locY1;
        
        % segment out pectoral muscle
        [x,y] = size(img);
        for u = 1:x
            for v = 1:y
                if u <= m*v + b
                    img(u,v) = 0;
                    imgThresh(u,v) = 0;
                end
            end
        end
        
    end
    
    if mod(i,2) == 0
        img = flip(img,2);
        imgThresh = flip(imgThresh,2);
    end
    processedList(:,:,i) = img;
    imgThreshList(:,:,i) = imgThresh;
    
end
fprintf("\n")

%% Extract density features

fprintf("FEATURE EXTRACTION\n")
ratioList = zeros(num,1);
centerList = zeros(num,2);
spreadList = zeros(num,2);
gaborMeanList = [];
gaborStdList = [];
for i = 1:num
    
    fprintf("extracting features from image " + i + "\n")
    
    img = processedList(:,:,i);
    imgThresh = imgThreshList(:,:,i);
    
    % extract muscle to fat ratio
    muscleCount = length(find(imgThresh == 1));
    fatCount = (length(find(imgThresh > 0)) - length(find(imgThresh == 1)));
    muscleFatRatio = muscleCount ./ fatCount;
    ratioList(i) = muscleFatRatio;
    
    % extract center of muscle tissue
    [row,col] = find(imgThresh == 1);
    rowMean = mean(row);
    colMean = mean(col);
    center = [rowMean colMean];
    centerList(i,:) = center;
    
    % extract spread of muscle tissue
    rowStd = std(row);
    colStd = std(col);
    spread = [rowStd colStd];
    spreadList(i,:) = spread;
    
    % Gabor filter features
    wavelength = 20; % choose angular orientations
    endLength = 180 - wavelength; % choose end orientation
    numLength = endLength./wavelength + 1;
    orientation = [0:wavelength:endLength];
    g = gabor(wavelength,orientation); % generate Gabor filter bank
    mag = imgaborfilt(img,g); % apply Gabor filter bank
    magMean = [];
    magStd = [];
    for z = 1:numLength % extract features
        magMean = [magMean mean(mag(:,:,z),'all')];
        magStd = [magStd std(mag(:,:,z),[],'all')];
    end
    gaborMeanList = [gaborMeanList; magMean];
    gaborStdList = [gaborStdList; magStd];
    
end
fprintf("\n")

%% Train density classifier

fprintf("DENSITY CLASSIFICATION\n")

iter = 10; % choose iteration number
resub = zeros(iter,1);
gen = zeros(iter,1);
correct = zeros(iter,1);
confusion = zeros(3);

for i = 1:iter
    
    fprintf("iteration " + i + "\n")
    dFeatures = [ratioList centerList spreadList gaborMeanList gaborStdList]; % se
    train = round(0.7.*num); % choose to train using 70% of data
    randTrain = randperm(num,train); % randomly select training set
    randTest = [];
    for j = 1:num % assign rest of data to testing set.
        if ~ismember(j,randTrain)
            randTest = [randTest j];
        end
    end
    dFeaturesTrain = dFeatures(randTrain,:);
    densityTrain = density(randTrain);
    dLDA = fitcdiscr(dFeaturesTrain,densityTrain); % train LDA classifier
    resubError = resubLoss(dLDA); % calculate resubstitution error
    resub(i) = resubError; 
    CVdSVM = crossval(dLDA); % cross-validate classifier
    genError = kfoldLoss(CVdSVM); % calculate generalized error
    gen(i) = genError;
    [label,~,~] = predict(dLDA,dFeatures(randTest,:)); % generate predicted labels
    correct = label == density(randTest); % determine correct labels
    countCorrect = length(find(correct == 1)); % calculate number of correct labels
    percentCorrect = countCorrect ./ length(correct); % calculate percent correct
    correct(i) = percentCorrect;
    C = confusionmat(density(randTest),label); % generate confusion matrix
    confusion = confusion + C;
    
end

% calculate statistics over all iterations
resub = mean(resub);
gen = mean(gen);
correct = mean(correct);
confusion = round(confusion ./ iter);

fprintf("\nTotal iterations: " + iter + "\n")
fprintf("\nAverage resubstitution error: " + resub + "\n")
fprintf("\nAverage generalized error: " + gen + "\n")
fprintf("\nAverage percent correct: " + correct + "\n")
fprintf("\n")

figure(1)
CC = confusionchart(confusion); % visualize confusion matrix
CC.Title = ['Breast Density Classification Using LDA Over ' num2str(iter) ' Iterations'];
CC.RowSummary = 'row-normalized';
CC.ColumnSummary = 'column-normalized';

%% Display images

fprintf("DISPLAY IMAGES\n")
figure(2), sgtitle('Breast Mammograms')

subplot(211)
montage(imgList,'Size',[7 46]), title('Original Mammograms')
colormap gray, axis tight equal

subplot(212)
montage(processedList,'Size',[7 46]), title('Processed Mammograms')
colormap gray, axis tight equal

fprintf("\n")
