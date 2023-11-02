clear, clc, close all

folder = 'C:\Users\cege-user\Documents\PMO\screenStimuli\images\';
% d = dir([folder,'IMG*.tif']);

d = dir([folder,'IMG_5656.tif']);


%%

Lstar = [0,0,-20];%[73,73,53];
radius = [32,52,32];

nPoints = 8;

for stimset = 1
    requestedLocations = zeros(nPoints,3);
    requestedLocations(:,1) = ones(nPoints,1)*Lstar(stimset);
    [requestedLocations(:,2),requestedLocations(:,3)] = pol2cart(0:2*pi/nPoints:2*pi-(2*pi/nPoints),radius(stimset));
    requestedLocations(end+1,:) = [Lstar(stimset),0,0];

    for im = 1:length(d)
        for filter = 1:size(requestedLocations,1)
            sRGBinput = imread([folder,d(im).name]);
            outputIm = imageProcessingPipeline(sRGBinput,requestedLocations(filter,:));
            for repeat = 1%:10
                rng(repeat)
                randomNumber = randperm(9);
                imwrite(outputIm,[folder,'filteredImages\',d(im).name(1:end-4),'_',...
                    num2str(repeat),'_',num2str(randomNumber(filter)),'_',num2str(filter),'_',...
                    num2str(Lstar(stimset)),'_',num2str(radius(stimset)),'.tif'])
            end
        end
    end
end

%