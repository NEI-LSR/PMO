clear, clc, close all

folder = 'C:\Users\cege-user\Documents\PMO\screenStimuli\images\';
d = dir([folder,'IMG*.tif']);

%%

for im = [1,3]%1:length(d)
    sRGBinput = imread([folder,d(im).name]);
    outputIm = imageProcessingPipeline(sRGBinput,[0,0,0]);
    figure, imshow(outputIm)
    randomPixels(:,:,:,im) = outputIm(randsample(size(outputIm,1),100),randsample(size(outputIm,2),100),:);
end

randomPixels_rsh = [randomPixels(:,:,:,1),randomPixels(:,:,:,3)];
figure,imshow(randomPixels_rsh)
randomPixels_rsh = reshape(randomPixels_rsh,100*200,3);

%%

% folder2 = 'C:\Users\cege-user\Documents\PMO\screenStimuli\images\filteredImages\';
folder2 = 'C:\Users\cege-user\Documents\PMO\screenStimuli\images\filteredImages\stimset3\';
d2 = dir([folder2,'IMG*.tif']);

for i = 1:length(d2)
    rng(i)
    randomNumber(i,:) = randperm(5);
end

for k = 1:5

    figure('Position', [50, 50, 2000/3, 1333/3]);

    axes('Position', [0, 0, 1, 1]);
    axis off;
    hold on

    for j = 1:2000
        drawElipse(randomPixels_rsh(randsample(size(randomPixels_rsh,1),1),:))
    end

    % axis image
    xlim([0.5,1.5])
    ylim([0.5,1.5])
    frame = getframe(gca);

    for i = 1:length(d2)

        % Save the figure as a TIFF image
        filename = [folder2,d2(i).name(1:end-4),'_random_ellipse_',num2str(randomNumber(i,k)),'.tif'];
        imwrite(permute(frame.cdata,[2,1,3]), filename, 'tif');

        close all
    end
end



%%
function drawElipse(col)
% Generate random parameters for the ellipse
a = rand()/7; % Semi-major axis length
b = rand()/7; % Semi-minor axis length
theta = rand() * 2 * pi; % Random rotation angle in radians

% Generate a random center point
center_x = rand()*2;
center_y = rand()*2;

% Generate points on the ellipse
t = linspace(0, 2 * pi, 100); % Angle values from 0 to 2*pi
x = a * cos(t);
y = b * sin(t);

% Rotate the points by the random angle
x_rotated = x * cos(theta) - y * sin(theta);
y_rotated = x * sin(theta) + y * cos(theta);

% Translate the points to the random center
x_translated = x_rotated + center_x;
y_translated = y_rotated + center_y;

% Plot the filled ellipse
fill(x_translated, y_translated, col,'LineStyle','none');

end



%