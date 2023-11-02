clc, clear, close all

%% Target colors

XYZ = zeros(3,44);
XYZ(1,1:11) = 0:0.1:1;
XYZ(2,12:22) = 0:0.1:1;
XYZ(3,23:33) = 0:0.1:1;
XYZ(:,34:44) = repmat(0:0.1:1,3,1);

% XYZ(3,34:44) = XYZ(3,34:44) + 0.1;
% XYZ(XYZ>1) = 1;

figure, plot(XYZ')
axis tight

%% Convert from XYZ to display space

load T_xyz1931.mat T_xyz1931 S_xyz1931 % Requires PsychToolbox

% Load screen characterization data

% screenSPD = readmatrix("../displayCharacterization/measurements/data_init_0.csv");
% S_screenSPD = [384,4,100]; % This is slightly wacky - it should have a value at 380, but it seems to get lost somewhere along the conversion... (TODO Look into this)

% load("SpectralMeasurement231030-120356.mat",'SPD','S_SPD'); % old characterization
% load("SpectralMeasurement231030-161536.mat",'SPD','S_SPD'); % personal machine
load("SpectralMeasurement231031-151235.mat",'SPD','S_SPD');

screenSPD = SPD;
S_screenSPD = S_SPD;

% Convert to XYZ

for i = 1:4
    screenSPDint(:,:,i) = SplineSpd(S_screenSPD,screenSPD(:,:,i),S_xyz1931); % should this be SplineSpd/SplineSrf/SplineRaw (what's the difference?)
end

figure, hold on
for i = 1:4
    plot(SToWls(S_screenSPD),screenSPD(:,:,i),'k:');
    plot(SToWls(S_xyz1931),screenSPDint(:,:,i),'r:');
end
xlabel('Wavelength (nm)')
ylabel('Transmisson')

for i = 1:4
    screenXYZ(:,:,i) = T_xyz1931*screenSPDint(:,:,i);
end

screenXYZNormFactor = screenXYZ(2,end,end);
save('screenXYZNormFactor','screenXYZNormFactor')
screenXYZ = (screenXYZ/screenXYZNormFactor);

for i = 1:4
    screenxyY(:,:,i) = XYZToxyY(screenXYZ(:,:,i)); % Note, Y is not normalised here
end

figure, hold on
scatter3(screenXYZ(1,:),screenXYZ(2,:),screenXYZ(3,:))
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
plot3([0,1],[0,1],[0,1])
% figure,
% scatter3(screenXYZ(1,37:54),screenXYZ(2,37:54),screenXYZ(3,37:54)) %just blue

figure, hold on
DrawChromaticity()
for i = 1:4
    scatter3(screenxyY(1,:,i),screenxyY(2,:,i),screenxyY(3,:,i),'k','filled')
end

%% Create conversion matrix

% xr = screenxyY(1,18,1);
% yr = screenxyY(2,18,1);
% xg = screenxyY(1,18,2);
% yg = screenxyY(2,18,2);
% xb = screenxyY(1,18,3);
% yb = screenxyY(2,18,3);
% xw = screenxyY(1,18,4);
% yw = screenxyY(2,18,4);
% 
% XYZToRGBMat = XYZToRGBMatrix(xr, yr, xg, yg, xb, yb, xw, yw);

% XYZblack = squeeze(screenXYZ(:,1,[1,2,3]));

RGBToXYZMat = squeeze(screenXYZ(:,end,[1,2,3]));% - XYZblack;
XYZToRGBMat = inv(RGBToXYZMat); 
% I prefer this way (it's cleaner)

disp(XYZToRGBMat);
disp(XYZToRGBMat*squeeze(screenXYZ(:,end,[1,2,3])))

% figure,
% imagesc(M)

% figure,
% imagesc(XYZToRGBMat)

%% XYZ to RGB

RGBlin = XYZToRGBMat*XYZ;

%% Compute linearization LUT

screenxyY_reshape = reshape(screenxyY,[3,18,4]);
for i = 1:4
    screenxyY_reshape_norm(:,:,i) = screenxyY_reshape(:,:,i)./screenxyY_reshape(:,end,i);
end

cols = {'r','g','b'};

% figure, hold on
for i = 1:3
    figure(100+i), hold on
    plot(0:15:255,screenxyY_reshape_norm(3,:,i),...
        'o-','Color',cols{i})
    xlabel('gunVal')
    ylabel('Y')
end

x2 = 0:255;
for i = 1:3
    figure(100+i)
    % p(i,:) = polyfit(0:15:255,screenxyY_reshape_norm(3,:,i),2);
    LUT(:,i) = spline(0:15:255,screenxyY_reshape_norm(3,:,i),x2);
    plot(x2,LUT(:,i),...
        '--','Color',cols{i},'LineWidth',2)
end

% x2 = linspace(0,255,1000);
% y2 = polyval(p,x2);
% plot(x2,y2,'k')

%% Test

XYZrequested = [0.5,0.5,0.5];

for i = 1:3
    [~,minloc(i)] = min(abs(XYZrequested(i) - LUT(:,i)));
    minloc(i) = minloc(i)-1; % to correct for indexing starting at 1
    % disp(minloc(i))
    figure(100+i)
    xline(minloc(i))
    yline(0.5)
end

%% RGBlin to gunVals, using LUT

gunVals = zeros(size(RGBlin));
for i = 1:size(RGBlin,2)
    for j = 1:3
        [~,gunVals(j,i)] = min(abs(RGBlin(j,i) - LUT(:,j)));
    end
end
gunVals = gunVals-1; % to correct for indexing starting at 1

for i = 1:3
    figure,
    histogram(gunVals(i,:))
end

figure,
histogram(gunVals)

%%
figure,
imshow(uint8(permute(gunVals,[3,2,1])))


save('gunVals','gunVals')

