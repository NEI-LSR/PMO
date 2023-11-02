function outputIm = imageProcessingPipeline(sRGBinput,LabOffset)

% clc, clear, close all

%%

% sRGBinput = imread('C:\Users\cege-user\Documents\PMO\screenStimuli\_MG_9169.tif');
% sRGBinput = imread('C:\Users\cege-user\Documents\PMO\screenStimuli\IMG_5644.tif');
% sRGBinput(1,1,:) = 255; % trying forcing the top value to be 255 because it seems like the sRGB conversion scales to max, and so assumes that there will be pixel values at max

% figure,imshow(sRGBinput)
sRGBreshape = reshape(sRGBinput,[size(sRGBinput,1)*size(sRGBinput,2),size(sRGBinput,3)]);

sRGBlin = SRGBGammaUncorrect(sRGBreshape);
XYZ = SRGBPrimaryToXYZ(sRGBlin');

% % test hack
% hack = zeros(size(XYZ));
% hack(3,:) = 0.05;
% XYZ = XYZ + hack;

xy = XYZToxyY(XYZ);

sRGBgamut_sRGB = [1,0,0;0,1,0;0,0,1];
sRGBgamut_sRGBlin = SRGBGammaUncorrect(sRGBgamut_sRGB);
sRGBgamut_XYZ = SRGBPrimaryToXYZ(sRGBgamut_sRGBlin');
sRGBgamut_xy = XYZToxyY(sRGBgamut_XYZ);

% figure, hold on
% DrawChromaticity
% % scatter(xy(1,1:100:end),xy(2,1:100:end),'k','filled',...
% %     'MarkerEdgeAlpha',0.1,'MarkerFaceAlpha',0.1)
% scatter3(xy(1,1:100:end),xy(2,1:100:end),xy(3,1:100:end),...
%     [],double(sRGBreshape(1:100:end,:))/255,'filled')
% scatter(sRGBgamut_xy(1,:),sRGBgamut_xy(2,:),'r*')
% plot([sRGBgamut_xy(1,:),sRGBgamut_xy(1,1)],[sRGBgamut_xy(2,:),sRGBgamut_xy(2,1)],'r:')

%% Convert from XYZ to display space

load T_xyz1931.mat T_xyz1931 S_xyz1931 % Requires PsychToolbox

% Load screen characterization data

% screenSPD = readmatrix("../displayCharacterization/measurements/data_init_0.csv");
% S_screenSPD = [384,4,100]; % This is slightly wacky - it should have a value at 380, but it seems to get lost somewhere along the conversion... (TODO Look into this)

load("SpectralMeasurement231031-151235",'SPD','S_SPD');
screenSPD = SPD;
S_screenSPD = S_SPD;

% Convert to XYZ

for i = 1:4
    screenSPDint(:,:,i) = SplineSpd(S_screenSPD,screenSPD(:,:,i),S_xyz1931); % should this be SplineSpd/SplineSrf/SplineRaw (what's the difference?)
end

% figure, hold on
% for i = 1:4
%     plot(SToWls(S_screenSPD),screenSPD(:,:,i),'k');
%     plot(SToWls(S_xyz1931),screenSPDint(:,:,i),'r--');
% end
% xlabel('Wavelength (nm)')
% ylabel('Transmisson')
% axis tight

for i = 1:4
    screenXYZ(:,:,i) = T_xyz1931*screenSPDint(:,:,i);
end

screenXYZNormFactor = screenXYZ(2,end,end);
% save('screenXYZNormFactor','screenXYZNormFactor')
screenXYZ = (screenXYZ/screenXYZNormFactor);

for i = 1:4
    screenxyY(:,:,i) = XYZToxyY(screenXYZ(:,:,i));
end

% figure,
% scatter3(screenXYZ(1,:),screenXYZ(2,:),screenXYZ(3,:))
% figure,
% scatter3(screenXYZ(1,37:54),screenXYZ(2,37:54),screenXYZ(3,37:54)) %just blue

% figure, hold on
% for i = 1:4
%     scatter3(screenxyY(1,:,i),screenxyY(2,:,i),screenxyY(3,:,i),'k','filled')
% end

%% Convert into Lab, modify, and convert back to XYZ

testingRoomWall_XYZ = [0.978261507670876, 1, 0.727287896622727]';

Lab = XYZToLab(XYZ,testingRoomWall_XYZ); 

% LabNorm = Lab;
% LabNorm(1,:) = LabNorm(1,:) - mean(LabNorm(1,:)); % normalise L* (but preserve a* and b*)
% LabShift = LabNorm + repmat(LabOffset',1,size(XYZ,2));
LabShift = Lab + repmat(LabOffset',1,size(XYZ,2));


% figure, hold on
% scatter(LabNorm(2,:),LabNorm(3,:),'k.')
% scatter(LabShift(2,:),LabShift(3,:),'r.')
% axis equal
% xline(0)
% yline(0)

XYZ = LabToXYZ(LabShift,testingRoomWall_XYZ);

% xy_shifted = XYZToxyY(XYZ);

% figure, hold on
% DrawChromaticity
% scatter3(xy_shifted(1,1:100:end),xy_shifted(2,1:100:end),xy_shifted(3,1:100:end))
% scatter(sRGBgamut_xy(1,:),sRGBgamut_xy(2,:),'r*')
% plot([sRGBgamut_xy(1,:),sRGBgamut_xy(1,1)],[sRGBgamut_xy(2,:),sRGBgamut_xy(2,1)],'r:')


%% Create conversion matrix

% -- % PTB method

% xr = screenxyY(1,18,1);
% yr = screenxyY(2,18,1);
% xg = screenxyY(1,18,2);
% yg = screenxyY(2,18,2);
% xb = screenxyY(1,18,3);
% yb = screenxyY(2,18,3);
% xw = screenxyY(1,18,4);
% yw = screenxyY(2,18,4);
% 
% XYZToRGBMat1 = XYZToRGBMatrix(xr, yr, xg, yg, xb, yb, xw, yw);

% -- % Own method

% XYZblack = squeeze(screenXYZ(:,1,[1,2,3]));

RGBToXYZMat = squeeze(screenXYZ(:,end,[1,2,3]));% - XYZblack;
XYZToRGBMat = inv(RGBToXYZMat); 
% I prefer this way (it's cleaner)

% -- % Bruce Lindbloom method (haven't got this to work yet)

% screenXYZgamut = squeeze(screenXYZ(:,end,[1,2,3]));
% S = screenXYZgamut^-1 * squeeze(screenXYZ(:,end,4));

% XYZToRGBMat3 = squeeze(screenXYZ(:,end,[1,2,3])).*S';

% Katie's method:
% XYZToRGBMat3 = [S(1)*screenXYZgamut(1,1) S(2)*screenXYZgamut(1,2) S(3)*screenXYZgamut(1,3);...
% S(1)*screenXYZgamut(2,1) S(2)*screenXYZgamut(2,2) S(3)*screenXYZgamut(2,3);...
% S(1)*screenXYZgamut(3,1) S(2)*screenXYZgamut(3,2) S(3)*screenXYZgamut(3,3)];

% Direct from Bruce website:

% Xr = xr/yr;
% Yr = 1;
% Zr = (1 - xr - yr)/yr;
% Xg = xg/yg;
% Yg = 1;
% Zg = (1 - xg - yg)/yg;
% Xb = xb/yb;
% Yb = 1;
% Zb = (1 - xb - yb)/yb;
% 
%  RGBToXYZMat3= [...
%     S(1)*Xr, S(2)*Xg, S(3)*Xb;...
%     S(1)*Yr, S(2)*Yg, S(3)*Yb;...
%     S(1)*Zr, S(2)*Zg, S(3)*Zb];
% 
% XYZToRGBMat3 = inv(RGBToXYZMat3);

% disp(XYZToRGBMat);
% disp(XYZToRGBMat*squeeze(screenXYZ(:,end,[1,2,3])))

% try
%     figure,
%     imagesc(XYZToRGBMat1)
%     colorbar
% 
%     figure,
%     imagesc(XYZToRGBMat2)
%     colorbar
% 
%     figure,
%     imagesc(XYZToRGBMat3)
%     colorbar
% catch
% end

%% XYZ to RGB

RGBlin = XYZToRGBMat*XYZ;

% figure,
% histogram(RGBlin)
% title('RGBlin')
% 
% figure,
% histogram(sRGBinput)
% title('sRGB input')

%% Compute linearization LUT

screenxyY_reshape = reshape(screenxyY,[3,18,4]);
for i = 1:4
    screenxyY_reshape_norm(:,:,i) = screenxyY_reshape(:,:,i)./screenxyY_reshape(:,end,i);
end

cols = {'r','g','b'};

% figure, hold on
% for i = 1:3
%     figure(100+i), hold on
%     plot(0:15:255,screenxyY_reshape_norm(3,:,i),...
%         'o-','Color',cols{i})
%     xlabel('gunVal')
%     ylabel('Y')
% end

x2 = 0:255;
for i = 1:3
    % figure(100+i)
    LUT(:,i) = spline(0:15:255,screenxyY_reshape_norm(3,:,i),x2);
    % plot(x2,LUT(:,i),...
    %     '--','Color',cols{i},'LineWidth',2)
end

% x2 = linspace(0,255,1000);
% y2 = polyval(p,x2);
% plot(x2,y2,'k')

%% Test

% XYZrequested = [0.5,0.5,0.5];
% 
% for i = 1:3
%     [~,minloc(i)] = min(abs(XYZrequested(i) - LUT(:,i)));
%     minloc(i) = minloc(i)-1; % to correct for indexing starting at 1
%     % disp(minloc(i))
%     figure(100+i)
%     xline(minloc(i))
%     yline(0.5)
% end

%% RGBlin to gunVals, using LUT


% for i = 1:3
%     figure,
%     histogram(RGBlin(i,:))
% end

gunVals = zeros(size(RGBlin));
for i = 1:size(RGBlin,2)
    for j = 1:3
        [~,gunVals(j,i)] = min(abs(RGBlin(j,i) - LUT(:,j)));
    end
end
gunVals = gunVals-1; % to correct for indexing starting at 1

% for i = 1:3
%     figure,
%     histogram(gunVals(i,:))
% end

% figure,
% histogram(gunVals)

%%

outputIm = uint8(reshape(gunVals',[size(sRGBinput)]));

% figure,
% imshow(outputIm)



