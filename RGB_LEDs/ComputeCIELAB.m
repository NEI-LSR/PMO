clear, clc, close all

% testingRoomWall_SPD     = load("SpectralMeasurement231016-162128.mat",'SPD');
% testingRoomWall_SPD     = testingRoomWall_SPD.SPD;
% testingRoomWall_S_SPD    = load("SpectralMeasurement231016-162128.mat",'S_SPD');
% testingRoomWall_S_SPD    = testingRoomWall_S_SPD.S_SPD;
% 
% % today, for more direct comparison (pointing at different bit of wall etc)
% testingRoomWall_SPD2     = load("SpectralMeasurement231102-124627",'SPD');
% testingRoomWall_SPD2     = testingRoomWall_SPD2.SPD;
%
% % an hour or so later...
% testingRoomWall_SPD3     = load("SpectralMeasurement231102-160304",'SPD');
% testingRoomWall_SPD3     = testingRoomWall_SPD3.SPD;
% 
% testingRoomWall_SPD4     = load("SpectralMeasurement231102-160527",'SPD');
% testingRoomWall_SPD4     = testingRoomWall_SPD4.SPD;
% 
% testingRoomWall_SPD5     = load("SpectralMeasurement231102-160749",'SPD');
% testingRoomWall_SPD5     = testingRoomWall_SPD5.SPD;
% 
% testingRoomWall_SPD6     = load("SpectralMeasurement231102-160839",'SPD');
% testingRoomWall_SPD6     = testingRoomWall_SPD6.SPD;

testingRoomWall_SPD     = load("SpectralMeasurement231102-161024.mat",'SPD');
testingRoomWall_SPD     = testingRoomWall_SPD.SPD;

testingRoomWall_S_SPD    = load("SpectralMeasurement231102-161024.mat",'S_SPD');
testingRoomWall_S_SPD    = testingRoomWall_S_SPD.S_SPD;

figure, hold on
% 
% plot(SToWls(testingRoomWall_S_SPD),testingRoomWall_SPD)
% plot(SToWls(testingRoomWall_S_SPD),testingRoomWall_SPD2)
% plot(SToWls(testingRoomWall_S_SPD),testingRoomWall_SPD3)
% plot(SToWls(testingRoomWall_S_SPD),testingRoomWall_SPD4)
% plot(SToWls(testingRoomWall_S_SPD),testingRoomWall_SPD5)
% plot(SToWls(testingRoomWall_S_SPD),testingRoomWall_SPD6,'--')
plot(SToWls(testingRoomWall_S_SPD),testingRoomWall_SPD)

% legend
%% Load measurements

% TODO REPLACE THIS WITH REAL MEASUREMENTS
d = dir('SpectralMeasurement231101*.mat'); % For 53,32 and 73,52

%% Convert to XYZ

load T_xyz1931.mat T_xyz1931 S_xyz1931 % Requires PsychToolbox

SPDint = SplineSpd(S_SPD,SPD,S_xyz1931); % should this be SplineSpd/SplineSrf/SplineRaw (what's the difference?)
testingRoomWall_SPDint = SplineSpd(S_SPD,testingRoomWall_SPD,S_xyz1931);

XYZ = T_xyz1931*SPDint;
testingRoomWall_XYZ = T_xyz1931*testingRoomWall_SPDint;

sRGBlin = XYZToSRGBPrimary(XYZ./testingRoomWall_XYZ(2)); % TODO Consider what the implied white point is
sRGB = uint8(SRGBGammaCorrect(sRGBlin,0)');

%% Convert to CIELAB

Lab = XYZToLab(XYZ,testingRoomWall_XYZ);

figure,
scatter3(Lab(2,:),Lab(3,:),Lab(1,:),...
    [],double(sRGB)/255,'filled','MarkerEdgeColor','k')
daspect([1,1,1])
xlabel('a*')
ylabel('b*')
view(2)
title('CIELab')

%%
