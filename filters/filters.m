clear, 
%clc, 
close all

saveFigs = true;

%% Load data

data = readtable('C:\Users\cege-user\Documents\PMO\filters\CORRECT ORDER GELS MEASUREMENT SPECTRUMS 3200 5600 - DBR V1 - 21 12 2020.xlsx');
% This by default loads the first sheet only (which is the 3200K measurements)

lightSource = readtable('C:\Users\cege-user\Documents\PMO\filters\3200k AND 5600k REFERENCE SOURCES SPECTRUM - DBR V1.xlsx',...
    'VariableNamingRule','preserve');

% testingRoomWall_SPD     = load("SpectralMeasurement231011-114850.mat",'SPD');
% testingRoomWall_SPD     = testingRoomWall_SPD.SPD;
% testingRoomWall_S_PD    = load("SpectralMeasurement231011-114850.mat",'S_SPD');
% testingRoomWall_S_PD    = testingRoomWall_S_PD.S_SPD;

testingRoomWall_SPD     = load("SpectralMeasurement231016-162128.mat",'SPD');
testingRoomWall_SPD     = testingRoomWall_SPD.SPD;
testingRoomWall_S_PD    = load("SpectralMeasurement231016-162128.mat",'S_SPD');
testingRoomWall_S_PD    = testingRoomWall_S_PD.S_SPD;

%% Normalise and plot SPDs

SPD_raw = table2array(data(2:end,13:end));
S_SPD = WlsToS(table2array(data(1,13:end))');

% figure, hold on
% plot(SPD_raw(195,:))
% plot(table2array(lightSource(:,2)))

NormaliseByClear = false;
if NormaliseByClear
    SPD = (SPD_raw./SPD_raw(195,:))';

    figure, hold on
    plot(SPD)
    plot(SPD(:,195),'k','LineWidth',2)
    axis tight
end

normaliseByLightSource = true;
if normaliseByLightSource
    SPD = (SPD_raw./table2array(lightSource(:,2))')';

    figure, hold on
    plot(SPD)
    plot(SPD(:,195),'k','LineWidth',2)
    axis tight
end

%%

SPD = [SPD,...
    SPD.*SPD(:,48),... % light grey
    SPD.*SPD(:,49),... % pale grey
    SPD.*SPD(:,50),... % medium grey
    SPD.*SPD(:,530)]; % neutral grey

% figure, 
% % plot(SToWls(S_SPD),SPD);
% plot(SToWls(S_SPD),SPD(:,1:100)); % subset, easier to see
% xlabel('Wavelength (nm)')
% ylabel('Transmisson')

%% Compute effect in testing room
% Take into account the light source, 
% and use the back wall as a test reflector

SPD = testingRoomWall_SPD.*SPD;

figure, hold on
plot(SPD)
plot(SPD(:,195),'k','LineWidth',2)
axis tight

% load spd_houser
% 
% SPD = SplineSpd(S_houser,spd_houser(:,22),S_SPD).*SPD;
% 
% figure, hold on
% plot(SPD)
% plot(SPD(:,195),'k','LineWidth',2)
% axis tight

%% Compute CIE xy_1931

load T_xyz1931.mat T_xyz1931 S_xyz1931 % Requires PsychToolbox

SPDint = SplineSpd(S_SPD,SPD,S_xyz1931); % should this be SplineSpd/SplineSrf/SplineRaw (what's the difference?)

% figure, 
% % plot(SToWls(S_SPD),SPD);
% plot(SToWls(S_xyz1931),SPDint(:,1:100)); % subset, easier to see
% xlabel('Wavelength (nm)')
% ylabel('Transmisson')

XYZ = T_xyz1931*SPDint;
xyY = XYZToxyY(XYZ);

% Compute sRGB (for plotting)
sRGBlin = XYZToSRGBPrimary(XYZ./max(XYZ(2,:))); % TODO Consider what the implied white point is
sRGB = uint8(SRGBGammaCorrect(sRGBlin,0)');

figure,
DrawChromaticity
scatter3(xyY(1,:),xyY(2,:),xyY(3,:),...
    [],double(sRGB)/255,'filled','MarkerEdgeColor','k')
daspect([1,1,200])
xlim([0,0.8])
ylim([0,0.9])
zlabel('Y_{1931}')
view(2)
title('CIE xy')

if saveFigs
    print(gcf,'-vector','-dsvg',['1931xy','.svg'])
end

% Plot SPDs again, but this time using the sRGB colors for plotting
figure, hold on
for i = 1:size(SPD,2)
    % plot(SToWls(S_SPD),SPD(:,i),'Color',double(sRGB(i,:))/255)
    plot3(SToWls(S_SPD),SPD(:,i),ones(size(SPD(:,i)))*i,...
        'Color',double(sRGB(i,:))/255)
end
xlabel('Wavelength (nm)')
ylabel('Transmisson')
axis tight
title('SPD')

if saveFigs
    print('SPD.png','-dpng')
end

%% u'v'

upvp = xyTouv(xyY(1:2,:));

figure,
DrawChromaticity('upvp')
scatter3(upvp(1,:),upvp(2,:),xyY(3,:),...
    [],double(sRGB)/255,'filled','MarkerEdgeColor','k')
daspect([1,1,200])
xlim([0,0.6])
ylim([0,0.6])
xlabel('u''')
ylabel('v''')
zlabel('Y_{1931}')
view(2)
title('CIE upvp')

if saveFigs
    print(gcf,'-vector','-dsvg',['upvp','.svg'])
end

%% CIELUV

Luv = XYZToLuv(XYZ,XYZ(:,195)); % taking the clear filter as the white point - right decision? (TODO)
Lab = XYZToLab(XYZ,XYZ(:,195));

figure,
scatter3(Luv(2,:),Luv(3,:),Luv(1,:),...
    [],double(sRGB)/255,'filled','MarkerEdgeColor','k')
% zlim([60,70])
daspect([1,1,1])
xlabel('u*')
ylabel('v*')
view(2)
title('CIELuv')

% zlim([45,55])

if saveFigs
    print(gcf,'-vector','-dsvg',['CIELuv','.svg'])
end

figure,
scatter3(Lab(2,:),Lab(3,:),Lab(1,:),...
    [],double(sRGB)/255,'filled','MarkerEdgeColor','k')
daspect([1,1,1])
xlabel('a*')
ylabel('b*')
view(2)
title('CIELab')

% zlim([60,70])

if saveFigs
    print(gcf,'-vector','-dsvg',['CIELab','.svg'])
end

%% DKL

load T_cones_sp.mat T_cones_sp S_cones_sp
load T_xyzJuddVos.mat T_xyzJuddVos S_xyzJuddVos
T_Y = 683*T_xyzJuddVos(2,:);
T_Y = SplineCmf(S_xyzJuddVos,T_Y,S_cones_sp);

LMS     = T_cones_sp * SplineSpd(S_SPD,SPD,S_cones_sp);
bgLMS   = mean(LMS,2); % TODO Replace this
LMSinc = LMS - bgLMS;

[M_ConeIncToDKL,LMLumWeights] = ComputeDKL_M(bgLMS,T_cones_sp,T_Y);

% DKL = M_ConeIncToDKL*LMS;
DKL = M_ConeIncToDKL*LMSinc;

figure,
scatter3(DKL(2,:),DKL(3,:),DKL(1,:),...
    [],double(sRGB)/255,'filled','MarkerEdgeColor','k')
xlabel('L-M')
ylabel('S-(L+M)')
zlabel('L+M')
title('DKL')
view(2)

% note: relative axis scaling is a choice

if saveFigs
    print(gcf,'-vector','-dsvg',['DKL','.svg'])
end

% zlim([0, 2])

%% Select colors

% Convert angles that isolate DKL mechanisms into LMS
% Convter from LMS to XYZ, then CIELUV

% Or... just pick specific filters from DKL, and see where they are in
% CIELUV?

% Or both?

% Let's start with just picking from CIELUV, for simplicity

radius = 70;
Lstar = 60;

requestedLocations = [Lstar,0,0;...
    Lstar,radius,0;...
    Lstar,0,radius;...
    Lstar,-radius,0;...
    Lstar,0,-radius];


nPoints = 12;
requestedLocations = zeros(nPoints,3);

requestedLocations(:,1) = ones(nPoints,1)*Lstar;
[requestedLocations(:,2),requestedLocations(:,3)] = pol2cart(0:2*pi/nPoints:2*pi-(2*pi/nPoints),radius);
requestedLocations(end+1,:) = [Lstar,0,0];    

figure,  hold on
scatter3(requestedLocations(:,2),requestedLocations(:,3),requestedLocations(:,1))

closestInd = zeros(size(requestedLocations,1),1);

for i = 1:length(closestInd)
    [~,closestInd(i)] = min(sqrt(sum((requestedLocations(i,:)'-Luv).^2)));
end

scatter3(Luv(2,closestInd),Luv(3,closestInd),Luv(1,closestInd),...
    [],double(sRGB(closestInd,:))/255,'filled','MarkerEdgeColor','k')
daspect([1,1,1])

xlabel('u*')
ylabel('v*')
view(2)
title('CIELuv')

% Which filter(s) are those?

figure, hold on
for j = 1:length(closestInd)
    i = closestInd(j);
    % plot(SToWls(S_SPD),SPD(:,i),'Color',double(sRGB(i,:))/255)
    plot3(SToWls(S_SPD),SPD(:,i),ones(size(SPD(:,i)))*i,...
        'Color',double(sRGB(i,:))/255)
end
xlabel('Wavelength (nm)')
ylabel('Transmisson')
axis tight
title('SPD')



NDfilterInd = floor(closestInd/size(SPD_raw,1)) - 1
filterInd = mod(closestInd,size(SPD_raw,1));

data.Var1(filterInd+1)
data.Var4(filterInd+1)
data.Var5(filterInd+1)

%%

