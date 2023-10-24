clear,
clc,
close all

saveFigs = false;

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
    % plot(SPD(:,195),'k','LineWidth',2)
    axis tight
end

%% Exclude SPDs that we don't have access to

SPD = SPD(:,        ~strcmp('E COLOR',table2array(data(2:end,4))) & ~strcmp('GAM',table2array(data(2:end,4))));
data_mod = data(    ~strcmp('E COLOR',table2array(data(:    ,4))) & ~strcmp('GAM',table2array(data(:    ,4))),:);
% data_mod = data;

figure, hold on
plot(SPD)
axis tight

%%

SPD = [SPD,...
    SPD.*SPD(:,find(strcmp('Light Grey',table2array(data_mod(:,1))))-1),... % light grey
    SPD.*SPD(:,find(strcmp('Pale Grey',table2array(data_mod(:,1))))-1),... % pale grey
    SPD.*SPD(:,find(strcmp('Medium Grey',table2array(data_mod(:,1))))-1),... % medium grey
    SPD.*SPD(:,find(strcmp('Neutral Grey',table2array(data_mod(:,1))))-1)]; % neutral grey

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
plot(testingRoomWall_SPD,'k','LineWidth',2)
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
testingRoomWall_SPDint = SplineSpd(S_SPD,testingRoomWall_SPD,S_xyz1931);

% figure,
% % plot(SToWls(S_SPD),SPD);
% plot(SToWls(S_xyz1931),SPDint(:,1:100)); % subset, easier to see
% xlabel('Wavelength (nm)')
% ylabel('Transmisson')

XYZ = T_xyz1931*SPDint;
xyY = XYZToxyY(XYZ); % Note, Y is not normalised here

testingRoomWall_XYZ = T_xyz1931*testingRoomWall_SPDint;
testingRoomWall_xyY = XYZToxyY(testingRoomWall_XYZ); % Note, Y is not normalised here

% Compute sRGB (for plotting)
sRGBlin = XYZToSRGBPrimary(XYZ./testingRoomWall_XYZ(2)); % TODO Consider what the implied white point is
sRGB = uint8(SRGBGammaCorrect(sRGBlin,0)');

figure(421)
DrawChromaticity
scatter3(xyY(1,:),xyY(2,:),xyY(3,:),...
    [],double(sRGB)/255,'filled','MarkerEdgeColor','k')
daspect([1,1,0.3])
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

Luv = XYZToLuv(XYZ,testingRoomWall_XYZ); 
Lab = XYZToLab(XYZ,testingRoomWall_XYZ);

figure(420), hold on
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

%% CIELUV and CIELAB subset (pilot colors)

whichFilters = {'CALCOLOR 15 GREEN','CALCOLOR 30 GREEN','CALCOLOR 60 GREEN','CALCOLOR 90 GREEN','CALCOLOR 15 RED','CALCOLOR 30 RED','CALCOLOR 60 RED','CALCOLOR 90 RED'};
for i = 1:length(whichFilters)
    whichFiltersInd(i) = find(strcmp(whichFilters{i},table2array(data_mod(:,1))))-1;
end

figure, hold on
scatter3(Luv(2,whichFiltersInd),Luv(3,whichFiltersInd),Luv(1,whichFiltersInd),...
    [],double(sRGB(whichFiltersInd,:))/255,'filled','MarkerEdgeColor','k')
% zlim([60,70])
daspect([1,1,1])
xlabel('u*')
ylabel('v*')
view(2)
title('CIELuv')

figure, hold on
scatter3(Lab(2,whichFiltersInd),Lab(3,whichFiltersInd),Lab(1,whichFiltersInd),...
    [],double(sRGB(whichFiltersInd,:))/255,'filled','MarkerEdgeColor','k')
% zlim([60,70])
daspect([1,1,1])
xlabel('a*')
ylabel('b*')
view(2)
title('CIELab')

%% DKL

load T_cones_sp.mat T_cones_sp S_cones_sp
load T_xyzJuddVos.mat T_xyzJuddVos S_xyzJuddVos
T_Y = 683*T_xyzJuddVos(2,:);
T_Y = SplineCmf(S_xyzJuddVos,T_Y,S_cones_sp);

LMS     = T_cones_sp * SplineSpd(S_SPD,SPD,S_cones_sp);
bgLMS   = T_cones_sp * SplineSpd(S_SPD,testingRoomWall_SPD,S_cones_sp);
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

figure(6523)
tiledlayout(3,3)

figure(6524), hold on
tiledlayout(3,3)

for Lstar = [75,65,55] %1:100 75 %
    for radius =  [32,42,52] %1:80 32 %
        % radius = 50; %[35,45,55]
        % Lstar = 65; %[50,65,80]

        nPoints = 8;
        requestedLocations = zeros(nPoints,3);

        requestedLocations(:,1) = ones(nPoints,1)*Lstar;
        [requestedLocations(:,2),requestedLocations(:,3)] = pol2cart(0:2*pi/nPoints:2*pi-(2*pi/nPoints),radius);
        requestedLocations(end+1,:) = [Lstar,0,0];

        figure(6524), hold on
        nexttile
        % figure(422), hold on
        scatter3(requestedLocations(:,2),requestedLocations(:,3),requestedLocations(:,1),'k')

        closestInd = zeros(size(requestedLocations,1),1);

        for i = 1:length(closestInd)
            [~,closestInd(i)] = min(sqrt(sum((requestedLocations(i,:)'-Lab).^2)));
        end

        figure(6524), hold on
        scatter3(Lab(2,closestInd),Lab(3,closestInd),Lab(1,closestInd),...
            [],double(sRGB(closestInd,:))/255,'filled','MarkerEdgeColor','k')
        daspect([1,1,1])
        view([-75,0])

        xlabel('a*')
        ylabel('b*')
        % view(2)
        title('CIELab')
        axis equal

        sq_er(Lstar,radius) = sum((requestedLocations' - Lab(:,closestInd)).^2,"all");

        figure(6523)
        nexttile
        hold on
        for j = 1:length(closestInd)
            plot(SToWls(S_SPD),SPD(:,closestInd(j)),...
                'Color',double(sRGB(closestInd(j),:))/255)
        end
        axis tight
        NDfilterInd = floor(closestInd/size(data_mod,1))
        filterInd = mod(closestInd,size(data_mod,1)-1);

        data_mod.Var1(filterInd+1)
        data_mod.Var4(filterInd+1)
        data_mod.Var5(filterInd+1)

        legend(data_mod.Var1(filterInd+1))
    end
end

figure, imagesc(sq_er)
colorbar
axis equal tight
caxis([0,1000])
xlabel('radius')
ylabel('L*')
set(gca,'YDir','normal')

hold on
scatter([sqrt((Lab(2,whichFiltersInd).^2)+(Lab(3,whichFiltersInd).^2))],...
    Lab(1,whichFiltersInd),...
    [],double(sRGB(whichFiltersInd,:))/255,'filled');

% scatter(repelem([32,42,52],3),repmat([55,65,75],1,3),'k*')

% Which filter(s) are those?

figure, hold on
for j = 1:length(closestInd)
    i = closestInd(j);
    % plot(SToWls(S_SPD),SPD(:,i),'Color',double(sRGB(i,:))/255)
    plot3(SToWls(S_SPD),SPD(:,i),ones(size(SPD(:,i)))*i,...
        'Color',double(sRGB(i,:))/255)
    % plot3(SToWls(S_SPD),SPD(:,mod(i,size(data_mod,1))),ones(size(SPD(:,i)))*i,...
    %     'Color',double(sRGB(i,:))/255)
end
xlabel('Wavelength (nm)')
ylabel('Transmisson')
axis tight
title('SPD')

NDfilterInd = floor(closestInd/size(data_mod,1))
filterInd = mod(closestInd,size(data_mod,1)-1);

data_mod.Var1(filterInd+1)
data_mod.Var4(filterInd+1)
data_mod.Var5(filterInd+1)

legend(data_mod.Var1(filterInd+1))

%% Screen measurements

screenSPD = readmatrix("../displayCharacterization/measurements/data_init_0.csv");
S_screenSPD = [384,4,100]; % This is slightly wacky - it should have a value at 380, but it seems to get lost somewhere along the conversion... (TODO Look into this)

% Convert to XYZ

screenSPDint = SplineSpd(S_screenSPD,screenSPD,S_xyz1931); % should this be SplineSpd/SplineSrf/SplineRaw (what's the difference?)
% TODO Flagging that this is an opportunity to get diverging values between the
% MATLAB vs python pipeline, that I might want to double check

% figure, hold on
% plot(SToWls(S_screenSPD),screenSPD,'k:');
% plot(SToWls(S_xyz1931),screenSPDint,'r:');
% xlabel('Wavelength (nm)')
% ylabel('Transmisson')

screenXYZ = T_xyz1931*screenSPDint;
screenxyY = XYZToxyY(screenXYZ); % Note, Y is not normalised here

% Convert to Luv

screenLuv = XYZToLuv(screenXYZ,testingRoomWall_XYZ); % taking the clear filter as the white point - right decision? (TODO)

figure(420)
scatter3(screenLuv(2,:),screenLuv(3,:),screenLuv(1,:),...
    [],'filled','r','MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.7)
% zlim([60,70])

figure(421)
scatter3(screenxyY(1,:),screenxyY(2,:),screenxyY(3,:),...
    [],'filled','r','MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.7)

% %% Real measurements comparison
% 
% d = dir('SpectralMeasurement231019*.mat');
% 
% for i = 1:length(d)
%     t = load(d(i).name,"SPD"); % doing it this way so I don't overwrite "SPD"
%     filterMeasurements_SPD(i,:) = t.SPD;
%     % 
% end
% 
% figure,
% plot(SToWls(S_SPD),SPD,'k')
% axis tight
% 
% whichFilterMeasurements = ...
%     [34:36;...
%     37:39;...
%     40:42;...
%     43:45;...
%     46:48;...
%     25:27;...
%     28:30;...
%     49:51;...
%     52:54];
% 
% % whichFilterMeasurements = reshape(1:54,3,[])';
% 
% %21 last ones
% %plus blue and purple from before
% 
% % figure, hold on
% % for i = 1:size(whichFilterMeasurements,1)
% %     plot(SToWls(S_SPD),...
% %         filterMeasurements_SPD(whichFilterMeasurements(i,:),:),...
% %         'Color',double(sRGB(closestInd(i),:))/255)
% % end
% 
% % figure, hold on
% % for i = 1:size(whichFilterMeasurements,1)
% %     plot(SToWls(S_SPD),...
% %         filterMeasurements_SPD(whichFilterMeasurements(i,:),:),'k')
% % end
% 
% filterMeasurements_SPDint = SplineSpd(S_SPD,filterMeasurements_SPD',S_xyz1931); % should this be SplineSpd/SplineSrf/SplineRaw (what's the difference?)
% filterMeasurements_XYZ = T_xyz1931*filterMeasurements_SPDint;
% filterMeasurements_Luv = XYZToLuv(filterMeasurements_XYZ,testingRoomWall_XYZ); 
% 
% figure(422)
% % for i = 1:size(whichFilterMeasurements,1)
% %     scatter3(filterMeasurements_Luv(2,whichFilterMeasurements(i,:)),...
% %         filterMeasurements_Luv(3,whichFilterMeasurements(i,:)),...
% %         filterMeasurements_Luv(1,whichFilterMeasurements(i,:)),...
% %         [],double(sRGB(closestInd(i),:))/255,'filled')
% % end
% 
% % figure,hold on
% for i = 1:size(whichFilterMeasurements,1)
%     scatter3(filterMeasurements_Luv(2,whichFilterMeasurements(i,:)),...
%         filterMeasurements_Luv(3,whichFilterMeasurements(i,:)),...
%         filterMeasurements_Luv(1,whichFilterMeasurements(i,:)),...
%         [],'k','filled')
%     text(filterMeasurements_Luv(2,whichFilterMeasurements(i,1))+5,...
%         filterMeasurements_Luv(3,whichFilterMeasurements(i,1)),...
%         filterMeasurements_Luv(1,whichFilterMeasurements(i,1)),...
%         num2str(i))
% end
% axis equal
% daspect([1,1,1])


