clear, clc, close all

%% Load data

data = readtable('C:\Users\cege-user\Documents\PMO\filters\CORRECT ORDER GELS MEASUREMENT SPECTRUMS 3200 5600 - DBR V1 - 21 12 2020.xlsx');
% This by default loads the first sheet only (which is the 3200K measurements)

%% Normalise and plot SPDs

SPD_raw = table2array(data(2:end,13:end));
S_SPD = WlsToS(table2array(data(1,13:end))');

SPD = (SPD_raw./SPD_raw(195,:))'; %normalisation by "clear"

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

%%

load T_xyz1931.mat T_xyz1931 S_xyz1931 % Requires PsychToolbox

SPDint = SplineSpd(S_SPD,SPD,S_xyz1931); % should this be SplineSpd/SplineSrf/SplineRaw (what's the difference?)

% figure, 
% % plot(SToWls(S_SPD),SPD);
% plot(SToWls(S_xyz1931),SPDint(:,1:100)); % subset, easier to see
% xlabel('Wavelength (nm)')
% ylabel('Transmisson')

XYZ = T_xyz1931*SPDint;
xyY = XYZToxyY(XYZ);

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

% print(gcf,'-vector','-dsvg',['1931xy','.svg'])

figure, hold on
for i = 1:size(SPD,2)
    % plot(SToWls(S_SPD),SPD(:,i),'Color',double(sRGB(i,:))/255)
    plot3(SToWls(S_SPD),SPD(:,i),ones(size(SPD(:,i)))*i,...
        'Color',double(sRGB(i,:))/255)
end
xlabel('Wavelength (nm)')
ylabel('Transmisson')
axis tight
title('CIE xy')

% print('SPD.png','-dpng')

%% u'v'

upvp = xyTouv(xyY(1:2,:));

figure,
DrawChromaticity('upvp')
scatter3(upvp(1,:),upvp(2,:),xyY(3,:),...
    [],double(sRGB)/255,'filled','MarkerEdgeColor','k')
daspect([1,1,200])
xlim([0,0.6])
ylim([0,0.6])
zlabel('Y_{1931}')
view(2)
title('CIE upvp')

% print(gcf,'-vector','-dsvg',['upvp','.svg'])

%% CIELUV

Luv = XYZToLuv(XYZ,XYZ(:,195)); % taking the clear filter as the white point - right decision?
Lab = XYZToLab(XYZ,XYZ(:,195));

figure,
scatter3(Luv(2,:),Luv(3,:),Luv(1,:),...
    [],double(sRGB)/255,'filled','MarkerEdgeColor','k')
% zlim([60,70])
daspect([1,1,1])
xlabel('u''')
ylabel('v''')
view(2)
title('CIELuv')

zlim([45,55])

print(gcf,'-vector','-dsvg',['CIELuv','.svg'])


figure,
scatter3(Lab(2,:),Lab(3,:),Lab(1,:),...
    [],double(sRGB)/255,'filled','MarkerEdgeColor','k')
% zlim([60,70])
daspect([1,1,1])
xlabel('a*')
ylabel('b*')
view(2)
title('CIELab')

% print(gcf,'-vector','-dsvg',['CIELab','.svg'])

%% DKL

load T_cones_sp.mat T_cones_sp S_cones_sp
load T_xyzJuddVos.mat T_xyzJuddVos S_xyzJuddVos
T_Y = 683*T_xyzJuddVos(2,:);
T_Y = SplineCmf(S_xyzJuddVos,T_Y,S_cones_sp);

LMS     = T_cones_sp * SplineSpd(S_SPD,SPD,S_cones_sp);
bgLMS   = mean(LMS,2); % !!!!!!!!!!!
LMSinc = LMS - bgLMS;

[M_ConeIncToDKL,LMLumWeights] = ComputeDKL_M(bgLMS,T_cones_sp,T_Y);

DKL = M_ConeIncToDKL*LMS;

figure,
scatter3(DKL(2,:),DKL(3,:),DKL(1,:),...
    [],double(sRGB)/255,'filled','MarkerEdgeColor','k')
xlabel('L-M')
ylabel('S-(L+M)')
zlabel('L+M')
title('DKL')
view(2)

% note: relative axis scaling is a choice

% print(gcf,'-vector','-dsvg',['DKL','.svg'])

%%
% 
% figure
% scatter3(Luv(2,:),Luv(3,:),Luv(1,:),...
%     [],double(sRGB)/255,'filled','MarkerEdgeColor','k')
% % zlim([60,70])
% daspect([1,1,1])
% xlabel('u''')
% ylabel('v''')
% view(2)
% title('CIELuv')
% 
% % Add a neutral density filter to all the colors,
% % but keep the white point the same...
% SPDint_nd = SPDint.*SPDint(:,530); % 530 is "neutral gray"
% XYZ_nd = T_xyz1931*SPDint_nd;
% Luv_nd = XYZToLuv(XYZ_nd,XYZ(:,195)); % taking the clear filter as the white point 
% sRGBlin_nd = XYZToSRGBPrimary(XYZ_nd./max(XYZ(2,:)));
% sRGB_nd = uint8(SRGBGammaCorrect(sRGBlin_nd,0)');
% 
% hold on
% 
% scatter3(Luv_nd(2,:),Luv_nd(3,:),Luv_nd(1,:),...
%     [],double(sRGB_nd)/255,'filled','MarkerEdgeColor','w')
% zlim([40,70])
% 
% % And another...
% SPDint_nd2 = SPDint.*SPDint(:,530).*SPDint(:,530); % 530 is "neutral gray"
% XYZ_nd2 = T_xyz1931*SPDint_nd2;
% Luv_nd2 = XYZToLuv(XYZ_nd2,XYZ(:,195)); % taking the clear filter as the white point 
% sRGBlin_nd2 = XYZToSRGBPrimary(XYZ_nd2./max(XYZ(2,:)));
% sRGB_nd2 = uint8(SRGBGammaCorrect(sRGBlin_nd2,0)');
% 
% scatter3(Luv_nd2(2,:),Luv_nd2(3,:),Luv_nd2(1,:),...
%     [],double(sRGB_nd2)/255,'filled','MarkerEdgeColor','r')
% zlim([45,55])
% 
% set(gca,'Color',[0.5,0.5,0.5])