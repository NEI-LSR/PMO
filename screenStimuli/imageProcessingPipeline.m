clc, clear, close all

%%

sRGBinput = imread('_MG_9169.tif');
% figure,imshow(sRGBinput)
sRGBreshape = reshape(sRGBinput,[size(sRGBinput,1)*size(sRGBinput,2),size(sRGBinput,3)]);

sRGBlin = SRGBGammaUncorrect(sRGBreshape);
XYZ = SRGBPrimaryToXYZ(sRGBlin');

xy = XYZToxyY(XYZ);

sRGBgamut_sRGB = [1,0,0;0,1,0;0,0,1];
sRGBgamut_sRGBlin = SRGBGammaUncorrect(sRGBgamut_sRGB);
sRGBgamut_XYZ = SRGBPrimaryToXYZ(sRGBgamut_sRGBlin');
sRGBgamut_xy = XYZToxyY(sRGBgamut_XYZ);

figure, hold on
DrawChromaticity
% scatter(xy(1,1:100:end),xy(2,1:100:end),'k','filled',...
%     'MarkerEdgeAlpha',0.1,'MarkerFaceAlpha',0.1)
scatter3(xy(1,1:100:end),xy(2,1:100:end),xy(3,1:100:end),...
    [],double(sRGBreshape(1:100:end,:))/255,'filled')
scatter(sRGBgamut_xy(1,:),sRGBgamut_xy(2,:),'r*')
plot([sRGBgamut_xy(1,:),sRGBgamut_xy(1,1)],[sRGBgamut_xy(2,:),sRGBgamut_xy(2,1)],'r:')

%% Convert from XYZ to display space

% Linearization
