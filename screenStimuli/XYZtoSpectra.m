clear, clc, close all

%% Load surface reflectances

load sur_vrhel
% first dim is wavelength

[coeff,score,latent,tsquared,explained,mu] = pca(sur_vrhel',...
    'NumComponents',50,...
    'Centered',false);

figure,
plot(SToWls(S_vrhel),coeff)
axis tight

%% Compute XYZ

load T_xyz1931.mat T_xyz1931 S_xyz1931 % Requires PsychToolbox

sur_vrhel_int = SplineSpd(S_vrhel,sur_vrhel,S_xyz1931,1);

% % checking that the extrapolation isn't a big problem
% figure, hold on
% plot(SToWls(S_xyz1931),sur_vrhel_int)
% plot(SToWls(S_xyz1931),T_xyz1931/4,'k','LineWidth',2)
% axis tight

XYZ = T_xyz1931*sur_vrhel_int;

%% Compute conversion matrix from XYZ to `score`

XYZtoScore = XYZ'\score;

%% Compute spectra from XYZ

reconScore = XYZtoScore'*XYZ;
reconSPD = reconScore'*coeff';

for i = 1%:10:170
    figure, hold on
    plot(sur_vrhel(:,i))
    plot(reconSPD(i,:))
    legend
end

figure, hold on
plot(SToWls(S_vrhel),reconSPD([81,82,84,85,86,87,88,91,92:95,98,100,101,103:108,110:117],:)')
axis tight

