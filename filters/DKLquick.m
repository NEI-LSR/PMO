clear, clc, close all

load('DKLsubset.mat', 'DKLsubset')
size(DKLsubset)

DKL_L3S1 = DKLsubset(:,:,73,32);
DKL_L3S3 = DKLsubset(:,:,73,52);
DKL_L1S1 = DKLsubset(:,:,53,32);

figure, scatter3(DKL_L3S1(2,:),DKL_L3S1(3,:),DKL_L3S1(1,:))
title('L3S1')
% saturation looks a little higher than anticipated for this one, compared
% to the others...

figure, scatter3(DKL_L3S3(2,:),DKL_L3S3(3,:),DKL_L3S3(1,:))
title('L3S3')

figure, scatter3(DKL_L1S1(2,:),DKL_L1S1(3,:),DKL_L1S1(1,:))
title('L1S1')

%%

DKL_hue_angles(1,:) = rad2deg(cart2pol(DKL_L3S1(2,:),DKL_L3S1(3,:))); 
DKL_hue_angles(2,:) = rad2deg(cart2pol(DKL_L3S3(2,:),DKL_L3S3(3,:))); 
DKL_hue_angles(3,:) = rad2deg(cart2pol(DKL_L1S1(2,:),DKL_L1S1(3,:))); 
DKL_hue_angles(DKL_hue_angles<0) = 360 + DKL_hue_angles(DKL_hue_angles<0);

figure, plot(DKL_hue_angles')
ylim([0,360])