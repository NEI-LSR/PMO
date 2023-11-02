clear, clc, close all

measure = 1;

if measure
    retval = PR655init('COM7');
    S_SPD = [380,1,401];

    load('C:\Users\cege-user\Documents\PMO\screenStimuli\gunVals.mat');

    figure, hold on
    for i = 34:44%1:size(gunVals,2)
        im = permute(repmat(gunVals(:,i),[1,1000,1000]),[3,2,1]);
        imshow(uint8(im))
        drawnow
        pause(0.1)
        % SPD = PR670measspd(S_SPD);
        SPD(:,i) = PR655measspd(S_SPD);
    end

    save(['SpectralMeasurement',char(datetime('now','format','yyMMdd-HHmmss'))]);
end

%% Compute XYZ

if ~measure
    load('SpectralMeasurement231030-113448.mat')
end

load T_xyz1931.mat T_xyz1931 S_xyz1931 % Requires PsychToolbox
screenSPD = reshape(SPD,[401,11,4]);
S_screenSPD = S_SPD;

% Convert to XYZ

for i = 1:4
    screenSPDint(:,:,i) = SplineSpd(S_screenSPD,screenSPD(:,:,i),S_xyz1931); % should this be SplineSpd/SplineSrf/SplineRaw (what's the difference?)
end

for i = 1:4
    screenXYZ(:,:,i) = T_xyz1931*screenSPDint(:,:,i);
end

load('screenXYZNormFactor.mat')
screenXYZ = (screenXYZ/screenXYZNormFactor); % pulled from screen calibration analsis (Y of white)

cols = {'r','g','b'};
names = {'X','Y','Z','all'};

% figure, hold on
% tiledlayout(2,2)
% for i = 1:4
%     nexttile
%     for j = 1:3
%         scatter(0:0.1:1,screenXYZ(j,:,i),cols{j})
%         hold on
%     end
%     title(names{i})
% 
% end

figure, hold on
for j = 1:3
    plot(0:0.1:1,screenXYZ(j,:,i),'-o','Color',cols{j})
    hold on
end
plot(0:0.1:1,0:0.1:1,'k')
% plot(0:0.1:1,[0.1:0.1:1,1],'k')
axis equal tight


% screenXYZ = (screenXYZ/screenXYZ(2,end,end));


