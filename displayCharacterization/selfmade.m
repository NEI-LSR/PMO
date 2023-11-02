clear, clc, close all

retval = PR655init('COM7');
S_SPD = [380,1,401];

steps = 0:15:255;
% steps = 255;

figure, hold on
for j = 1:3
    for i = 1:length(steps)
        im = zeros(1000,1000,3);
        im(:,:,j) = ones(1000)*steps(i);
        imshow(uint8(im))
        drawnow        
        % SPD = PR670measspd(S_SPD);
        SPD(:,i,j) = PR655measspd(S_SPD);
    end
end

figure, hold on
for i = 1:length(steps)
    im = ones(1000,1000,3)*steps(i);
    imshow(uint8(im))
    drawnow
    % SPD = PR670measspd(S_SPD);
    SPD(:,i,4) = PR655measspd(S_SPD);
end

%%
cols = {'r','g','b','k'};
figure, hold on
for i = 1:4
    plot(SPD(:,:,i),'Color',cols{i})
end

%%

save(['SpectralMeasurement',char(datetime('now','format','yyMMdd-HHmmss'))]);


