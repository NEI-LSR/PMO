clear, clc, close all

load("SpectralMeasurement231030-120356.mat",'SPD','S_SPD');

% SPD = reshape(SPD,[401,11,4])

cols = {'r','g','b','k'};
figure, hold on
for i = 1:4
    plot(SPD(:,:,i),'Color',cols{i})
end

%%

figure, hold on
for i = 1:18
plot(SToWls(S_SPD),SPD(:,i,4),'k','DisplayName','white')
plot(SToWls(S_SPD),SPD(:,i,1)+SPD(:,i,2)+SPD(:,i,3),'DisplayName','all added')
end
%legend('Location','best')
axis tight