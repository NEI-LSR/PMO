clc, clear, close all

% Take measurements with the PR670

% Spectral Resolution: 3.12 nm / pixel

%%

retval = PR670init('COM7');

%%
S_SPD = [380,1,401]; 
SPD = PR670measspd(S_SPD);

%%
figure, hold on
plot(SToWls(S_SPD),SPD)

save(['SpectralMeasurement',char(datetime('now','format','yyMMdd-HHmmss'))],...
    'S_SPD','SPD');


