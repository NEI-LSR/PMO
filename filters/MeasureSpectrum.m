clc, clear, close all

% Take measurements w ith the PR655/670

%%

% retval = PR670init('COM7');
retval = PR655init('COM7');

%%
S_SPD = [380,1,401]; 
% SPD = PR670measspd(S_SPD);
SPD = PR655measspd(S_SPD);

%%
% figure, hold on
% plot(SToWls(S_SPD),SPD)

save(['SpectralMeasurement',char(datetime('now','format','yyMMdd-HHmmss'))],...
    'S_SPD','SPD');


