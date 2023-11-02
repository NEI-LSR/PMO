

folder2 = 'C:\Users\cege-user\Documents\PMO\screenStimuli\images\filteredImages\';
d2 = dir([folder2,'IMG*.tif']);

for i = 1:length(d2)
        ind(i) = length(d2(i).name) > 30;
    % a(i) = length(d2(i).name);
end

d2(ind) = [];

% figure, 
% histogram(a)

% Then copy paste into excel and sort in excel