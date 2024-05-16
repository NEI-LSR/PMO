clc, clear, close all

T = readtable('live_stim_data.xlsx','Sheet','All resuts (V1 V2 faces)');
% T = readtable('combined_live_stim_data.csv');

hue_angle = [135,180,225,270,315,0,45,90];

L3S1 = table2array(T(2:end-1,3));
L3S3 = table2array(T(2:end-1,8));
L1S1 = table2array(T(2:end-1,13));

%% CIELAB

figure,
polarplot([deg2rad(hue_angle),deg2rad(hue_angle(1))],...
    [L3S1',L3S1(1)],'ko-','DisplayName','L3S1')
hold on

polarplot([deg2rad(hue_angle),deg2rad(hue_angle(1))],...
    [L3S3',L3S3(1)],'ko:','DisplayName','L3S3')

polarplot([deg2rad(hue_angle),deg2rad(hue_angle(1))],...
    [L1S1',L1S1(1)],'ko--','DisplayName','L1S1')

legend('AutoUpdate','off')

% polarplot(linspace(0,2*pi,1000),ones(1000,1)*3,'Color',[0,0,0,0.4])

thetaticks(0:45:315)

title('Report as a function of CIELAB hue angle')

% TODO What does it look like with the _actual_ hue angles rather than the
% target hue angles?

%% DKL

load('DKLsubset.mat', 'DKLsubset')

DKL_L3S1 = DKLsubset(:,:,73,32);
DKL_L3S3 = DKLsubset(:,:,73,52);
DKL_L1S1 = DKLsubset(:,:,53,32);

% DKL_L3S1 = [-0.856476597173577	-0.937278812544493	-0.928213781150367	-0.895139653946265	-1.03364836782582	-0.852719027448376	-1.11053102337354	-1.12948498155823	-0.940685221634920
% 0.0578819018756198	0.0589601427455918	0.0275326482784112	-0.0183270573549438	-0.0377142960052881	-0.0557364967490240	-0.0229625523307020	0.0124652479093524	-0.00527754014345971
% 0.0145238348180642	-0.157965476167355	-0.224799066066804	-0.149629553427485	0.0468584652941974	0.223811291554045	0.313042947659699	0.228846737723592	0.0350165912946425];
% 
% DKL_L3S3 = [-1.01543454521201	-1.14122934475943	-0.964190637406599	-0.736603796208872	-0.955169996782731	-1.12780220530950	-1.11577450917705	-1.01492927498581	-0.940685221634920
% 0.0920627928123777	0.0785085516869831	0.0307150157235090	-0.0403752115595957	-0.0801942521710604	-0.0774678524089014	-0.0408465829782501	0.0444155706306907	-0.00527754014345971
% -0.0585375192990609	-0.204682051114546	-0.335776966905387	-0.303823814213230	0.0216131592738020	0.414575029140712	0.398774373763702	0.207401538461126	0.0350165912946425];
% 
% DKL_L1S1 = [-1.32359687952135	-1.31483940975877	-1.41080008825850	-1.34739521429362	-1.42947444191624	-1.41966912622512	-1.36215985785706	-1.33035028650727	-1.32165660248284
% 0.0421852825649517	0.0292454767323462	0.00812251922556773	-0.0177991246206018	-0.0326326640512683	-0.0335070866889083	-0.0162809928997107	0.00893044520964481	-0.00242986842901338
% 0.00837468178722843	-0.0874327984709025	-0.117986559942733	-0.114656725721755	0.0103926879515793	0.131457499323839	0.212310369007810	0.123180158432990	-0.00598604857583987];

DKL_hue_angles(1,:) = rad2deg(cart2pol(DKL_L3S1(2,:),DKL_L3S1(3,:))); 
DKL_hue_angles(2,:) = rad2deg(cart2pol(DKL_L3S3(2,:),DKL_L3S3(3,:))); 
DKL_hue_angles(3,:) = rad2deg(cart2pol(DKL_L1S1(2,:),DKL_L1S1(3,:))); 
DKL_hue_angles(DKL_hue_angles<0) = 360 + DKL_hue_angles(DKL_hue_angles<0);

% DKL_hue_angles(:,1) = [];

% figure, plot(0:45:315,DKL_hue_angles(1:end-1))
figure, hold on
plot(DKL_hue_angles(1,1:end-1),'o-k')
plot(DKL_hue_angles(2,1:end-1),'o:k')
plot(DKL_hue_angles(3,1:end-1),'o--k')

figure,
polarplot([deg2rad(DKL_hue_angles(1,1:end-1)),deg2rad(DKL_hue_angles(1,1))],...
    [L3S1([6:8,1:5])',L3S1(6)],'k-o','DisplayName','L3S1')
hold on

polarplot([deg2rad(DKL_hue_angles(2,1:end-1)),deg2rad(DKL_hue_angles(2,1))],...
    [L3S3([6:8,1:5])',L3S3(6)],'k:o','DisplayName','L3S3')

polarplot([deg2rad(DKL_hue_angles(3,1:end-1)),deg2rad(DKL_hue_angles(3,1))],...
    [L1S1([6:8,1:5])',L1S1(6)],'k--o','DisplayName','L1S1')

legend('AutoUpdate','off')

% polarplot(linspace(0,2*pi,1000),ones(1000,1)*3,'Color',[0,0,0,0.4])

thetaticks(0:45:315)

title('Report as a function of DKL hue angle')


