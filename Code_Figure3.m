clc; clear all;

%%  CODE TO PLOT FIGURE 3
%
%   Whole-brain multimodal neuroimaging model using serotonin receptor maps explain non-linear functional effects of LSD
%   Deco,G., Cruzat,J., Cabral, J., Knudsen,G.M., Carhart-Harris,R.L., Whybrow,P.C., 
%       Logothetis,N.K. & Kringelbach,M.L. (2018) Current Biology
%
%
%   Written by Josephine Cruzat josephine.cruzat@upf.edu
%   October 2017, Barcelona-Spain
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DEFINE COLORS

lineProps1.col{1} = 'b'; 
lineProps2.col{1} = [0 0.5 0];
lineProps3.col{1} = [1 0.84 0];
lineProps4.col{1} = 'r'; 


%% FIGURE 3A: Whole-brain fitting of the placebo condition 

load( 'fneuro.mat' );

mFCDfitt5   = mean(FCDfitt5,2);
stdFCDfitt5 = std(FCDfitt5,[],2);
mfitting5   = mean(fitting5,2);
stdfitting5 = std(fitting5,[],2);

figure
mseb(WE,mFCDfitt5',stdFCDfitt5',lineProps2)
hold on;
mseb(WE,mfitting5',stdfitting5',lineProps4)
legend('FCD placebo','FC Placebo')
xlabel('Global Coupling')
ylabel('Fitting')


%% FIGURE 3B: Neuronal gain for the placebo and LSD condition 

load( 'fgain_plbo.mat' )
load( 'fgain_LSD.mat' )

mFCDfitt_LSD    = mean(FCDfitt_LSD,2);
stdFCDfitt_LSD  = std(FCDfitt_LSD,[],2);
mFCDfitt_plbo   = mean(FCDfitt_plbo,2);
stdFCDfitt_plbo = std(FCDfitt_plbo,[],2);

figure
mseb(WG,mFCDfitt_plbo',stdFCDfitt_plbo',lineProps2)
hold on;
mseb(WG,mFCDfitt_LSD',stdFCDfitt_LSD',lineProps1)
legend('Placebo','LSD')
xlabel('Excitatory Gain Modulation')
ylabel('FCD Fitting')


%% FIGURE 3C: Reshuffling the 5-HT2A receptor densities randomly

load( 'fgain_rnd.mat' )

vec_FCDfitt_rnd = FCDfitt_rnd(100,:);     % 100 is the position of the min value
vec_FCDfitt_LSD = FCDfitt_LSD(100,:);

[pma hma] = ranksum(vec_FCDfitt_rnd,vec_FCDfitt_LSD); %i del wge LSD

figure
boxplot([vec_FCDfitt_LSD',vec_FCDfitt_rnd'],'Notch','on','Labels',{'LSD','Reshuffled'},'Whisker',1)

%% FIGURE 3D: Specificity of other receptor binding maps

load( 'fgain_rcp1A.mat' )
load( 'fgain_rcp1B.mat' )
load( 'fgain_rcpT4.mat' )
load( 'fgain_rcpTT.mat' )

% plot(mean(FCDfitt_LSD,2))         % To find the min

LSD     = FCDfitt_LSD(100,:);       % 100 is the position of the min value of the FCDfitt_LSD
random  = FCDfitt_rnd(100,:);
uniform = FCDfitt_rnd(1,:);
rcp1A   = FCDfitt_rcp1A(100,:);
rcp1B   = FCDfitt_rcp1B(100,:);
rcpTT   = FCDfitt_rcpTT(100,:);
rcpT4   = FCDfitt_rcpT4(100,:);


load( 'LSD_training.mat' )
% % To find the min
% plot(mean(FCDfitt_train,2)) 
LSD_training = FCDfitt_train(92,:); % 92 is the position of the minimun

load( 'LSD_test.mat' )
LSD_test = FCDfitt_test(92,:);

% P VALUES

[pma_uniform hma] = ranksum(uniform,LSD);
[pma_recp1A hma]  = ranksum(rcp1A,LSD); 
[pma_recp1B hma]  = ranksum(rcp1B,LSD); 
[pma_recpT4 hma]  = ranksum(rcpT4,LSD); 
[pma_recpTT hma]  = ranksum(rcpTT,LSD); 
[pma_LSD_test hma]= ranksum(LSD_test,LSD); 
[pma_LSD_training hma] = ranksum(LSD_training,LSD); 

figure
boxplot([LSD',LSD_training',LSD_test',rcp1A',rcp1B',rcpT4',rcpTT',uniform'],'Notch','on','Labels',...
    {'LSD','LSD training','LSD test','5HT1A','5HT1B','5HT4','5HTT','Uniform'},'Whisker',1)
ylabel('FCD Kolmogorov-Smirnov Distance')


