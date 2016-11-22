% run the 3-D Questionnaire on flexible trees using several speakers saying
% selected words

%% Initialization
close all;
clear;

addpath(genpath('./3D_Questionnaire'));

%% perm data

load DataTimeVaryingSpring.mat

%% Run Questionnaire 3D
params = SetQuestParamsPen(ndims(data), false);
[ Trees, dual_aff, init_aff, embedding ] = RunGenericDimsQuestionnaire( params, data );

%%
[vD1, vD2, vD3, vD4] = size(data);
% PlotTreesAndEmbedding(dual_aff, Trees, {1:vD1, 1:vD2, 1:vD3}, ....
%                       {'Frame', 'Time', 'Params'});

%%
mPhi = embedding{2};

%%
f1 = sqrt(g / L) / (2 * pi);
f2 = sqrt( (2 * L * vK + g * m) / (L * m) ) / (2 * pi);
 
%%
figure; hold on;
Win = 60;
spectrogram(mPhi(:,1), blackman(Win), Win - 1, 512, Fs, 'yaxis');
colormap gray
colormap(flipud(colormap))
 
%%
xLim = [0.5 4.5];
yLim = [0, 30];
vF = sqrt( (2 * L * vK + g * m) / (L * m) ) / (2 * pi);
plot(T, vF, '--y', 'LineWidth', 1); grid on;
plot(T, f1*ones(1,length(T)), '--r', 'LineWidth', 1); grid on;
axis([xLim, yLim]);
