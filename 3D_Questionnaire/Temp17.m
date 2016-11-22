close all;
clear;

addpath(genpath('./'));

D = rand(512,512,96);

params = SetQuestParams2(ndims(D), false);
RunGenericDimsQuestionnaire(params, D);
