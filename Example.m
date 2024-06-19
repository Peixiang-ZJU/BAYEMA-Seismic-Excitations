%% This script provides guidance for the use of the program.
clear;clc;close all;

% load data
in = load('Synthetic data\two-story shear-type building structure\2SB_resp_El.mat','tdata','fdata','fs');

% set initial guesses and selected freq. band
in.f0 = {0.706
         1.222
         [3.055,3.200]};
in.f1f2 = [0.601	0.841
           0.997	1.498
           2.918	3.396];


% optional fields
in.myoptions(1) = 1;


out = BayEma_SE(in);

% For more information, please refer to the comments in 'BayEma_SE.m'

