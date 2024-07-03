%% This script provides guidance for the use of the program.
clear;clc;close all;

%% modal identification

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

%% make seismic responses if necessary
load('Synthetic data\two-story shear-type building structure\2SB_resp_Lp.mat','fdata','fs'); % load new ground motion

nsample = 200;
Nmode = length(out.f);

% mean predicted structural resposnes
in.f = out.f;
in.z = out.z;
in.phi = out.Phi;
in.Gamma = out.Gamma;

in.fdata = fdata;
in.fs = fs;
in.eplison_k = normrnd(0,mean(out.coefv.Se)*mean(out.Se));

mean_pred_tdata = Pred_freq(in); % predicted tdata (structural responses)

% credible interval
for jj = 1:nsample
    in.f = out.pred.f(jj,:);
    in.z = out.pred.z(jj,:);
    in.phi = reshape(out.pred.phi(jj,:)',[],Nmode);
    in.Gamma = reshape(out.pred.Gamma(jj,:)',[],Nmode);
    in.eplison_k = normrnd(0,mean(out.coefv.Se)*mean(out.Se));

    pred_tdata(:,:,jj) = Pred_freq(in); % one of predicted tdata
end

%% For more information, please refer to the comments in 'BayEma_SE.m'

