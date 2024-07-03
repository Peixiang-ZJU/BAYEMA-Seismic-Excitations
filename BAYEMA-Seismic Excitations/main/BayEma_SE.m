function out = BayEma_SE(in)
% *************************************************************************
%
%         BAYesian Modal Analysis (Known Seismic excitations)
%     ================================================================
%             Peixiang WANG |Zhejiang University| June 2024
%     ================================================================
% If you use this code, please cite:
% Wang P, Li B, Zhang F, Chen X, Ni Y. Fast Bayesian modal identification with known
% seismic excitations. Earthquake Engng Struct Dyn. 2024;1-30. https://doi.org/10.1002/eqe.4181
% *************************************************************************
%
% OUT = BayEma_EQ(IN)
%
% IN contains the following mandatory fields:
%    tdata = (nt,nch) measured structural responses, relative to the ground;
%            nt = no. of sample points;
%            nch = total no. of measured dofs
%    fdata = (nt,nch) measured ground acceleration
%            nt = no. of time points;
%            nch = no. of seismic directions
%       fs = scalar, sampling rate (Hz)
%     f1f2 = (nb,2), f1f2(i,1) and f1f2(i,2) give the lower and upper bound
%            frequency of the i-th band; nb = no. of bands
%       f0 = {nb,1} cell array, f0{i} = 1-D array of initial guesses of natural
%            frequency in the i-th band
%
% The following optional fields can be supplied in IN:
%    iband = 1-D array of indices of bands to be considered
%            If iband=[] or not proivded, all bands will be considered
%    options(1) = 0: do NOT make responses predction
%               = 1: make responses predction
%    nsample = the number of sampling for response prediction
%
% OUT is a structure with the following fields:
%
%    f = (1,m) mpv of natural frequencies (Hz), m = no. of modes
%    z = (1,m) mpv of damping ratios
%    Se = (1,m) mpv of PSD of prediction error (e.g., channel noise)
%    phi = (n,m) mpv of mode shape, each column gives a mode shape with unit noramlisation
%    Gamma = (d,m) modal participation factor, d = no. of excitation direction
%    coefv = a structure containing the posterior coefficient of variation
%            (c.o.v.) of f, z, Se, phi, Gamma. (Unit: 1)
%    pred = a structure containing the modal parameters obtained from
%    sampling
%
% *************************************************************************
%**************************************************************************

% copyright statement
disp('============================================');
disp('BAYEMA-SE --- BAYesian Experimental Modal Analysis - Seismic Excitations');
disp('Written by Peixiang Wang, ZJU, 2024');
disp('============================================');


in_default.maxiter = 1000;
in_default.tol_cvg = 1e-6;
in_default.iband = [];     
in_default.alg = 'Direct';
in_default.myoptions(1) = 0; 
in_default.nsample = 200;


in = optfield(in,in_default);

switch in.alg
    case 'Direct'
        out = BayEma_SE_proc(in,'BayEma_SE_direct');  % direct optimization 
end
