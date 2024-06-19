function Hessian = BayEma_SE_HessAssm(Hess,Nmode,Ndof,nd)
% Function to assembles the Hessian matrix.

% set parameter indices first
Efz = 1:2*Nmode;
EPhi = 2*Nmode+1:2*Nmode+Nmode*Ndof;
EGamma = 2*Nmode+Nmode*Ndof+1:2*Nmode+Nmode*Ndof+Nmode*nd;
ESe = 2*Nmode+Nmode*Ndof+Nmode*nd+1;

% pivot terms
Hessian(Efz,Efz) = Hess.ddL_theta_m;
Hessian(EPhi,EPhi) = Hess.ddL_Phi;
Hessian(EGamma,EGamma) = Hess.ddL_Gamma;
Hessian(ESe,ESe) = Hess.ddL_Se;

% cross terms
Hessian(Efz,EPhi) = Hess.ddL_Phi_theta_m;

Hessian(Efz,ESe) = Hess.ddL_Se_theta_m;

Hessian(Efz,EGamma) = Hess.ddL_Gamma_theta_m;

Hessian(EPhi,EGamma) = Hess.ddL_Phi_Gamma;

Hessian(EPhi,ESe) = Hess.ddL_Phi_Se;

Hessian(EGamma,ESe) = Hess.ddL_Se_Gamma;


Hess = real(Hessian);
Hessian = triu(Hessian) + tril(Hessian.',-1);
