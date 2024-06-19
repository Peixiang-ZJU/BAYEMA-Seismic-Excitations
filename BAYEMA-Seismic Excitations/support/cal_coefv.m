function coefv = cal_coefv(VarMat,Nmode,Ndof,nd,Gamma_norm)
% Function to calculate the c.o.v.s of modal parameters of interest

% % c.o.v of frequency
var_f = diag( VarMat(1:2:2*Nmode,1:2:2*Nmode) )';
coefv.f = sqrt(var_f);


% % c.o.v of damping ratio
var_z = diag( VarMat(2:2:2*Nmode,2:2:2*Nmode) )';
coefv.z = sqrt(var_z);

% % c.o.v of mode shapes
for ii = 1:Nmode
    var_phi = VarMat(2*Nmode+1:2*Nmode+Nmode*Ndof,2*Nmode+1:2*Nmode+Nmode*Ndof);
    coefv.phi(ii) = sqrt( trace(var_phi((ii-1)*Ndof+1:ii*Ndof,(ii-1)*Ndof+1:ii*Ndof) ) );
end

% % c.o.v of modal participation
for ii = 1:Nmode
    var_Gamma = VarMat(2*Nmode+Nmode*Ndof+1:2*Nmode+Nmode*Ndof+Nmode*nd,2*Nmode+Nmode*Ndof+1:2*Nmode+Nmode*Ndof+Nmode*nd);
end
for jj = 1:nd
    coefv.gamma(jj) = sqrt( var_Gamma(jj,jj) ./ Gamma_norm(jj,1)^2 );
end

% % c.o.v of PSD of channels noise
var_Se = diag( VarMat(2*Nmode+Nmode*Ndof+Nmode*nd+1,2*Nmode+Nmode*Ndof+Nmode*nd+1) )';
coefv.Se = sqrt(var_Se);

