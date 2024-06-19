function [Hessian] = BayEma_SE_hessian(MPV,tdata_fft_band,fdata_fft_band,freq_band)
% Function to analytically calculate the Hessian matrix

f = MPV.f'; z = MPV.z'; Phi = MPV.Phi;
Gamma = MPV.Gamma; Se = MPV.Se;

Nmode = length(f);
nd = size(fdata_fft_band,2);
[Nf,Ndof] = size(tdata_fft_band);


tdata_fft_row = tdata_fft_band.'; % (Ndof,Nf)
fdata_fft_row = fdata_fft_band.'; % (nd,Nf)

% generate Ld
Ld = zeros(Nmode^2,Nmode);
Ld(1:Nmode+1:end,:) = eye(Nmode);

if Ld'*Ld ~= eye(Nmode)
    error('transpose(Ld)*Ld should be equal to eye(Nmode)')
end

PGamma = Perm(Gamma); % Permutation matrix, vec(Gamma') = PGamma*vec(Gamma)
PPhi = Perm(Phi);   % Permutation matrix, vec(Phi') = PGamma*vec(Phi)

beta = f'./freq_band; % Row vector (Number of modes)
h = ( 1-(beta.^2) - 1i*(2.*z'.*beta) ).^-1; % Row vector (Number of modes)


% set the data in selected band in a page-wise manner;
page_tdata_fft_row(:,1,:) = tdata_fft_row;
page_fdata_fft_row(:,1,:) = fdata_fft_row;

% set some commonly used fields for convenience
iSe = 1/Se; iSe2 = 1/Se.^2; iSe3 = 1/Se.^3;

pffdata = pagemtimes(page_fdata_fft_row,'none',page_fdata_fft_row,'ctranspose');
ptfdata = pagemtimes(page_tdata_fft_row,'none',page_fdata_fft_row,'ctranspose');

for ii = 1:Nf
    H = diag(h(ii,:));
    % ------------------------- Pivot terms -------------------------
    % (1) theta_m (i.e., f and z)
    L2vecHk = -Phi'*( conj(ptfdata(:,:,ii)) - Phi*conj(H)*Gamma'*conj(pffdata(:,:,ii)) )*Gamma;
    L2vecHk = vec(L2vecHk).'*iSe;
    L2vecHkconj = conj(L2vecHk);

    Pre_vec_dHk2theta_m = [fun_h(f,z,freq_band(ii),'f') fun_h(f,z,freq_band(ii),'z')];
    for mm = 1:Nmode
        vec_dHk2theta_m(:,2*(mm-1)+1:2*mm) = Pre_vec_dHk2theta_m(:,mm:Nmode:end);
    end
    vec_dconjHk2theta_m = conj(vec_dHk2theta_m);

    % first term 
    ddL_theta_m_first(:,:,ii) = 2*iSe*real( vec_dHk2theta_m.' * Ld' * kron(Gamma'*pffdata(:,:,ii)*Gamma,Phi'*Phi) * Ld * vec_dconjHk2theta_m );
    

    ddH_ff = [fun_h(f,z,freq_band(ii),'ff')]; % H2ff
    ddH_fz = [fun_h(f,z,freq_band(ii),'fz')]; % H2fz;
    ddH_zf = [fun_h(f,z,freq_band(ii),'zf')]; % H2zf;
    ddH_zz = [fun_h(f,z,freq_band(ii),'zz')]; % H2zz;
    % second term
    for mm = 1:Nmode
        Pre_G(1,1,mm) = [L2vecHk L2vecHkconj] * kron(eye(2),Ld) *  [ddH_ff(:,mm) ; conj(ddH_ff(:,mm))];
        Pre_G(1,2,mm) = [L2vecHk L2vecHkconj] * kron(eye(2),Ld) *  [ddH_fz(:,mm) ; conj(ddH_fz(:,mm))];
        Pre_G(2,1,mm) = [L2vecHk L2vecHkconj] * kron(eye(2),Ld) *  [ddH_zf(:,mm) ; conj(ddH_zf(:,mm))];
        Pre_G(2,2,mm) = [L2vecHk L2vecHkconj] * kron(eye(2),Ld) *  [ddH_zz(:,mm) ; conj(ddH_zz(:,mm))];

        G(2*(mm-1)+1:2*mm,2*(mm-1)+1:2*mm,ii) = squeeze(Pre_G(:,:,mm));
    end
    ddL_theta_m_second(:,:,ii) = G(:,:,ii);

    % (2) Phi
    B = H*Gamma'*pffdata(:,:,ii)*Gamma*H';
    ddL_Phi(:,:,ii) = 2*iSe*kron( real(B),eye(Ndof) );

    % (3) Gamma
    C = H'*(Phi'*Phi)*H;
    ddL_Gamma(:,:,ii) = 2*iSe*real( kron(C.',pffdata(:,:,ii)) );

    % (4) Se
    A = (tdata_fft_row(:,ii) - Phi*H*Gamma'*fdata_fft_row(:,ii));
    ddL_Se(:,:,ii) = -Ndof*1*iSe2+2*iSe3*(A'*A) ;
    % ********************* Pivot terms End *******************************

    % ------------------------- cross terms -------------------------
    % (1) theta_m and Phi
    pre_ddL_Phi_theta_m = -Se^(-1)*kron(Gamma'*ptfdata(:,:,ii)',eye(Nmode))*PPhi + iSe*kron(Gamma'*pffdata(:,:,ii)*Gamma*conj(H),Phi') + ...
        iSe*kron(Gamma'*pffdata(:,:,ii)*Gamma*conj(H)*Phi',eye(Nmode))*PPhi;
    pre_ddL_Phi_theta_mconj = conj(pre_ddL_Phi_theta_m);
    ddL_Phi_theta_m(:,:,ii) = [vec_dHk2theta_m.'   vec_dconjHk2theta_m.'] * kron(eye(2),Ld') * [pre_ddL_Phi_theta_m; pre_ddL_Phi_theta_mconj];

    % (2) theta_m and Gamma
    pre_ddL_Gamma_theta_m = -iSe*kron(eye(Nmode),Phi'*conj(ptfdata(:,:,ii))) + iSe*kron(eye(Nmode),Phi'*Phi*conj(H)*Gamma'*conj(pffdata(:,:,ii))) + ...
        iSe*kron(Gamma'*pffdata(:,:,ii),Phi'*Phi*conj(H))*PGamma;
    pre_ddL_Gamma_theta_mconj = conj(pre_ddL_Gamma_theta_m);
    ddL_Gamma_theta_m(:,:,ii) = [vec_dHk2theta_m.'   vec_dconjHk2theta_m.'] * kron(eye(2),Ld') * [pre_ddL_Gamma_theta_m; pre_ddL_Gamma_theta_mconj];

    % (3) theta_m and Se
    ddL_Se_vecHk = iSe2*vec( Phi'*(conj(ptfdata(:,:,ii)) - Phi*conj(H)*Gamma'*conj(pffdata(:,:,ii)) )*Gamma );
    ddL_Se_vecHkconj = iSe2*vec( Phi'*(tdata_fft_row(:,ii)*fdata_fft_row(:,ii)' - Phi*H*Gamma'*pffdata(:,:,ii) )* Gamma );
    ddL_Se_theta_m(:,:,ii) = [vec_dHk2theta_m.'  vec_dconjHk2theta_m.' ] * kron(eye(2),Ld') * [ddL_Se_vecHk ; ddL_Se_vecHkconj];

    % (4) phi and Gamma
    ddL_Phi_Gamma_first = -2*iSe*real(kron(conj(H),ptfdata(:,:,ii)));
    ddL_Phi_Gamma_second = 2*iSe*real(kron(conj(H),Phi*H*Gamma'*pffdata(:,:,ii)) + kron(conj(H)*Gamma'*conj(pffdata(:,:,ii)),Phi*H) * PGamma );
    ddL_Phi_Gamma(:,:,ii) = ddL_Phi_Gamma_first + ddL_Phi_Gamma_second;

    % (5) phi and Se
    ddL_Phi_Se(:,:,ii) = 2*iSe2*vec(real(ptfdata(:,:,ii)*Gamma*H')) - 2*iSe2*kron(real(B),eye(Ndof)) *vec(Phi);

    % (6) Gamma and Se
    ddL_Se_Gamma(:,:,ii) = 2*iSe2*vec(real(ptfdata(:,:,ii)'*Phi*H)) - 2*iSe2*real( kron(C.',pffdata(:,:,ii)) ) *vec(Gamma);

end
% *********** begin to assemble pivot terms into Hessian matrix ***********

% set parameter indices first
Efz = 1:2*Nmode;
EPhi = 2*Nmode+1:2*Nmode+Nmode*Ndof;
EGamma = 2*Nmode+Nmode*Ndof+1:2*Nmode+Nmode*Ndof+Nmode*nd;
ESe = 2*Nmode+Nmode*Ndof+Nmode*nd+1;

% pivot terms
Hess.ddL_theta_m = sum(ddL_theta_m_first,3) + sum(ddL_theta_m_second,3);
Hess.ddL_Phi = sum(ddL_Phi,3);
Hess.ddL_Gamma = sum(ddL_Gamma,3);
Hess.ddL_Se = sum(ddL_Se,3);

% cross terms
Hess.ddL_Phi_theta_m = sum(ddL_Phi_theta_m,3);
Hess.ddL_Gamma_theta_m = sum(ddL_Gamma_theta_m,3);
Hess.ddL_Se_theta_m = sum(ddL_Se_theta_m,3);
Hess.ddL_Phi_Gamma = sum(ddL_Phi_Gamma,3);
Hess.ddL_Phi_Se = sum(ddL_Phi_Se,3);
Hess.ddL_Se_Gamma = sum(ddL_Se_Gamma,3);


Hessian = BayEma_SE_HessAssm(Hess,Nmode,Ndof,nd);
