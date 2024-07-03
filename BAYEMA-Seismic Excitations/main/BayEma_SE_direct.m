function out = BayEma_SE_direct(in)
%%
%             Main function for calculating the MPV and PCM
%     ================================================================
%             Peixiang WANG |Zhejiang University| June 2024
%     ================================================================
%%
tic; % start counting computational time

%% extract fields within single band to be identified
f0 = cell2mat(in.f0); fs = in.fs; f1f2 = in.f1f2;
tdata = in.tdata; fdata = in.fdata;

maxiter = in.maxiter;
tol_cvg = in.tol_cvg;


Nmode = length(f0); % no. of modes in a band
nd = size(fdata,2); % No. of input excitation direction
dt = 1/fs; % Time step
[nt,Ndof] = size(tdata);
if nt < Ndof
    error('tdata should be set as (Nt,Ndof)')
end

%% prepare data
% Calcualte the scaled fft value of tdata
[tdata_fft,~,~] = fft_scale(tdata,dt);
tdata_fft = tdata_fft(2:end,:);

% Calculate the scaled fft value of fdata
[fdata_fft,~,freq] = fft_scale(fdata,dt);
% Ignore the 1st fft data
freq = freq(2:end);
fdata_fft = fdata_fft(2:end,:);

% set frequency indices based on the lower and upper frequency bound
II = find(freq >= f1f2(1) & freq <= f1f2(2));
% Extract tdata within the selected band
tdata_fft_band = tdata_fft(II,:);
% Extract fdata within the selected band
fdata_fft_band = fdata_fft(II,:);

freq_band = freq(II);
Nf = length(II);

%% set initial values of parameter to start the iteration
x0 = repmat([1;1],Nmode,1);
z0 = 0.01*ones(1,Nmode);
% initial guess of mode shape
phi0 = cal_phi0(fdata,tdata,fs,f0);
phi0 = normc(phi0);

%% Do minimazation to find MPV based on coordinate descend algorithm
for iter = 1:maxiter

    tdata_fft_row = transpose(tdata_fft_band); % (Number of DoFs,Nf)
    fdata_fft_row = transpose(fdata_fft_band); % (Number of DoFs,Nf)

    beta = f0./freq_band;
    h = ( 1-(beta.^2) - 1i*2.*z0.*beta ).^(-1); % FRF hk

    % --------------------- calculate Gamma --------------------
    Pre_Gamma_Left = zeros(nd*Nmode,nd*Nmode);
    Pre_Gamma_middle = zeros(nd*Nmode,Ndof*Nmode);

    for ii = 1:Nf
        H = diag(h(ii,:));  % Hk = diag(hk)
        C = H'*(phi0'*phi0)*H;
        Pre_Gamma_Left = Pre_Gamma_Left + kron(C.',fdata_fft_row(:,ii)*fdata_fft_row(:,ii)');
        Pre_Gamma_middle = Pre_Gamma_middle + kron(H,fdata_fft_row(:,ii)*tdata_fft_row(:,ii)');
    end
    Gamma = real(Pre_Gamma_Left) \ ( real( Pre_Gamma_middle ) * vec(phi0) );
    Gamma = reshape(Gamma,[],Nmode);
    % ***********************************************************

    % --------------------- calculate Phi -----------------------
    B = zeros(Nmode,Nmode);
    Pre_Phi_Left = zeros(Ndof,Nmode);
    for ii = 1:Nf
        H = diag(h(ii,:));  % Hk = diag(hk)
        B = B + H*Gamma'*fdata_fft_row(:,ii)*fdata_fft_row(:,ii)'*Gamma*H';
        Pre_Phi_Left = Pre_Phi_Left + tdata_fft_row(:,ii)*fdata_fft_row(:,ii)'*Gamma*H';
    end
    Phi = real( Pre_Phi_Left ) / real( B );
    % ***********************************************************

    % ---- re-normlize Phi and Gamma to ensure Phi'*Phi = 1 -----
    Phi_norm = Phi/diag(sqrt(diag(Phi'*Phi)));
    Gamma_norm = Gamma*diag(sqrt(diag(Phi'*Phi)));
    % ***********************************************************

    % --------------------- calculate Se ------------------------
    for ii = 1:Nf
        H = diag(h(ii,:));  % Hk = diag(hk)
        Pre_A = tdata_fft_row(:,ii) - Phi_norm*H*Gamma_norm'*fdata_fft_row(:,ii);
        A(ii) = Pre_A'*Pre_A;
    end
    Se = sum(A)/Ndof/Nf;
    % ***********************************************************

    % calculate the likelihood function value
    L(iter) = sum(A)/Se+Ndof*Nf*log(pi)+Ndof*Nf*log(Se);

    % --------------- numerically optimize f and z ---------------
    options = optimset('TolFun',1e-6,'TolX',1e-6,'MaxFunEvals',10000,'MaxIter',10000);
    [x,~,~] = fminsearch('BayEma_SE_obj',x0,options,tdata_fft_band,fdata_fft_band,freq_band,f0,z0,Phi_norm,Gamma_norm);

    f0 = f0.*x(1:2:end)';
    z0 = z0.*x(2:2:end)';
    phi0 = Phi_norm;

    fprintf('Iter: %1.0f   ', iter)
    fprintf('freq(Hz): %4.3f   ', f0)
    fprintf('damping: %4.2e   ', z0)
    fprintf('Se: %4.2e   \n', Se(:))


    % convergence criteria
    if iter > 1 && abs((L(iter)-L(iter-1))/L(iter)) < tol_cvg
        break
    end

end

% 'MPV' contains the MPV of modal parameters
MPV.f = f0;
MPV.z = z0;
MPV.Phi = Phi_norm;
MPV.Gamma = Gamma_norm;
MPV.Se = Se;
out = MPV;

%% Uncertainty quantification (analytically gives c.o.v.)

% calculate the Hessian matrix using the obtained MPV
[raw_Hessian] = BayEma_SE_hessian(MPV,tdata_fft_band,fdata_fft_band,freq_band);

% eliminate the effect of different orders of magnitude
fz = vecT([f0(:),z0(:)]');
MPVs = blkdiag(diag(fz), eye(Ndof*Nmode), eye(nd*Nmode), Se );
Hessian = MPVs*raw_Hessian*MPVs;

% when calculating PCM, take mode shape norm constraint into account
g2theta = zeros(Nmode,size(Hessian,1));
for mm = 1:Nmode
    g2theta(mm,2*Nmode+1+(mm-1)*Ndof:2*Nmode+mm*Ndof) = 2*Phi_norm(:,mm)';
end
U_theta = null(g2theta);
CoeMat = U_theta*((U_theta'*Hessian*U_theta)\U_theta');


coefv = cal_coefv(CoeMat,Nmode,Ndof,nd,Gamma_norm);
out.coefv = coefv;

toc; % record computation time


%% make response prediction
if in.myoptions(1) == 1

    nsample = in.nsample;

    MPVs = blkdiag(eye(2*Nmode,2*Nmode), eye(Ndof*Nmode), eye(nd*Nmode), Se );
    Hessian = MPVs*raw_Hessian*MPVs;

    VarMat = U_theta*((U_theta'*Hessian*U_theta)\U_theta');

    mpv = [vec(fz);vec(Phi_norm);vec(Gamma_norm)];
    VarMat(1:end-1,1:end-1) = 0.5*(VarMat(1:end-1,1:end-1) + VarMat(1:end-1,1:end-1)');
    pred_mp = mvnrnd(mpv, VarMat(1:end-1,1:end-1),nsample);

    pred.f = pred_mp(:,1:2:2*Nmode);
    pred.z = pred_mp(:,2:2:2*Nmode);
    pred.phi = pred_mp(:,2*Nmode+1:2*Nmode+Nmode*Ndof);
    pred.Gamma = pred_mp(:,2*Nmode+Nmode*Ndof+1:2*Nmode+Nmode*Ndof+Nmode*nd);


    out.pred = pred;


end




end



