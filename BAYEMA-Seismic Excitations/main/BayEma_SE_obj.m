function [A] = BayEma_SE_obj(x0,tdata_fft,fdata_fft,freq_band,f0,z0,Phi0,Gamma)
% % Function to numerically optimize f and z

f_tmp = x0(1:2:end).*f0';
z_tmp = x0(2:2:end).*z0';
Phi = Phi0;

Nmode = length(f_tmp);
[Nf,~] = size(tdata_fft);

tdata_fft_row = transpose(tdata_fft); % (Number of DoFs,Nf)
fdata_fft_row = transpose(fdata_fft); % (Nd,Nf)

% generate Ld, which is the transfer matrix between vecd(Hk) and vec(Hk)
Ld = zeros(Nmode^2,Nmode);
Ld(1:Nmode+1:end,:) = eye(Nmode);

if Ld'*Ld ~= eye(Nmode)
    error('transpose(Ld)*Ld should be equal to eye(Nmode)')
end

beta = f_tmp'./freq_band; % Row vector (Number of modes)
h = ( 1-(beta.^2) - 1i*(2.*z_tmp'.*beta) ).^-1; % Row vector (Number of modes)

for ii = 1:Nf
    H = diag(h(ii,:)); 
    Pre_A(:,ii) = tdata_fft_row(:,ii) - Phi*H*Gamma'*fdata_fft_row(:,ii);
    A(ii) = Pre_A(:,ii)'*Pre_A(:,ii);
end
A = sum(A);

end

