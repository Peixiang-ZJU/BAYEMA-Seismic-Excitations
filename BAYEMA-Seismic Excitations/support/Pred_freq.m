function [pred_tdata] = Pred_freq(in)
% function to make structural responses using mode superposition rule in
% the frequency domain

% load sampled modal parameters
Phi = in.phi;
Gamma = in.Gamma;
f = in.f;
z = in.z;
eplison_k = in.eplison_k;

% load new ground motion
fdata = in.fdata;
fs = in.fs;

dt = 1/fs;
[Nt,~] = size(fdata);

[scale_fft,w,~] = fft_scale(fdata,dt);
fdata_fft = scale_fft;
freq = w(1:end)/2/pi;

beta = f./freq; 
Hk = ( 1-(beta.^2) - 1i*2.*z.*beta ).^(-1); 

% calculate responses in the frequncy domain
for ii = 1:length(freq)
     Yk(ii,:) = Phi * diag(Hk(ii,:)) * Gamma' * fdata_fft(ii,:).' + eplison_k;
end
tdata_fft = [Yk(1:end,:);flip(conj(Yk(2:end-1,:)))]./sqrt(dt/Nt);

% get time history using iFFT
pred_tdata = real(ifft(tdata_fft));


