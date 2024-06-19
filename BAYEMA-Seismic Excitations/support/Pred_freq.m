function [pred_tdata] = Pred_freq(in)


Phi = in.phi';
Gamma = in.Gamma;
f = in.f;
z = in.z;

fdata = in.fdata;
fs = in.fs;

% eplison_k = in.eplison_k;


dt = 1/fs;
[Nt,~] = size(fdata);



[scale_fft,w,~] = fft_scale(fdata,dt);
fdata_fft = scale_fft;
freq = w(1:end)/2/pi;

% tdata = in.tdata;
% [scale_tdata_fft,w,~] = fft_scale(tdata,dt);
% scale_tdata_fft = scale_tdata_fft(1:end,:);
% freq = w(1:end)/2/pi;

beta = f./freq; % Row vector (Number of modes)
Hk = ( 1-(beta.^2) - 1i*2.*z.*beta ).^(-1); % Row vector (Number of modes)

for ii = 1:length(freq)
    Yk(ii,:) = Phi * diag(Hk(ii,:)) * Gamma' * fdata_fft(ii,:).' ;
     % Yk(ii,:) = Phi * diag(Hk(ii,:)) * Gamma * fdata_fft(ii,:).' + eplison_k;
end


tdata_fft = [Yk(1:end,:);flip(conj(Yk(2:end-1,:)))]./sqrt(dt/Nt);


% plot(abs(scale_tdata_fft(1:200,end))./sqrt(dt/Nt),'r'); 
% hold on
% plot(abs(tdata_fft(1:200,end)),'k'); 



pred_tdata = real(ifft(tdata_fft));


