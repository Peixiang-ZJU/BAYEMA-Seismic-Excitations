function [scale_fft,w,freq] = fft_scale(x,dt)
% Function to calculate the scaled FFTs
% Input:  x = data in time domain, (nt,nch)
%             nt = no. of sample points;
%             nch = total no. of measured dofs
%        dt = sampling interval

[Nt,Nx] = size(x);

if Nt < Nx
    errorbox('The data is not arranged in column-by-column');
end

Nw = floor((Nt+2)/2);

dw = 2*pi/Nt/dt;
fftx = sqrt(dt/Nt)*fft(x); % process by column
scale_fft = fftx(1:Nw,:);

w = [0:dw:pi/dt]';
freq = w/2/pi;

