function [y,H] = fun_h(f,z,fk,x)
% y = fft02_funh(f,z,ff,x)
% returns up to the second derivatives of the complex-valued function
% h = 1./((b^2-1)-i*(2*b*z)), b = f/fk(.)
% f = scalar, z = scalar, ff = (nf,1) or (1,nf),
% x = 'f','z','ff','zz','fz','zf'
% [y,h] = ... also returns h

if size(f(:))~=1
    error('f must be scalar');
end
if size(z(:))~=1
    error('z must be scalar');
end

y = size(fk);
beta = f./fk;
h = 1./(1-beta.^2-1i*(2*z.*beta));
H = diag(h);

switch x
    case ''
        y = H;
    case 'f'
        y = H.^2.*(2*beta+2*z*1i);
        y = y./fk;
    case 'z'
        y = H.^2.*(2*beta*1i);
    case 'ff'
        y = 2*H.^3.*(3*beta.^2+1-4*z.^2+6*z.*beta*1i);
        y = y./fk.^2;
    case 'zz'
        y = -8*H.^3.*beta.^2;
    case {'fz','zf'}
        y = 2*1i*H.^3.*(3*beta.^2+1+2.*z.*beta*1i);
        y = y./fk;
    otherwise
        error('Unknown x');
end