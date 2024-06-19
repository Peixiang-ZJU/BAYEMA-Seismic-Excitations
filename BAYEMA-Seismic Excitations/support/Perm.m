function P = Perm(A)
% Function to calculate the permutation matrix for 
%       vec(Gamma') = K*vec(Gamma)  Gamma is d-by-m 
%       vec(Phi') = K*vec(Phi)      Phi is n-by-m   
% n is the no. of DoFs, m is the no. of modes in the selected band, d is
% the no. of excitation direction


[N,M] = size(A); 
Num = 1:N*M;
NNum = vec(reshape(Num,M,N)');

[~,Locb] = ismember(NNum,Num );
P = zeros(N*M,N*M);
sz = [N*M N*M];

ind = sub2ind(sz,Num',Locb);
P(ind) = 1;
P = P';
