function A = optfield(B,C)
% Function to assign the union of structures B and C to structure A.

A = C;
fB = fieldnames(B);
nB = length(vec(fB));

for ii=1:nB
  A.(fB{ii}) = B.(fB{ii});
end

