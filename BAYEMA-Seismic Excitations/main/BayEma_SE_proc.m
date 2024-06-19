function out = BayEma_SE_proc(in,fun)

% check mandatory fields from IN
if ~isfield(in,'tdata')
  error('tdata is required field.');
end
if ~isfield(in,'fdata')
  error('fdata is required field.');
end
if ~isfield(in,'fs')
  error('fs is required field');
end
if ~isfield(in,'f0')
  error('f0 is required field.');
end
if ~isfield(in,'f1f2')
  error('f1f2 is required field.');
end

% run corresponding fun for ID
if isempty(in.iband)
    in.iband = 1:length(in.f0);
else
    in.f0 = in.f0(in.iband);
    in.f1f2 = in.f1f2(in.iband,:);
end

nb = size(in.f1f2,1); % no. of freq. bands
o = cell(1,nb);
in0 = in;
for ii = 1:nb
  in = in0;
  in.f0 = in0.f0(ii);
  in.f1f2 = in0.f1f2(ii,:);
  o{ii} = feval(fun,in);
end


out = o{1};
varx = fieldnames(o{1}).'; % cell array of field names
for ii=2:nb
  for xx = varx
    if isstruct(out.(xx{:}))
      vary = fieldnames(out.(xx{:})).';
      for yy=vary
        out.(xx{:}).(yy{:}) = [out.(xx{:}).(yy{:}),o{ii}.(xx{:}).(yy{:})];
      end
    else
      out.(xx{:}) = [out.(xx{:}),o{ii}.(xx{:})];
    end
  end
end





