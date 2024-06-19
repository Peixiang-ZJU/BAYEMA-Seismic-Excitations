function [phi0] = cal_phi0(fdata,tdata,fs,f0)
% Function to get the initial values of mode shape from FRF of data

nd = size(fdata,2);

if nd == 1

     [FRF_values,f] = tfestimate(fdata,tdata,[],[],[],fs);
        [~,f0_index] = min( abs(f - f0) );

        for mm = 1:length(f0)
            [U,~,~] = svd(squeeze(  imag( FRF_values(f0_index(mm),:,:) )'  ));
            phi0(:,mm) = U(:,1); % initial guess of mode shapes
        end

else 
        [FRF_values,f] = tfestimate(fdata,tdata,[],[],[],fs,'mimo');
        [~,f0_index] = min( abs(f - f0) );

        for mm = 1:length(f0)
            [U,~,~] = svd(squeeze(  imag( FRF_values(f0_index(mm),:,:) )  ));
            phi0(:,mm) = U(:,1); % initial guess of mode shapes
        end
       

end



