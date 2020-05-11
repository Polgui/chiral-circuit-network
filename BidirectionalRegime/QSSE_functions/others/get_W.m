function [Gamma_target,Lambda_target]=get_W(N_sys,dimS,amplitudes)

Gamma_target=cell(1,N_sys);
Lambda_target=cell(1,N_sys);

Gamma_target{1}=zeros(1,2,dimS);
Gamma_target{1}(1,1,1)=1;
Gamma_target{1}(1,2,2)=amplitudes(1)/abs(amplitudes(1));

Lambda_target{1}=zeros(2,1);
Lambda_target{1}(2)=abs(amplitudes(1));
Lambda_target{1}(1)=sqrt(1-abs(amplitudes(1))^2);

Gamma_target{end}=zeros(2,1,dimS);
Gamma_target{end}(2,1,1)=1;
Gamma_target{end}(1,1,2)=amplitudes(N_sys)/abs(amplitudes(N_sys));

Lambda_target{end}=zeros(1,1);
Lambda_target{end}(1)=1;

for iters=2:N_sys-1
    
    Lambda_target{iters}=zeros(2,1);
    Lambda_target{iters}(2)=sqrt((Lambda_target{iters-1}(2))^2+abs(amplitudes(iters))^2);
    Lambda_target{iters}(1)=sqrt(1-abs(Lambda_target{iters}(2))^2);
    
    Gamma_target{iters}=zeros(2,2,dimS);
    Gamma_target{iters}(1,1,1)=1/Lambda_target{iters-1}(1);
    Gamma_target{iters}(2,2,1)=1/Lambda_target{iters}(2);
    Gamma_target{iters}(1,2,2)=amplitudes(iters)/(Lambda_target{iters-1}(1)*Lambda_target{iters}(2));
    
end

end

