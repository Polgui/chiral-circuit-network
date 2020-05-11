function [Gamma_target,Lambda_target]=get_cluster(N_sys)

dimS=2;

Gamma_target=cell(1,N_sys+1);
Lambda_target=cell(1,N_sys+1);

for iters=1:N_sys
    Gamma_target{iters}=reshape([1/sqrt(2),-1/sqrt(2),zeros(1,dimS-2)],[1,1,dimS]);
    Lambda_target{iters}=1;
end

idS=eye(dimS);

P1=zeros(dimS,dimS);
P2=zeros(dimS,dimS);
P1(1,1)=1;
P2(2,2)=1;

CZ=reshape(KroneckerProduct(idS,idS)-2*KroneckerProduct(P1,P1),[dimS,dimS,dimS,dimS]);

for iters=1:N_sys-1
    [Gamma_target,Lambda_target]=LocalUnitaryEvolution(Gamma_target,Lambda_target,iters,2,CZ,256,1e-8);
end

Lambda_target{N_sys+1}=1;
Gamma_target{N_sys+1}=reshape([1,0],[1,1,2]);

end

