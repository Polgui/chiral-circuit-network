function [Gamma_target,Lambda_target]=get_1Dcluster_new(N_sys)

dimS=2;

Gamma_target=cell(1,N_sys);
Lambda_target=cell(1,N_sys);

for iters=1:N_sys
    Gamma_target{iters}=reshape([1/sqrt(2),1/sqrt(2),zeros(1,dimS-2)],[1,1,dimS]);
    Lambda_target{iters}=1;
end

idS=eye(dimS);

P2=zeros(dimS,dimS);
P2(2,2)=1;

CZ=reshape(KroneckerProduct(idS,idS)-2*KroneckerProduct(P2,P2),[dimS,dimS,dimS,dimS]);

for iters=1:N_sys-1
    [Gamma_target,Lambda_target]=LocalUnitaryEvolution(Gamma_target,Lambda_target,iters,2,CZ,256,1e-8);
end


end

