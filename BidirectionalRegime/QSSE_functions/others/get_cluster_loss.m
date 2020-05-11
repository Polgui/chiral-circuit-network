function [Gamma_target,Lambda_target]=get_cluster_loss(N_sys)

dimS=2;

dimB=3;

Gamma_target=cell(1,N_sys+1+N_sys);
Lambda_target=cell(1,N_sys+1+N_sys);

for iters=1:N_sys
    Gamma_target{2*(iters-1)+1}=reshape([1/sqrt(2),-1/sqrt(2),zeros(1,dimS-2)],[1,1,dimS]);
    Lambda_target{2*(iters-1)+1}=1;
    Gamma_target{2*(iters-1)+2}=reshape([1,0,0],[1,1,dimB]);
    Lambda_target{2*(iters-1)+2}=1;
end

idS=eye(dimS);

P1=zeros(dimS,dimS);
P2=zeros(dimS,dimS);
P1(1,1)=1;
P2(2,2)=1;

CZ=reshape(KroneckerProduct(idS,idS)-2*KroneckerProduct(P1,P1),[dimS,dimS,dimS,dimS]);

for iters=1:N_sys-1
    [Gamma_target,Lambda_target]=Shift_site1_to_siteN(Gamma_target,Lambda_target,2*(iters)+1,2*(iters-1)+2,256,1e-8);
    [Gamma_target,Lambda_target]=LocalUnitaryEvolution(Gamma_target,Lambda_target,2*(iters-1)+1,2,CZ,256,1e-8);
    [Gamma_target,Lambda_target]=Shift_site1_to_siteN(Gamma_target,Lambda_target,2*(iters-1)+2,2*(iters)+1,256,1e-8);
end

Lambda_target{2*N_sys+1}=1;
Gamma_target{2*N_sys+1}=reshape([0,1,0],[1,1,dimB]);

end

