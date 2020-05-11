function [Gamma_target,Lambda_target]=get_2Dcluster_evenlinks(Nqubits,Ncluster,iterqubit)
%% Cluster state in another basis

ground=[1,0]';
sx=[0,1;1,0];

dimS=2;

Gamma_target=cell(1,Ncluster*(Nqubits+1));
Lambda_target=cell(1,Ncluster*(Nqubits+1));

for iters=1:Ncluster*(Nqubits+1)
    Gamma_target{iters}=reshape(ground,[1,1,dimS]);
    Lambda_target{iters}=1;
end

Uhad=-expm(-1j*pi*sx/4);
for itercluster=1:Ncluster
    for iterq=1:Nqubits
        Gamma_target{(Nqubits+1)*(itercluster-1)+iterq}=tensor_contraction(Gamma_target{(Nqubits+1)*(itercluster-1)+iterq},Uhad,3,2);
    end
end

idS=eye(dimS);

P2=zeros(dimS,dimS);
P2(2,2)=1;

CZ=reshape(KroneckerProduct(idS,idS)-2*KroneckerProduct(P2,P2),[dimS,dimS,dimS,dimS]);

for itercluster=1:Ncluster
    for iterq=1:Nqubits-1
        [Gamma_target,Lambda_target]=LocalUnitaryEvolution(Gamma_target,Lambda_target,(Nqubits+1)*(itercluster-1)+iterq,2,CZ,256,1e-8);
    end
end

CZ=reshape(KroneckerProduct(idS,idS)-2*KroneckerProduct(P2,P2),[dimS,dimS,dimS,dimS]);

for iterq=1:iterqubit
    for iterclusterpair=1:floor(Ncluster/2)
        [Gamma_target,Lambda_target]=Shift_site1_to_siteN(Gamma_target,Lambda_target,(iterclusterpair-1)*2*(Nqubits+1)+iterq,(iterclusterpair-1)*2*(Nqubits+1)+iterq+Nqubits,256, 1e-8);
        [Gamma_target,Lambda_target]=LocalUnitaryEvolution(Gamma_target,Lambda_target,(iterclusterpair-1)*2*(Nqubits+1)+iterq+Nqubits,2,CZ,256,1e-8);
        [Gamma_target,Lambda_target]=Shift_site1_to_siteN(Gamma_target,Lambda_target,(iterclusterpair-1)*2*(Nqubits+1)+iterq+Nqubits,(iterclusterpair-1)*2*(Nqubits+1)+iterq,256, 1e-8);
    end
end

end

