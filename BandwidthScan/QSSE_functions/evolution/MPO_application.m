function [Gamma,Lambda]=MPO_application(Gamma,Lambda,MPO,site1,N)

if length(MPO)~=N
    error 'length(MPO)~=N'
end

tens1=tensor_contraction(Gamma{site1},MPO{1},3,1);

if site1>1
    tens1=Lambda_multiplication(tens1,Lambda{site1-1},1);
end
tens1=Lambda_multiplication(tens1,Lambda{site1},2);

for itern=2:N
    tens2=tensor_contraction(Lambda_multiplication(Gamma{site1+itern-1},Lambda{site1+itern-1},2),MPO{itern},3,1);
    tens=tensor_contraction(tens1,tens2,[2,4],[1,4]);
    dimension=size(tens);
    tensmerge=reshape(tens,[dimension(1)*dimension(2),dimension(3)*dimension(4)*dimension(5)]);
    [U,L1,V,~]=MySVD(tensmerge,1000,1e-6);
    Lambda{site1+itern-2}=L1;
    Gamma{site1+itern-2}=permute(reshape(U,[dimension(1),dimension(2),length(L1)]),[1,3,2]);
    if site1+itern-1>2
        Gamma{site1+itern-2}=Lambda_multiplication(Gamma{site1+itern-2},1./Lambda{site1+itern-3},1);
    end
    tens1=reshape(diag(L1)*(V'),[length(L1),dimension(3),dimension(4),dimension(5)]);
end

Gamma{site1+N-1}=Lambda_multiplication(tens1,1./Lambda{site1+N-1},2);
if N>1
    Gamma{site1+N-1}=Lambda_multiplication(Gamma{site1+N-1},1./Lambda{site1+N-2},1);
end

end

