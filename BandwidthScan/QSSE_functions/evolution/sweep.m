function [Gamma,Lambda]=sweep(Gamma,Lambda,start,finish,maxSchmidtrank, Eps)

for iter=start:finish-1
    size1=size(Gamma{iter},3);
    size2=size(Gamma{iter+1},3);
    id=reshape(eye(size1*size2),[size1,size2,size1,size2]);
    [Gamma,Lambda]= LocalUnitaryEvolution(Gamma,Lambda,iter,2,id,maxSchmidtrank, Eps);
end
for iter=finish:-1:start+1
    size1=size(Gamma{iter-1},3);
    size2=size(Gamma{iter},3);
    id=reshape(eye(size1*size2),[size1,size2,size1,size2]);
    [Gamma,Lambda]= LocalUnitaryEvolution(Gamma,Lambda,iter-1,2,id,maxSchmidtrank, Eps);
end

end

