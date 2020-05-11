function Fid=state_fidelity(Gamma,Lambda,Lambda_left,Gamma_target,Lambda_target,Lambda_left_target)

if size(Gamma)~=size(Gamma_target)
    error 'Mismatching sizes of Gamma'
elseif size(Lambda)~=size(Lambda_target)
    error 'Mismatching sizes of Lambda' 
end

Ntot=length(Gamma);

if length(Gamma)<3
    error 'length(Gamma)<3'
end

tens=Lambda_multiplication(Gamma{1},Lambda_left,1);
tens=Lambda_multiplication(tens,Lambda{1},2);
tens=tensor_contraction(tens,conj(tens),1 ,1);

temp=Lambda_multiplication(Gamma{2},Lambda{2},2);
tens=tensor_contraction(tens,temp,1,1);
tens=tensor_contraction(tens,conj(temp),2,1);

temp=Lambda_multiplication(Gamma_target{1},Lambda_left_target,1);
temp=Lambda_multiplication(temp,Lambda_target{1},2);
tens=tensor_contraction(tens,conj(temp),1,3);
tens=tensor_contraction(tens,temp,[1,6],[3,1]);

for iters=3:Ntot-1
    temp=Lambda_multiplication(Gamma{iters},Lambda{iters},2);
    tens=tensor_contraction(tens,temp,1,1);
    tens=tensor_contraction(tens,conj(temp),2,1);
    
    temp=Lambda_multiplication(Gamma_target{iters-1},Lambda_target{iters-1},2);
    tens=tensor_contraction(tens,conj(temp),[1,3],[3,1]);
    tens=tensor_contraction(tens,temp,[1,2],[3,1]);
    
end

temp=Lambda_multiplication(Gamma{Ntot},Lambda{Ntot},2);
tens=tensor_contraction(tens,temp,1,1);
tens=tensor_contraction(tens,conj(temp),[2,6],[1,2]);

temp=Lambda_multiplication(Gamma_target{Ntot-1},Lambda_target{Ntot-1},2);
tens=tensor_contraction(tens,conj(temp),[1,3],[3,1]);
tens=tensor_contraction(tens,temp,[1,2],[3,1]);
    
temp=Lambda_multiplication(Gamma_target{Ntot},Lambda_target{Ntot},2);
temp=tensor_contraction(conj(temp),temp,2,2);
tens=tensor_contraction(tens,temp,[1,2,3,4],[2,4,1,3]);

Fid=tens;

end

