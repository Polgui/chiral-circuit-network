function Fid=state_overlap_v2(Gamma,Lambda,Gamma_target,Lambda_target)
%% DO NOT USE, DOES NOT WORK

if size(Gamma)~=size(Gamma_target)
    error 'Mismatching sizes of Gamma'
elseif size(Lambda)~=size(Lambda_target)
    error 'Mismatching sizes of Lambda' 
end

Nqubits=length(Gamma)-1;

tens=Lambda_multiplication(Gamma{1},Lambda{1},2);

tens=tensor_contraction(tens,conj(Lambda_multiplication(Gamma_target{1},Lambda_target{1},2)),3,3);

tens=permute(tens,[1,3,2,4]);

for iters=2:Nqubits+1
   
tens=tensor_contraction(tens,Lambda_multiplication(Gamma{iters},Lambda{iters},2),length(size(tens))-1,1);
tens=tensor_contraction(tens,conj(Lambda_multiplication(Gamma_target{iters},Lambda_target{iters},2)),[length(size(tens))-2 length(size(tens))],[1 3]);
    
end

Fid=tens;

end

