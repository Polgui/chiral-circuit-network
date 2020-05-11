function Fid=state_overlap(Gamma,Lambda,Gamma_target,Lambda_target)
%% WARNING: WORKS ONLY FOR PURE STATES

if size(Gamma)~=size(Gamma_target)
    error 'Mismatching sizes of Gamma'
elseif size(Lambda)~=size(Lambda_target)
    error 'Mismatching sizes of Lambda' 
end

Nqubits=length(Gamma)-1;

tens=Lambda_multiplication(Gamma{1},Lambda{1},2);

tens=tensor_contraction(tens,conj(Lambda_multiplication(Gamma_target{1},Lambda_target{1},2)),[1 3],[1 3]);

for iters=2:Nqubits+1
   
tens=tensor_contraction(tens,Lambda_multiplication(Gamma{iters},Lambda{iters},2),1,1);
tens=tensor_contraction(tens,conj(Lambda_multiplication(Gamma_target{iters},Lambda_target{iters},2)),[1 3],[1 3]);
    
end

Fid=tens;

end

