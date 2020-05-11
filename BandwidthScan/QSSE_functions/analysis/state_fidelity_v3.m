function Fid=state_fidelity_v3(Gamma,Lambda,Lambda_left,Gamma_target,Lambda_target,Lambda_left_target)
if size(Gamma)~=size(Gamma_target)
    error 'Mismatching sizes of Gamma'
elseif size(Lambda)~=size(Lambda_target)
    error 'Mismatching sizes of Lambda' 
end

Ntot=length(Gamma);

if length(Gamma)<2
    tens=Lambda_multiplication(Gamma{1},Lambda_left,1);
    tens=Lambda_multiplication(tens,Lambda{1},2);

    temp=Lambda_multiplication(Gamma_target{1},Lambda_left_target,1);
    temp=Lambda_multiplication(temp,Lambda_target{1},2);
    
    Fid=abs(tensor_contraction(tens,conj(temp),[1,2,3],[1,2,3]))^2;
    
else

    tens=Lambda_multiplication(Gamma{1},Lambda_left,1);
    tens=Lambda_multiplication(tens,Lambda{1},2);

    temp=Lambda_multiplication(Gamma_target{1},Lambda_left_target,1);
    temp=Lambda_multiplication(temp,Lambda_target{1},2);
    temp2=Lambda_multiplication(Gamma_target{2},Lambda_target{2},2);
    temp=tensor_contraction(temp,temp2,2,1);

    tens=tensor_contraction(tens,conj(temp),3,2);

    for iters=2:Ntot-1

        if mod(iters,2)==0

            temp=Lambda_multiplication(Gamma{iters},Lambda{iters},2);
            temp2=Lambda_multiplication(Gamma{iters+1},Lambda{iters+1},2);
            temp=tensor_contraction(temp,temp2,2,1);
            tens=tensor_contraction(tens,temp,[2,5],[1,2]);

        elseif mod(iters,2)==1

            temp=Lambda_multiplication(Gamma_target{iters},Lambda_target{iters},2);
            temp2=Lambda_multiplication(Gamma_target{iters+1},Lambda_target{iters+1},2);
            temp=tensor_contraction(temp,temp2,2,1);
            tens=tensor_contraction(tens,conj(temp),[3,5],[1,2]);

            tens=permute(tens,[1,3,2,4,5]);

        else
            error 'mod(iters,2) is neither 0 nor 1'
        end

    end

    if mod(Ntot,2)==0
        temp=Lambda_multiplication(Gamma{Ntot},Lambda{Ntot},2);
        tens=tensor_contraction(tens,temp,[2,5],[1,3]);
    elseif mod(Ntot,2)==1
        temp=Lambda_multiplication(Gamma_target{Ntot},Lambda_target{Ntot},2);
        tens=tensor_contraction(tens,conj(temp),[3,5],[1,3]);
    else 
        error 'mod(Ntot,2) is neither 0 nor 1'
    end

    Fid=tensor_contraction(tens,conj(tens),[1:length(size(tens))],[1:length(size(tens))]);

    if imag(Fid) > 1e-2
        error 'imag(Fid)>1e-2'
    end
    
end

end

