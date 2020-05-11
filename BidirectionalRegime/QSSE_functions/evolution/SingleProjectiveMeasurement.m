function [Gamma,Lambda]=SingleProjectiveMeasurement(Gamma,Lambda,index,maxSchmidtrank,Eps)

rho=single_site_reduced_state_Gamma_efficient(Gamma,Lambda,index);
projE=rho(2,2)/trace(rho);
r=rand;
if r>projE
    Gamma{index}=tensor_contraction(Gamma{index},[1,0;0,0],3,2);
elseif r<=projE
    Gamma{index}=tensor_contraction(Gamma{index},[0,1;0,0],3,2);
    Gamma{index-1}=tensor_contraction(Gamma{index-1},[1,0;0,-1],3,2);
    if size(Gamma{index+1},3)==2
        Gamma{index+1}=tensor_contraction(Gamma{index+1},[1,0;0,-1],3,2);
    end
else
    error('wrong r')
end

[Gamma,Lambda]= sweep(Gamma,Lambda,1,length(Gamma),maxSchmidtrank, Eps);

end
