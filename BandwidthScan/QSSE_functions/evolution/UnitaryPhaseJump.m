function [Gamma,Lambda,njump]=UnitaryPhaseJump(Gamma,Lambda,index,U_phase,U_jump,coop,njump,maxSchmidtrank, Eps)

tens=Lambda_multiplication(Gamma{index},Lambda{index},2);
if index>1
    tens=Lambda_multiplication(tens,Lambda{index-1},1);
end
tens=tensor_contraction(tens,Lambda_multiplication(Gamma{index+1},Lambda{index+1},2),2,1);
rho=tensor_contraction(tens,conj(tens),[1,3],[1,3]);
Probjump=rho(3,2,3,2)*4*coop/(1+coop)^2;
r=rand;
if r<=Probjump && njump==0
    [Gamma,Lambda] = LocalUnitaryEvolution(Gamma,Lambda,index,2,U_jump,maxSchmidtrank, Eps);
    [Gamma,Lambda]= sweep(Gamma,Lambda,1,length(Gamma),maxSchmidtrank, Eps);
    njump=njump+1;
elseif r>Probjump && njump==0
    [Gamma,Lambda] = LocalUnitaryEvolution(Gamma,Lambda,index,2,U_phase,maxSchmidtrank, Eps);
end

end

