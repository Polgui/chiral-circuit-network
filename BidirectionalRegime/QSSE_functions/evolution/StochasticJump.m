function [Gamma,Lambda,njump]=StochasticJump(Gamma,Lambda,index,jumpG,jumpE,P_loss,njump,maxSchmidtrank,Eps)

rho=single_site_reduced_state_Gamma_efficient(Gamma,Lambda,index);
probJumpG=trace(jumpG*rho*jumpG')/trace(rho);
probJumpE=trace(jumpE*rho*jumpE')/trace(rho);
r=rand;
if r<=1-P_loss && njump==0
    r2=rand;
    if r2>=probJumpE/(probJumpG+probJumpE)
        Gamma{index}=tensor_contraction(Gamma{index},jumpG,3,2);
    elseif r2<probJumpE/(probJumpG+probJumpE)
        Gamma{index}=tensor_contraction(Gamma{index},jumpE,3,2);
    else 
        error 'invalid rand'
    end
    njump=njump+1;
    [Gamma,Lambda]= sweep(Gamma,Lambda,1,length(Gamma),maxSchmidtrank, Eps);
end

end

