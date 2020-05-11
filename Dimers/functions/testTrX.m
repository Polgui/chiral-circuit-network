psi1=[1;1]/sqrt(2);
psi2=[1;0];
psi3=[0;1];


rho=KroneckerProduct(psi1*psi1',psi2*psi2',psi3*psi3');

rhopost=TrX(TrX(rho,3,[2,2,2]),1,[2,2])


rhopost=TrX(rho,[3,1],[2,2,2])