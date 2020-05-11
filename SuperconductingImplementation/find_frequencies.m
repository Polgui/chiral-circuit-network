function [J,Delta1,Delta2,chi1,chi2,chi12]=find_frequencies(EJ1,EJ2,EJ12,EC1,EC2,r12,omegaph,dimB,dphi)

EJ1=EJ1*(1+dphi);
EJ2=EJ2*(1+0*dphi);
EJ12=EJ12*(1+0*dphi);

Cmat=[1/EC1+r12/EC2, -r12/EC2;-r12/EC2, 1/EC2+r12/EC2];
Cmat_inv=pinv(Cmat);

a=diag(sqrt(1:dimB-1),1);
idB=eye(dimB);

N1=1j*((EJ1+EJ12)/(32*Cmat_inv(1,1)))^(1/4)*KroneckerProduct(a'-a,idB);
psi1=((EJ1+EJ12)/(2*Cmat_inv(1,1)))^(-1/4)*KroneckerProduct(a'+a,idB);
N2=1j*((EJ2+EJ12)/(32*Cmat_inv(2,2)))^(1/4)*KroneckerProduct(idB,a'-a);
psi2=((EJ2+EJ12)/(2*Cmat_inv(2,2)))^(-1/4)*KroneckerProduct(idB,a'+a);

temp=expm(1j*(psi1));
cospsi1=(temp+temp')/2;
temp=expm(1j*(psi2));
cospsi2=(temp+temp')/2;
temp=expm(1j*(psi1-psi2));
cospsi12=(temp+temp')/2;

Hs=4*Cmat_inv(1,1)*N1*N1+4*Cmat_inv(2,2)*N2*N2+4*Cmat_inv(1,2)*N1*N2+4*Cmat_inv(2,1)*N2*N1...
    -EJ1*cospsi1-EJ2*cospsi2-EJ12*cospsi12;

Nexc=diag(KroneckerProduct(a'*a,idB)+KroneckerProduct(idB,a'*a));

groundenergy=Hs(1,1);
Hs=Hs-groundenergy*eye(dimB*dimB);

% Second order perturbation on counter-rotating terms

Hs_clear=Hs;
for index1=1:size(Hs_clear,1)
    for index2=1:size(Hs_clear,2)
        if abs(Nexc(index1)-Nexc(index2))>0.5
            Hs_clear(index1,index2)=0;
        end
    end
end
Hs_pert=Hs;
for index1=1:size(Hs_pert,1)
    for index2=1:size(Hs_pert,2)
        if abs(Nexc(index1)-Nexc(index2))<0.5
            Hs_pert(index1,index2)=0;
        end
    end
end

Hs_tot=Hs_clear;

for n1=0:2*(dimB-1)
    for n2=n1+1:2*(dimB-1)
        Hs_tot(abs(Nexc-n1)<1e-1,abs(Nexc-n1)<1e-1)=Hs_tot(abs(Nexc-n1)<1e-1,abs(Nexc-n1)<1e-1)...
            +Hs_pert(abs(Nexc-n1)<1e-1,abs(Nexc-n2)<1e-1)*Hs_pert(abs(Nexc-n2)<1e-1,abs(Nexc-n1)<1e-1)/(omegaph*(n1-n2));
        Hs_tot(abs(Nexc-n2)<1e-1,abs(Nexc-n2)<1e-1)=Hs_tot(abs(Nexc-n2)<1e-1,abs(Nexc-n2)<1e-1)...
            -Hs_pert(abs(Nexc-n2)<1e-1,abs(Nexc-n1)<1e-1)*Hs_pert(abs(Nexc-n1)<1e-1,abs(Nexc-n2)<1e-1)/(omegaph*(n1-n2));
    end
end

groundenergy_tot=Hs_tot(1,1);
Hs_tot=Hs_tot-groundenergy_tot*eye(dimB*dimB);

dvec_tot=diag(Hs_tot);

temp=dvec_tot(abs(Nexc-1)<1e-1);
Delta1=temp(2)-omegaph;
Delta2=temp(1)-omegaph;

temp=dvec_tot(abs(Nexc-2)<1e-1);
chi1=2*Delta1-temp(3)+2*omegaph;
chi2=2*Delta2-temp(1)+2*omegaph;
chi12=Delta1+Delta2-temp(2)+2*omegaph;

temp=Hs_tot(abs(Nexc-1)<1e-1,abs(Nexc-1)<1e-1);
J=temp(1,2);


end

