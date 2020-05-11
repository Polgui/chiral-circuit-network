clearvars
addpath('functions')

gammaloss=1e-2;

gamma=1;
gamma1=gamma;
gamma2=gamma;
r=0.2;

phi=imag(log(-1j*(-1j+r)/(1j+r)));
J=-(1+r^2)*sin(phi)*gamma;
        
gammar=2*gamma*(1+2*r*cos(phi)+r^2);

chi12vec=gamma*(10.^linspace(1,3,31));

Omegavec=gamma*(10.^linspace(0,1,21));

overlap=zeros(length(Omegavec),length(chi12vec));

for iterchi12=1:length(chi12vec)
    iterchi12
    
    chi12=chi12vec(iterchi12);

    for iterOmega=1:length(Omegavec)

        Omega=Omegavec(iterOmega);

        gammap=gammaloss;
        gammaphi=gammaloss;

        chi=300;

        delta=2*r*gamma*sin(phi);

        phitilde=0;

        dimB=2;

        a=diag(sqrt(1:dimB-1),1);
        idB=eye(dimB);

        L1=sqrt(gamma1)*(KroneckerProduct(a,idB)+r*KroneckerProduct(idB,a));
        L2=sqrt(gamma2)*(r*KroneckerProduct(a,idB)+KroneckerProduct(idB,a));

        LR=exp(1j*phi)*L1+L2;
        LL=L1+exp(1j*phi)*L2;

        H1=-chi12*KroneckerProduct(a'*a,a'*a)-(chi/2)*KroneckerProduct(a'*a'*a*a,idB)-(chi/2)*KroneckerProduct(idB,a'*a'*a*a)...
            -delta*KroneckerProduct(a'*a,idB)-delta*KroneckerProduct(idB,a'*a)...
            +J*KroneckerProduct(a',a)+J*KroneckerProduct(a,a')...
            +sin(phi)*(L2'*L1+L1'*L2)...
            +1j*(Omega/sqrt(gammar))*LR-1j*(Omega/sqrt(gammar))*LR';

        H2=-chi12*KroneckerProduct(a'*a,a'*a)-(chi/2)*KroneckerProduct(a'*a'*a*a,idB)-(chi/2)*KroneckerProduct(idB,a'*a'*a*a)...
            -delta*KroneckerProduct(a'*a,idB)-delta*KroneckerProduct(idB,a'*a)...
            +J*KroneckerProduct(a',a)+J*KroneckerProduct(a,a')...
            +sin(phi)*(L2'*L1+L1'*L2)...
            +1j*(Omega/sqrt(gammar))*exp(-1j*phitilde)*LR-1j*(Omega/sqrt(gammar))*exp(1j*phitilde)*LR';

        HL=-0.5*1j*KroneckerProduct(LL',LL)*exp(1j*phitilde);
        HL=HL+HL';

        HR=-0.5*1j*KroneckerProduct(LR,LR')*exp(1j*phitilde);
        HR=HR+HR';

        Htot=KroneckerProduct(H1,idB,idB)+KroneckerProduct(idB,idB,H2)+HL+HR;

        LLtot=KroneckerProduct(LL,idB,idB)+exp(1j*phitilde)*KroneckerProduct(idB,idB,LL);
        LRtot=KroneckerProduct(idB,idB,LR)+exp(1j*phitilde)*KroneckerProduct(LR,idB,idB);

        InitialState=KroneckerProduct([1,zeros(1,dimB-1)],4)';

        rhoS0=InitialState*InitialState';
        rhoSvec0=mat2vec(rhoS0);

        DL=SOpre(LLtot)*SOpost(LLtot')-1/2*SOpre(LLtot'*LLtot)-1/2*SOpost(LLtot'*LLtot);
        DR=SOpre(LRtot)*SOpost(LRtot')-1/2*SOpre(LRtot'*LRtot)-1/2*SOpost(LRtot'*LRtot);

        sp1=sqrt(gammap)*KroneckerProduct(a,idB,idB,idB);
        sp2=sqrt(gammap)*KroneckerProduct(idB,a,idB,idB);
        sp3=sqrt(gammap)*KroneckerProduct(idB,idB,a,idB);
        sp4=sqrt(gammap)*KroneckerProduct(idB,idB,idB,a);

        Dp1=SOpre(sp1)*SOpost(sp1')-1/2*SOpre(sp1'*sp1)-1/2*SOpost(sp1'*sp1);
        Dp2=SOpre(sp2)*SOpost(sp2')-1/2*SOpre(sp2'*sp2)-1/2*SOpost(sp2'*sp2);
        Dp3=SOpre(sp3)*SOpost(sp3')-1/2*SOpre(sp3'*sp3)-1/2*SOpost(sp3'*sp3);
        Dp4=SOpre(sp4)*SOpost(sp4')-1/2*SOpre(sp4'*sp4)-1/2*SOpost(sp4'*sp4);

        sphi1=sqrt(2*gammaphi)*KroneckerProduct(a'*a,idB,idB,idB);
        sphi2=sqrt(2*gammaphi)*KroneckerProduct(idB,a'*a,idB,idB);
        sphi3=sqrt(2*gammaphi)*KroneckerProduct(idB,idB,a'*a,idB);
        sphi4=sqrt(2*gammaphi)*KroneckerProduct(idB,idB,idB,a'*a);

        Dphi1=SOpre(sphi1)*SOpost(sphi1')-1/2*SOpre(sphi1'*sphi1)-1/2*SOpost(sphi1'*sphi1);
        Dphi2=SOpre(sphi2)*SOpost(sphi2')-1/2*SOpre(sphi2'*sphi2)-1/2*SOpost(sphi2'*sphi2);
        Dphi3=SOpre(sphi3)*SOpost(sphi3')-1/2*SOpre(sphi3'*sphi3)-1/2*SOpost(sphi3'*sphi3);
        Dphi4=SOpre(sphi4)*SOpost(sphi4')-1/2*SOpre(sphi4'*sphi4)-1/2*SOpost(sphi4'*sphi4);

        Liouvillian=-1i*SOpre(Htot)+1i*SOpost(Htot)+DL+DR+Dp1+Dp2+Dp3+Dp4+Dphi1+Dphi2+Dphi3+Dphi4;

        Nop=KroneckerProduct(a'*a,idB,idB,idB)+KroneckerProduct(idB,a'*a,idB,idB)+KroneckerProduct(idB,idB,a'*a,idB)+KroneckerProduct(idB,idB,idB,a'*a);

        [V,D]=eig(Liouvillian);
        d=diag(D);
        dr=real(d);
        [~,index]=max(real(dr));

        rhoSvec=V(:,index);
        rho=vec2mat(rhoSvec);
        rho=(rho+rho')/2;
        rho=rho/trace(rho);
        
        gammar=2*gamma*(1+2*r*cos(phi)+r^2);
        smr=(1j*KroneckerProduct(a,idB)+KroneckerProduct(idB,a))/sqrt(2);
        sml=(KroneckerProduct(a,idB)+1j*KroneckerProduct(idB,a))/sqrt(2);

        psidimer=(eye(dimB^4)-2*sqrt(2)*exp(-1j*angle(1-r*1j*exp(-2*1j*atan(1/r))))*(Omega/gammar)*(KroneckerProduct(smr',idB,idB)-exp(1j*phitilde)*KroneckerProduct(idB,idB,smr'))/sqrt(2))*InitialState;
        psidimer=psidimer/norm(psidimer);
        overlap(iterOmega,iterchi12)=trace(psidimer*psidimer'*rho);
    end
end

%%

myBlue=[0.4,0.6,0.8];
figure
box on
for iterOmega=1:length(Omegavec)
    h=loglog(chi12vec*gammar/(Omegavec(iterOmega)^2),1-overlap(iterOmega,:));
    hold on
    h.LineWidth=3;
    h.Color=myBlue*(length(Omegavec)-iterOmega+1)/length(Omegavec);  
end
% h=loglog([1e0,1e2],30*[1e0,1e2].^(-2),'--r');
% h.LineWidth=3;

set(gca,'FontSize',30)
ax=gca;
ax.TickLabelInterpreter='latex';
ax.LineWidth=1;
ax.XLim=[1e-1,1e2];
ax.XLabel.Interpreter='latex';
ax.XLabel.String='$\chi\gamma_r/\Omega^2$';
ax.YLabel.Interpreter='latex';
ax.YLabel.String='$1-\langle D|\hat \rho |D\rangle$';