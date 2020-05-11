clearvars
addpath('functions')

gamma=1;
gamma1=gamma;
gamma2=gamma;
r=0.2;

phi=imag(log(-1j*(-1j+r)/(1j+r)));
J=-(1+r^2)*sin(phi)*gamma;
gammar=2*gamma*(1+2*r*cos(phi)+r^2);

Omegavec=gamma*(10.^linspace(-1,1,50));

chi12vec=10.^linspace(0,4,51);
Plvec=zeros(length(Omegavec),length(chi12vec));
Prvec=zeros(length(Omegavec),length(chi12vec));
ratio_nl=zeros(length(Omegavec),length(chi12vec));

gammap=1e-2;
gammaphi=1e-2;


delta1=0;
delta2=0;

chi=50;

dimB=4;

a=diag(sqrt(1:dimB-1),1);
idB=eye(dimB);

L1=sqrt(gamma1)*(KroneckerProduct(a,idB)+r*KroneckerProduct(idB,a));
L2=sqrt(gamma2)*(r*KroneckerProduct(a,idB)+KroneckerProduct(idB,a));

LR=exp(1j*phi)*L1+L2;
LL=L1+exp(1j*phi)*L2;


for iterOmega=1:length(Omegavec)
    iterOmega/length(Omegavec)
    
    Omega=Omegavec(iterOmega);

    sml=(KroneckerProduct(a,idB)+1j*KroneckerProduct(idB,a))/sqrt(2);
    smr=(1j*KroneckerProduct(a,idB)+KroneckerProduct(idB,a))/sqrt(2);
    
    
    DL=SOpre(LL)*SOpost(LL')-1/2*SOpre(LL'*LL)-1/2*SOpost(LL'*LL);
    DR=SOpre(LR)*SOpost(LR')-1/2*SOpre(LR'*LR)-1/2*SOpost(LR'*LR);


    for iterchi12=1:length(chi12vec)
        
        chi12=chi12vec(iterchi12);
        
        Htot=-chi12*KroneckerProduct(a'*a,a'*a)-chi*KroneckerProduct(a'*a'*a*a,idB)-chi*KroneckerProduct(idB,a'*a'*a*a)...
        -delta1*KroneckerProduct(a'*a,idB)-delta2*KroneckerProduct(idB,a'*a)...
        +J*KroneckerProduct(a',a)+J*KroneckerProduct(a,a')...
        +sin(phi)*(L2'*L1+L1'*L2)...
        +(Omega/sqrt(gammar))*(LR+LR');
    


        sz1=KroneckerProduct(diag([-1,ones(1,dimB-1)]),idB);
        sz2=KroneckerProduct(idB,diag([-1,ones(1,dimB-1)]));
        DPhi=gammaphi*(SOpre(sz1)*SOpost(sz1)-1/2*SOpre(sz1'*sz1)-1/2*SOpost(sz1'*sz1)...
            +SOpre(sz2)*SOpost(sz2)-1/2*SOpre(sz2'*sz2)-1/2*SOpost(sz2'*sz2));

        sm1=KroneckerProduct(a,idB);
        sm2=KroneckerProduct(idB,a);
        Dp=gammap*(SOpre(sm1)*SOpost(sm1)-1/2*SOpre(sm1'*sm1)-1/2*SOpost(sm1'*sm1)...
            +SOpre(sm2)*SOpost(sm2)-1/2*SOpre(sm2'*sm2)-1/2*SOpost(sm2'*sm2));

        Liouvillian=-1i*SOpre(Htot)+1i*SOpost(Htot)+DL+DR+DPhi+Dp;


        [V,D]=eig(Liouvillian);
        d=diag(D);
        dr=real(d);
        [drsort,index]=sort(dr);
        if abs(drsort(end)-drsort(end-1))<0.2
            error 'many steady-states?'
        end

        rhoSvec=V(:,index(end));
        rho=vec2mat(rhoSvec);
        rho=(rho+rho')/2;
        rho=rho/trace(rho);
        
        Prvec(iterOmega,iterchi12)=trace(rho*LR'*LR);
        Plvec(iterOmega,iterchi12)=trace(rho*LL'*LL);
        
        op0=zeros(dimB,dimB);
        op0(1,1)=1;
        op1=zeros(dimB,dimB);
        op1(2,2)=1;
        ratio_nl(iterOmega,iterchi12)=1-trace(KroneckerProduct(op1,op0)*rho+KroneckerProduct(op0,op1)*rho+KroneckerProduct(op0,op0)*rho);
    end
end

%%

figure
[C1,h1]=contourf(chi12vec/gamma,Omegavec/gamma,log10(real(Plvec./(Prvec))),1000);
colormap(flipud(pink))
set(h1,'LineColor','none')
caxis([-2,0])
hold on
%[C2,h2]=contour(chi12vec,Omegavec,real(Plvec./(Prvec)),[1e-2,1e-1],'-k','Linewidth',1.5);
[C3,h3]=contour(chi12vec,Omegavec,real(ratio_nl),[1e0,1e-2,1e-3],'--r','Linewidth',1.5);
h4=plot([100,100],[.1,10],'--','Linewidth',1.5,'color',[0,0,0]+0.5);
cbar=colorbar;
cbar.LineWidth=1.5;
cbar.TickLength=0.02;
cbar.Ticks=[-2,-1,0];
cbar.TickLabels{1}='$10^{-2}$';
cbar.TickLabels{2}='$10^{-1}$';
cbar.TickLabels{3}='$10^{0}$';
cbar.TickLabelInterpreter='latex';
set(gca,'FontSize',30,'LineWidth',1.5)
ax=gca;
ax.TickLength=[0.02,0.2];
ax.TickLabelInterpreter='latex';
ax.XLabel.Interpreter='latex';
ax.YLabel.Interpreter='latex';
ax.YLabel.String='$\Omega/\gamma$';
ax.XLabel.String='$\chi/\gamma$';
ax.XScale='log';
ax.YScale='log';

%%

op=zeros(dimB,dimB);
op(dimB,dimB)=1;
trace(KroneckerProduct(op,idB)*rho)