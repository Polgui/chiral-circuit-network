clearvars
addpath('functions')

Omega=1;

gamma=1;
gamma1=gamma;
gamma2=gamma;
r=0.2;

phi=imag(log(-1j*(-1j+r)/(1j+r)));
J=-(1+r^2)*sin(phi)*gamma;

gammar=2*gamma*(1+2*r*cos(phi)+r^2);
    
gammap=1e-2;
gammaphi=1e-2;

chi=300;
chi12=50;

delta=2*r*gamma*sin(phi);

phitilde=0;

dt=0.1/chi;
dt=10/500;
tspan=0:dt:20;
Nt=length(tspan);

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

PL=zeros(1,Nt);
PR=zeros(1,Nt);
PL(1)=trace(LLtot'*LLtot*rhoS0);
PR(1)=trace(LRtot'*LRtot*rhoS0);

SubPurity2=zeros(1,Nt);
SubPurity2(1)=trace(TrX(rhoS0,[1,2],[dimB,dimB,dimB,dimB])^2);
SubPurity1=zeros(1,Nt);
SubPurity1(1)=trace(TrX(rhoS0,[3,4],[dimB,dimB,dimB,dimB])^2);

Purity=zeros(1,Nt);
Purity(1)=trace(rhoS0^2);

Nop=KroneckerProduct(a'*a,idB,idB,idB)+KroneckerProduct(idB,a'*a,idB,idB)+KroneckerProduct(idB,idB,a'*a,idB)+KroneckerProduct(idB,idB,idB,a'*a);
Nexc=zeros(1,Nt);
Nexc(1)=trace(Nop*rhoS0);

gammar=2*gamma*(1+2*r*cos(phi)+r^2);
smr=(1j*KroneckerProduct(a,idB)+KroneckerProduct(idB,a))/sqrt(2);
sml=(KroneckerProduct(a,idB)+1j*KroneckerProduct(idB,a))/sqrt(2);

psidimer=(eye(dimB^4)-1j*2*sqrt(2)*exp(-1j*angle(1-r*1j*exp(-2*1j*atan(1/r))))*Omega/sqrt(gammar)*(KroneckerProduct(smr',idB,idB)-exp(1j*phitilde)*KroneckerProduct(idB,idB,smr'))/sqrt(2))*InitialState;
psidimer=psidimer/norm(psidimer);
overlap=zeros(1,Nt);
overlap(1)=trace(psidimer*psidimer'*rhoS0);

rhovect=zeros(length(rhoSvec0),Nt+1);
rhovect(:,1)=rhoSvec0;

counter=0;
for k=1:Nt-1
    counter=counter+1;
    if counter>Nt/100
        k/Nt
        counter=0;
    end
        
    rhovect(:,k+1)=expv(dt,Liouvillian,rhovect(:,k));
    rho=vec2mat(rhovect(:,k+1));
    rho=(rho+rho')./2;
    rho=rho./trace(rho);
    PL(k+1)=trace(LLtot'*LLtot*rho);
    PR(k+1)=trace(LRtot'*LRtot*rho);
    Purity(k+1)=trace(rho^2);
    SubPurity2(k+1)=trace(TrX(rho,[1,2],[dimB,dimB,dimB,dimB])^2);
    SubPurity1(k+1)=trace(TrX(rho,[3,4],[dimB,dimB,dimB,dimB])^2);
    Nexc(k+1)=trace(rho*Nop);
    
   
    smr=(1j*KroneckerProduct(a,idB)+KroneckerProduct(idB,a))/sqrt(2);
    sml=(KroneckerProduct(a,idB)+1j*KroneckerProduct(idB,a))/sqrt(2);

    psidimer=(eye(dimB^4)-2*sqrt(2)*exp(-1j*angle(1-r*1j*exp(-2*1j*atan(1/r))))*(Omega/gammar)*(KroneckerProduct(smr',idB,idB)-exp(1j*phitilde)*KroneckerProduct(idB,idB,smr'))/sqrt(2))*InitialState;
    psidimer=psidimer/norm(psidimer);
    overlap(k+1)=trace(psidimer*psidimer'*rho);
end

%%

myBlue=[0.4,0.6,0.8];
myGreen=[0.4,0.8,0.6];

figure
box on
hold on

h=plot(tspan,SubPurity1,'k');
h.LineWidth=3;
h=plot(tspan,SubPurity2,'--k');
h.LineWidth=3;

h=plot(tspan,Nexc);
h.LineWidth=3;
h.Color=myGreen;

%h=plot(tspan,PL/2,'--');
%h.LineWidth=3;
%h.Color=myBlue;    
h=plot(tspan,PR);
h.LineWidth=3;
h.Color=myBlue;



h=plot(tspan,Purity,'r');
h.LineWidth=3;
h=plot(tspan,overlap,'--r');
h.LineWidth=3;




set(gca,'FontSize',30)
ax=gca;
ax.YLim=[0,1];
ax.TickLabelInterpreter='latex';
ax.LineWidth=1;


ax.XLabel.Interpreter='latex';
ax.XLabel.String='$\gamma t$';
ax.YLabel.Interpreter='latex';

pleg=legend;
pleg.Interpreter='latex';

pleg.String{1}='Tr$(\hat \rho_{(1)}^2)$';
pleg.String{2}='Tr$(\hat \rho_{(2)}^2)$';
pleg.String{3}='$\sum_{n,k}\langle (\hat a^n_k)^\dagger \hat a^n_k\rangle$';
pleg.String{4}='$\langle \hat L_R^\dagger \hat L_R \rangle$';
pleg.String{5}='Tr$(\hat \rho^2)$';
pleg.String{6}='$\langle{D}|\hat \rho|{D}\rangle$';


reorderLegend([1,2,4,6,5,3])
set(pleg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.9]));