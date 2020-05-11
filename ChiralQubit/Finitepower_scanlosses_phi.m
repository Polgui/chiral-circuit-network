clearvars
addpath('functions')


gamma=1;
gamma1=gamma;
gamma2=gamma;
r=0.2;

phi=imag(log(-1j*(-1j+r)/(1j+r)));
J=-(1+r^2)*sin(phi)*gamma;

gammar=2*gamma*(1+2*r*cos(phi)+r^2);
Omegavec=gammar*[0.1,0.5,1,2,5];

gammaphivec=10.^linspace(-4,0,51);
Plvec=zeros(length(Omegavec),length(gammaphivec));
Prvec=zeros(length(Omegavec),length(gammaphivec));
gammap=1e-2;

delta1=0;
delta2=0;

chi12=50;
chi=50;

dimB=4;

a=diag(sqrt(1:dimB-1),1);
idB=eye(dimB);

L1=sqrt(gamma1)*(KroneckerProduct(a,idB)+r*KroneckerProduct(idB,a));
L2=sqrt(gamma2)*(r*KroneckerProduct(a,idB)+KroneckerProduct(idB,a));

LR=exp(1j*phi)*L1+L2;
LL=L1+exp(1j*phi)*L2;


for iterOmega=1:length(Omegavec)
    
    Omega=Omegavec(iterOmega)

    sml=(KroneckerProduct(a,idB)+1j*KroneckerProduct(idB,a))/sqrt(2);
    smr=(1j*KroneckerProduct(a,idB)+KroneckerProduct(idB,a))/sqrt(2);
    
    Htot=-chi12*KroneckerProduct(a'*a,a'*a)-chi*KroneckerProduct(a'*a'*a*a,idB)-chi*KroneckerProduct(idB,a'*a'*a*a)...
        -delta1*KroneckerProduct(a'*a,idB)-delta2*KroneckerProduct(idB,a'*a)...
        +J*KroneckerProduct(a',a)+J*KroneckerProduct(a,a')...
        +sin(phi)*(L2'*L1+L1'*L2)...
        +(Omega/sqrt(gammar))*(LR+LR');


    DL=SOpre(LL)*SOpost(LL')-1/2*SOpre(LL'*LL)-1/2*SOpost(LL'*LL);
    DR=SOpre(LR)*SOpost(LR')-1/2*SOpre(LR'*LR)-1/2*SOpost(LR'*LR);


    for itergammaphi=1:length(gammaphivec)
        gammaphi=gammaphivec(itergammaphi);


        sz1=KroneckerProduct(a'*a,idB);
        sz2=KroneckerProduct(idB,a'*a);
        DPhi=2*gammaphi*(SOpre(sz1)*SOpost(sz1)-1/2*SOpre(sz1'*sz1)-1/2*SOpost(sz1'*sz1)...
            +SOpre(sz2)*SOpost(sz2)-1/2*SOpre(sz2'*sz2)-1/2*SOpost(sz2'*sz2));

        sm1=KroneckerProduct(a,idB);
        sm2=KroneckerProduct(idB,a);
        Dp=gammap*(SOpre(sm1)*SOpost(sm1)-1/2*SOpre(sm1'*sm1)-1/2*SOpost(sm1'*sm1)...
            +SOpre(sm2)*SOpost(sm2)-1/2*SOpre(sm2'*sm2)-1/2*SOpost(sm2'*sm2));

        Liouvillian=-1i*SOpre(Htot)+1i*SOpost(Htot)+DL+DR+DPhi+Dp;


        [V,D]=eig(Liouvillian);
        d=diag(D);
        dr=real(d);
        [~,index]=max(real(dr));

        rhoSvec=V(:,index);
        rho=vec2mat(rhoSvec);
        rho=(rho+rho')/2;
        rho=rho/trace(rho);
        
        Prvec(iterOmega,itergammaphi)=trace(rho*LR'*LR);
        Plvec(iterOmega,itergammaphi)=trace(rho*LL'*LL);
    end
end

%%
myBlue=[0.4,0.6,0.8];

figure
box on
for iterOmega=1:length(Omegavec)
    h=loglog(gammaphivec/gammar,Plvec(iterOmega,:)./(Prvec(iterOmega,:)));
    h.LineWidth=3;
    h.Color=myBlue*(length(Omegavec)-iterOmega+1)/length(Omegavec);
    hold on
    
    set(gca,'FontSize',30)
    ax=gca;
    ax.YLim=[0,0.501];
    ax.TickLabelInterpreter='latex';
    ax.LineWidth=1;
end

h=loglog(gammaphivec,gammaphivec,'--r');
h.LineWidth=3;

% pleg=legend;
% pleg.Interpreter='latex';
% % pleg.FontSize=22;
% for iterOmega=1:length(Omegavec)
%     pleg.String{iterOmega}=['$\Omega/\gamma_r =$ ',num2str(Omegavec(iterOmega)/gammar)];
% end
% 
% pleg.String{end}='';

legstr=cell(1,length(Omegavec));
for iterOmega=1:length(Omegavec)
     legstr{iterOmega}=['$\Omega/\gamma_r =$ ',num2str(Omegavec(iterOmega)/gammar)];
end
legend(legstr,'Interpreter','latex')
pleg=legend;

set(pleg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.9]));

% [legend_h,object_h]=columnlegend(2, legstr, 'location','SouthEast','padding',-0.2);

ax.YLim=[1e-4,1e0];
ax.XLim=[1e-4,1e0];
ax.XLabel.Interpreter='latex';
ax.XLabel.String='$\gamma_\varphi/\gamma_r$';
ax.YLabel.Interpreter='latex';
ax.YLabel.String='$\langle \hat L_L^\dagger \hat L_L \rangle/\langle \hat L_R^\dagger \hat L_R \rangle$';


% drawnow;
% set(findall(object_h, 'type', 'text'), 'interpreter', 'latex')
% set(findall(legend_h, 'type', 'text'), 'interpreter', 'latex')