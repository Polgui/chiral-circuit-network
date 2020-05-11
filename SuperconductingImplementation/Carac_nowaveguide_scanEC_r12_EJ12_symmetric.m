clearvars
addpath('functions')

% Frequency units is 1MHZ

omegaph=8e3;

r12vec=linspace(0.1,0.5,21);
ECvec=linspace(200,400,21);

EJ12opt=zeros(length(ECvec),length(r12vec));
EJopt=zeros(length(ECvec),length(r12vec));
Jvec=zeros(length(ECvec),length(r12vec));
Delta=zeros(length(ECvec),length(r12vec));
chi=zeros(length(ECvec),length(r12vec));
chi12=zeros(length(ECvec),length(r12vec));
Jopt=zeros(length(ECvec),length(r12vec));
Deltaopt=zeros(length(ECvec),length(r12vec));

for iterEC=1:length(ECvec)
    
    iterEC
    
EC=ECvec(iterEC); % e^2/(2C1)

EJ=omegaph^2/(8*EC);

options = optimset('TolX',1e-3,'TolFun',1e-3);

EJ12=EJ*r12vec(1);
    
    for iter=1:length(r12vec)
        
        r12=r12vec(iter);

        dimB=5;

        gamma=1;

        Cmat=[1/EC+r12/EC, -r12/EC;-r12/EC, 1/EC+r12/EC];
        Cmat_inv=pinv(Cmat);

        reff=Cmat_inv(2,1)/Cmat_inv(1,1);

        phi=-pi/2-2*atan(1/reff);

        Jopt(iterEC,iter)=-gamma*(1+reff^2)*sin(phi);

        Deltaopt(iterEC,iter)=2*reff*gamma*sin(phi);

        fun = @(x)optimize_EJ12_symmetric(x,[EJ,EC,r12,omegaph,dimB,Deltaopt(iterEC,iter),Jopt(iterEC,iter)]);
        x0 = EJ12;
        x = fminsearch(fun,x0,options);
        EJ12opt(iterEC,iter)=x;
        EJ12=x;

        fun = @(x)optimize_frequencies_symmetric(x,[EJ12opt(iterEC,iter),EC,r12,omegaph,dimB,Deltaopt(iterEC,iter)]);
        x0 = EJ;
        x = fminsearch(fun,x0,options);
        EJopt(iterEC,iter)=x;
        EJ=x;

        [Jvec(iterEC,iter),Delta(iterEC,iter),~,chi(iterEC,iter),~,chi12(iterEC,iter)]=find_frequencies(x,x,EJ12opt(iterEC,iter),EC,EC,r12,omegaph,dimB,0);

    end
end


%%


figure
contourf(r12vec,ECvec,chi12)
colorbar
figure
contourf(r12vec,ECvec,chi12./chi)
colorbar

%% 

myBlue=[0.4,0.6,0.8];
myRed=[0.8,0.4,0.4];

figure
h=loglog(10.^linspace(1,3,100),0.25./(10.^linspace(1,3,100)).^(0.45),'--');
h.LineWidth=3;
h.Color=myRed;
hold on
h=loglog(10.^linspace(1,3,100),0.25./(10.^linspace(1,3,100)).^(0.65),'--');
h.LineWidth=3;
h.Color=myBlue;

for n=10:-1:5
    omegaph=n*1e3
    load(['omegaph',num2str(n),'e3.mat'])

    reffopt=zeros(1,length(ECvec));
    chi12max=zeros(1,length(ECvec));
    chimax=zeros(1,length(ECvec));
    EJeff=zeros(1,length(ECvec));
    ECeff=zeros(1,length(ECvec));


    for iterEC=1:length(ECvec)

        EC=ECvec(iterEC);


        [chi12max(iterEC),index]=max(chi12(iterEC,:));
        chimax(iterEC)=chi(iterEC,index);
        EJeff(iterEC)=EJopt(iterEC,index)+EJ12opt(iterEC,index);

        r12=r12vec(index);
        Cmat=[1/EC+r12/EC, -r12/EC;-r12/EC, 1/EC+r12/EC];
        Cmat_inv=pinv(Cmat);

        reffopt(iterEC)=Cmat_inv(2,1)/Cmat_inv(1,1);
        ECeff(iterEC)=Cmat_inv(1,1);
    end

    h=loglog(EJeff./ECeff,chi12max/omegaph);
    h.LineWidth=6;
    h.Color=myBlue*(6-(n-4)+1)/6;
    
    h=loglog(EJeff./ECeff,chimax/omegaph);
    h.LineWidth=6;
    h.Color=myRed*(6-(n-4)+1)/6;
 
end

set(gca,'FontSize',30)
ax=gca;
ax.XLim=[30,1000];
ax.YLim=[0,0.1];
ax.TickLabelInterpreter='latex';
ax.LineWidth=1;
ax.XLabel.Interpreter='latex';
ax.XLabel.String='$\overline E_J^k/E_C^k$';
ax.YLabel.Interpreter='latex';
ax.YLabel.String='$(U,\chi_k)/\omega_0$';

%%
myBlue=[0.4,0.6,0.8];
myRed=[0.8,0.4,0.4];

n=8;
omegaph=n*1e3;
load(['omegaph',num2str(n),'e3.mat'])

ECeff=zeros(length(ECvec),length(r12vec));
reff=zeros(length(ECvec),length(r12vec));
for iterEC=1:length(ECvec)
    for iterr12=1:length(r12vec)
        r12=r12vec(iterr12);
        EC=ECvec(iterEC);
        
        Cmat=[1/EC+r12/EC, -r12/EC;-r12/EC, 1/EC+r12/EC];
        Cmat_inv=pinv(Cmat);

        reff(iterEC,iterr12)=Cmat_inv(2,1)/Cmat_inv(1,1);
        ECeff(iterEC,iterr12)=Cmat_inv(1,1);
    end
end

figure
box on
for iterEC=length(ECvec):-3:1
    h=plot(reff(iterEC,:),chi12(iterEC,:)./ECeff(iterEC,:));
    hold on
    h.LineWidth=3;
    h.Color=myBlue*(length(ECvec)-iterEC+1)/length(ECvec);

    h=plot(reff(iterEC,:),chi(iterEC,:)./ECeff(iterEC,:));
    h.LineWidth=3;
    h.Color=myRed*(length(ECvec)-iterEC+1)/length(ECvec);
end

set(gca,'FontSize',30)
ax=gca;
ax.YLim=[0,1.001];
ax.TickLabelInterpreter='latex';
ax.LineWidth=1;
ax.XLabel.Interpreter='latex';
ax.XLabel.String='$r_k$';
ax.YLabel.Interpreter='latex';
ax.YLabel.String='$(U,\chi_k)/E_C^k$';
