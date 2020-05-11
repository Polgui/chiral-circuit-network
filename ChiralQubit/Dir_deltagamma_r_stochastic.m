Nreal=10000;
delta1=0;
delta2=0;
gamma=1;
rb=0.2;
phi=imag(log(-1j*(-1j+rb)/(1j+rb)));
J=-(1+rb^2)*sin(phi)*gamma;

beta_dir_list=cell(1,Nreal);
    
dgammavec=10.^linspace(-2,0,50);
rvec=10.^linspace(-2,0,51);

t=linspace(0,50,1000);
dt=t(2)-t(1);

parfor iterreal=1:Nreal
    iterreal
    
    beta_dir_stoch=zeros(length(dgammavec),length(rvec));
    
    for iterdgamma=1:length(dgammavec)
        
        for iterr=1:length(rvec)

            r1=rb+rvec(iterr)*(rand-0.5);
            r2=rb+rvec(iterr)*(rand-0.5);

            gamma1=gamma+dgammavec(iterdgamma)*(rand-0.5);
            gamma2=gamma+dgammavec(iterdgamma)*(rand-0.5);

            [ft1,ft2]=ft_function(r1,r2,gamma1,gamma2,delta1,delta2,J,phi,t);
            ftr=exp(1j*phi)*sqrt(gamma1)*(ft1+r2*ft2)+sqrt(gamma2)*(ft2+r1*ft1);
            ftl=sqrt(gamma1)*(ft1+r2*ft2)+exp(1j*phi)*sqrt(gamma2)*(ft2+r1*ft1);

            FR=trapz(abs(ftr).^2)*dt;
            FL=trapz(abs(ftl).^2)*dt;

            beta_dir_stoch(iterdgamma,iterr)=beta_dir_stoch(iterdgamma,iterr)+(1/Nreal)*(FR/(FR+FL));

        end
    end
    
    beta_dir_list{iterreal}=beta_dir_stoch;
    
end

%%
 beta_dir=zeros(length(dgammavec),length(rvec));
 
for iterreal=1:Nreal
    beta_dir=beta_dir+beta_dir_list{iterreal};
end

save('data')

%%

f=figure;
[C1,h1]=contourf(rvec/(2*sqrt(3)),dgammavec/(2*sqrt(3)),log10(1-beta_dir),500);
colormap(flipud(pink))
set(h1,'LineColor','none')
hold on
[C2,h2]=contour(rvec/(2*sqrt(3)),dgammavec/(2*sqrt(3)),1-beta_dir,10.^(-(2:4)),'-k','Linewidth',1.5);
cbar=colorbar;
cbar.LineWidth=1.5;
cbar.TickLength=0.02;
cbar.Ticks=-(4:-1:0);
cbar.TickLabels{1}='$10^{-4}$';
cbar.TickLabels{2}='$10^{-3}$';
cbar.TickLabels{3}='$10^{-2}$';
%title(cbar,'$1-\beta_{dir}$','Interpreter','latex');
cbar.TickLabelInterpreter='latex';
set(gca,'FontSize',30,'LineWidth',1.5)
ax=gca;
ax.TickLength=[0.02,0.2];
ax.XTick=10.^[-2,-1,0];
ax.TickLabelInterpreter='latex';
ax.XLabel.Interpreter='latex';
ax.XLabel.String='$\delta r$';
ax.YLabel.Interpreter='latex';
ax.YLabel.String='$\delta\gamma/\gamma$';
ax.XScale='log';
ax.YScale='log';
ax.XLim=[ax.XLim(1),0.1];
cbar.Ticks=[-4,-3,-2];
caxis([-4.6,-1.3]);