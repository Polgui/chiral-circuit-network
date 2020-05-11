r1=0.2;
r2=0.2;
gamma1=1;
gamma2=1;
delta1=0;
delta2=0;

Jvec=linspace(-2,2,100);
phivec=linspace(0,2*pi,101);

beta_dir=zeros(length(Jvec),length(phivec));

t=linspace(0,50,1000);
dt=t(2)-t(1);

for iterJ=1:length(Jvec)
    iterJ/length(Jvec)
    
    for iterphi=1:length(phivec)

        J=Jvec(iterJ);
        phi=phivec(iterphi);

        [ft1,ft2]=ft_function(r1,r2,gamma1,gamma2,delta1,delta2,J,phi,t);
        ftr=exp(1j*phi)*sqrt(gamma1)*(ft1+r2*ft2)+sqrt(gamma2)*(ft2+r1*ft1);
        ftl=sqrt(gamma1)*(ft1+r2*ft2)+exp(1j*phi)*sqrt(gamma2)*(ft2+r1*ft1);
        
        FR=trapz(abs(ftr).^2)*dt;
        FL=trapz(abs(ftl).^2)*dt;
        
        beta_dir(iterJ,iterphi)=FR/(FR+FL);

    end
end

%%

f=figure;
[C1,h1]=contourf(phivec,Jvec,beta_dir,1000);
colormap(pink)
set(h1,'LineColor','none')
hold on
[C2,h2]=contour(phivec,Jvec,beta_dir,[0.75,0.9,0.99],'-k','Linewidth',1.5);
caxis([0,1])
cbar=colorbar;
cbar.LineWidth=1.5;
cbar.TickLength=0.02;
cbar.Ticks=[0,0.5,0.75,0.9,0.99];
cbar.TickLabels{1}='$0$';
cbar.TickLabels{2}='$0.5$';
cbar.TickLabels{3}='$0.75$';
cbar.TickLabels{4}='$0.9$';
cbar.TickLabels{5}='$0.99$';
%title(cbar,'$\beta_{dir}$','Interpreter','latex');
cbar.TickLabelInterpreter='latex';
set(gca,'FontSize',30,'LineWidth',1.5)
ax=gca;
ax.TickLength=[0.02,0.2];
ax.TickLabelInterpreter='latex';
ax.XLabel.Interpreter='latex';
ax.XLabel.String='$\phi$';
ax.XTick=[0,pi/2,pi,3*pi/2,2*pi];
ax.XTickLabels{1}=0;
ax.XTickLabels{2}='$\pi/2$';
ax.XTickLabels{3}='$\pi$';
ax.XTickLabels{4}='$3\pi/2$';
ax.XTickLabels{5}='$2\pi$';
ax.YLabel.Interpreter='latex';
ax.YLabel.String='$J/\gamma$';

