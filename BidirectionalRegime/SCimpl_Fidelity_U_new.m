clearvars
addpath('QSSE_functions/others/')

tic

U1vec=linspace(0,4,501);
U2vec=linspace(0,4,501);

Fid=zeros(length(U1vec),length(U2vec));

for iterU1=1:length(U1vec)
    iterU1/length(U1vec)
    for iterU2=1:length(U2vec)

        phitilde=0;
        
        U1=U1vec(iterU1);
        U2=U2vec(iterU2);
        
        N=4;
        
        J=1;
        
        phi=-pi/2;
        
        id=eye(2*N,2*N);
        
        omega=1e-5;

        plus=[1;1]/sqrt(2);

        Plus=KroneckerProduct(plus,N);

        Plus_f0=(KroneckerProduct([1;1],N)-KroneckerProduct(-[-1;1],N)).*Plus;
        Plus_f0=Plus_f0/norm(Plus_f0);
        Plus_f1=(KroneckerProduct([1;1],N)+KroneckerProduct(-[-1;1],N)).*Plus;
        Plus_f1=Plus_f1/norm(Plus_f1);

        sz=[-1,0;0,1];


        gamma=1;
        
        delta=zeros(1,N)-gamma+1e-10;
        
        U=2*gamma;

        StotR0=zeros(2^N,1);
        StotR1=zeros(2^N,1);
        StotL0=zeros(2^N,1);
        StotL1=zeros(2^N,1);
        
        Uh=[1,-1;1,1]/sqrt(2);

        % A matrix part generic for all qubit states
        A=zeros(2*N,2*N);
        % B matrix part qubit-dependent frequency shift
        B=zeros(2*N,2*N);
        for j=1:N
            A(2*(j-1)+(1:2),2*(j-1)+(1:2))=[gamma-1j*delta(j),1j*J+gamma*exp(1j*phi);1j*J+gamma*exp(1j*phi),gamma-1j*delta(j)];
            B(2*(j-1)+(1:2),2*(j-1)+(1:2))=[-1j*U1,0;0,-1j*U2];

            for k=1:j-1
                A(2*(j-1)+(1:2),2*(k-1)+(1:2))=gamma*exp(1j*phitilde*(j-k))*[1,exp(-1j*phi);exp(1j*phi),1];
            end

             for k=j+1:N
                A(2*(j-1)+(1:2),2*(k-1)+(1:2))=gamma*exp(1j*phitilde*(k-j))*[1,exp(1j*phi);exp(-1j*phi),1];
            end

        end
        
        bin=zeros(2*N,1);
        for j=1:N
            bin(2*(j-1)+(1:2))=sqrt(gamma)*exp(1j*phitilde*(j-1))*[1;exp(1j*phi)];
        end
        bout=zeros(1,2*N);
        for j=1:N
            bout(2*(j-1)+(1:2))=sqrt(gamma)*exp(1j*phitilde*(N-j+1))*[1,exp(-1j*phi)];
        end

        for indexcomb=1:2^N
            u=dec2bin(indexcomb-1)-'0';
            s=zeros(1,N); % qubit state vector in computational basis (0 or 1)
            s(N-length(u)+1:N)=u;
            
            s=repelem(s,2);
            
            A_red=A+B.*(s'*s);
            bin_red=bin;
            bout_red=bout;

            invA=linsolve(1j*omega*id-A_red,bin_red);

            phicorr=-N*pi/2;
            SR=exp(1j*phitilde*N)*Uh*[1,0;0,exp(1j*phicorr)]*Uh;

            SR0=SR(1,1)+Uh(1,2)*exp(1j*phicorr)*bout_red*invA*Uh(2,1);
            StotR0(indexcomb)=SR0;
            SR1=SR(2,1)+Uh(2,2)*exp(1j*phicorr)*bout_red*invA*Uh(2,1);
            StotR1(indexcomb)=SR1;
        end
        Fid(iterU1,iterU2)=abs(Plus_f0'*(StotR0.*Plus))^2+abs(Plus_f1'*(StotR1.*Plus))^2;
    end
    
end

toc

%% Plots


myBlue=[0.4,0.6,0.8];

h=figure;
box on

[C1,h1]=contourf(U1vec/2,U2vec/2,Fid',500);
hold on
[C2,h2]=contour(U1vec/2,U2vec/2,Fid',[0.5,0.75,0.9,0.99],'-k','Linewidth',1.5);
colormap('pink')
set(h1,'LineColor','none')
cbar=colorbar;
caxis([0.1,1])
cbar.LineWidth=1.5;
cbar.TickLength=0.02;
cbar.Ticks=[.5,0.75,0.9,0.99];
title(cbar,'$\mathcal F_{\mathcal Z}$','Interpreter','latex');
cbar.TickLabelInterpreter='latex';
set(gca,'FontSize',30)
ax=gca;
ax.TickLabelInterpreter='latex';
ax.YLabel.Interpreter='latex';
ax.YLabel.String='$V_2^n/\gamma_r$';
ax.XLabel.Interpreter='latex';
ax.XLabel.String='$V_1^n/\gamma_r$';
ax.XTick=[0,1,2];
ax.YTick=[0,1,2];
