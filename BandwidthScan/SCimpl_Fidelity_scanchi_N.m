% Settings

clearvars
addpath('QSSE_functions/analysis/')
addpath('QSSE_functions/evolution/')
addpath('QSSE_functions/others/')
addpath('QSSE_functions/unitaries/')


maxSchmidtrank=256;
Eps=1e-8;

Nqubitsvec=[2,8,20,50];
chivec=10.^linspace(-2,1,101);
Fid=zeros(length(Nqubitsvec),length(chivec));
Fid0=zeros(length(Nqubitsvec),length(chivec));
Fid1=zeros(length(Nqubitsvec),length(chivec));

gammar=1;

ground=[1;0];
sm=[0,1;0,0];
sp=sm';
sx=[0,1;1,0];
sy=[0,1j;-1j,0];
sz=[-1,0;0,1];
id=[1,0;0,1];
Pe=[0,0;0,1];
Pg=[1,0;0,0];
U_had=[1,-1;1,1]/sqrt(2);

for iterNqubits=1:length(Nqubitsvec)
    Nqubits=Nqubitsvec(iterNqubits)
    U_scat=zeros(2,2,2,2);
    U_scat(:,:,1,1)=id;
    U_scat(:,:,2,2)=[1,0;0,-1];
    MPO_scat=cell(1,Nqubits);
    
    
    if Nqubits==1
        MPO_scat{1}=tensor_contraction(tensor_contraction(U_scat,U_had(:,1),3,1),U_had,3,2);
    else
        MPO_scat{1}=tensor_contraction(tensor_contraction(U_scat,U_had(:,1),3,1),id,3,2);
        for itern=2:Nqubits-1
            MPO_scat{itern}=tensor_contraction(U_scat,id,4,2);
        end
        MPO_scat{Nqubits}=tensor_contraction(U_scat,U_had,4,2);
    end

    Gamma=cell(1,Nqubits);
    Lambda=cell(1,Nqubits);
    for itern=1:Nqubits
        Gamma{itern}=reshape(U_had*ground,[1,1,2]);
        Lambda{itern}=1;
    end

    [GammaZ,LambdaZ]=MPO_application(Gamma,Lambda,MPO_scat,1,Nqubits);
    
    GammaZ0=GammaZ;
    GammaZ0{Nqubits}=tensor_contraction(GammaZ0{Nqubits},[1,0],4,2);
    GammaZ1=GammaZ;
    GammaZ1{Nqubits}=tensor_contraction(GammaZ1{Nqubits},[0,1],4,2);
    
%     rho0=N_sites_density_matrix(GammaZ0,LambdaZ,1,Nqubits);
%     [U0,V0]=eig(rho0);
%     rho1=N_sites_density_matrix(GammaZ1,LambdaZ,1,Nqubits);
%     [U1,V1]=eig(rho1);
        
    for iterchi=1:length(chivec)

        chi=gammar+chivec(iterchi)/(Nqubits);

   %     t0q=-(gamma+1j*(delta-chi))/(gamma-1j*(delta-chi));
   %     t1q=-(gamma+1j*(delta))/(gamma-1j*(delta));
   
        t0q=-(gammar+2*1j*(-gammar/2))/(gammar-2*1j*(-gammar/2));
        t1q=-(gammar+2*1j*(chi-gammar/2))/(gammar-2*1j*(chi-gammar/2));

        U_scat=zeros(2,2,2,2);
        U_scat(:,:,1,1)=id;
        U_scat(:,:,2,2)=[t0q,0;0,t1q];

        MPO_scat=cell(1,Nqubits);
        
        U_lastphase=[1,0;0,exp(-1j*pi*Nqubits/2)];
        
        if Nqubits==1
            MPO_scat{1}=tensor_contraction(tensor_contraction(U_scat,U_had(:,1),3,1),U_had*U_lastphase,3,2);
        else
            MPO_scat{1}=tensor_contraction(tensor_contraction(U_scat,U_had(:,1),3,1),id,3,2);
            for itern=2:Nqubits-1
                MPO_scat{itern}=tensor_contraction(U_scat,id,4,2);
            end
            MPO_scat{Nqubits}=tensor_contraction(U_scat,U_had*U_lastphase,4,2);
        end

        % Scattering

        [Gammatemp,Lambdatemp]=MPO_application(Gamma,Lambda,MPO_scat,1,Nqubits);
        
        Gammatemp0=Gammatemp;
        Gammatemp0{Nqubits}=tensor_contraction(Gammatemp0{Nqubits},[1,0],4,2);
        Gammatemp1=Gammatemp;
        Gammatemp1{Nqubits}=tensor_contraction(Gammatemp1{Nqubits},[0,1],4,2);

       % Fid(iterNqubits,iterdelta)=abs(state_overlap(Gammatemp0,Lambdatemp,GammaZ0,LambdaZ)+state_overlap(Gammatemp1,Lambdatemp,GammaZ1,LambdaZ))^2;
        Fid0(iterNqubits,iterchi)=state_fidelity_v3(Gammatemp0,Lambdatemp,1,GammaZ0,LambdaZ,1)/sqrt(state_fidelity_v3(GammaZ0,LambdaZ,1,GammaZ0,LambdaZ,1));
        Fid1(iterNqubits,iterchi)=state_fidelity_v3(Gammatemp1,Lambdatemp,1,GammaZ1,LambdaZ,1)/sqrt(state_fidelity_v3(GammaZ1,LambdaZ,1,GammaZ1,LambdaZ,1));

    end
end

%% Plots

myBlue=[0.4,0.6,0.8];

figure
box on
for iterNqubits=1:length(Nqubitsvec)
    h=loglog(abs(chivec),1-Fid0(iterNqubits,:)-Fid1(iterNqubits,:));
    hold on
    h.LineWidth=3;
    h.Color=myBlue*(length(Nqubitsvec)-iterNqubits+1)/length(Nqubitsvec);
    set(gca,'FontSize',30)
    ax=gca;
    ax.YLim=[0,1.00];
    ax.TickLabelInterpreter='latex';
    ax.LineWidth=1;
end

h=loglog(abs(chivec),(chivec.^2)/4,'--r');
h.LineWidth=3;

pleg=legend;
pleg.Interpreter='latex';
for iterNqubits=1:length(Nqubitsvec)
    pleg.String{iterNqubits}=['$n_G =$ ',num2str(Nqubitsvec(iterNqubits))];
end
pleg.String{end}='$n_G^2(V-\gamma_r)^2/(2\gamma_r)^2$';

ax.XLabel.Interpreter='latex';
ax.XLabel.String='$n_G(V-\gamma_r)/\gamma_r$';
ax.YLabel.Interpreter='latex';
ax.YLabel.String='$1-\mathcal F_{\mathcal Z}$';
ax.YTick=[1e-4,1e-2,1e0];
ax.XTick=[1e-2,1e-1,1e0,1e1];
