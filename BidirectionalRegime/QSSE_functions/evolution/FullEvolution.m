function [Gamma,Lambda,Population,index1,indexN]...
        = FullEvolution(Gamma0,Lambda0,N_wav,N_sys,Ugamma,Ug,Ukappa, Nt, maxSchmidtrank, Eps, showprogress)
    
warning('on','mytest:maxrank')

count=0;

if showprogress
    prog_bar = waitbar(0);
    titleHandle = get(findobj(prog_bar,'Type','axes'),'Title');
    set(titleHandle,'FontSize',20);
    waitbar(0,prog_bar,sprintf('0.00%',0));
    tic
end

Gamma=Gamma0;
Lambda=Lambda0;
  
index1=0;
indexN=(N_sys-1)*(2*N_wav+1)+3*N_wav+1;

Population=zeros(Nt,indexN);

%% Time evolution  
for k=1:Nt 
    
    for iter=1:indexN-index1
        index=index1+iter;
        dim=size(Gamma{index},3);
        rho=single_site_reduced_state_Gamma_efficient(Gamma,Lambda,index);
        Population(k,iter)=rho(dim,dim);
    end
     
    % g interactions
    
    for j=1:N_wav
        [Gamma,Lambda] = LocalUnitaryEvolution(Gamma,Lambda,index1+(j-1)*3+1,3,Ug,maxSchmidtrank, Eps);
    end
    for j=1:N_wav
        [Gamma,Lambda]=Shift_site1_to_siteN(Gamma,Lambda,index1+(j-1)*3+2,index1+3*N_wav+1+2*(j-1),maxSchmidtrank, Eps);
        [Gamma,Lambda] = LocalUnitaryEvolution(Gamma,Lambda,index1+3*N_wav+1+2*(j-1),3,Ug,maxSchmidtrank, Eps);
        [Gamma,Lambda]=Shift_site1_to_siteN(Gamma,Lambda,index1+3*N_wav+1+2*(j-1),index1+(j-1)*3+2,maxSchmidtrank, Eps);
    end
    
    for iters=2:N_sys-1
        
        
        for j=1:N_wav
            
            [Gamma,Lambda]=Shift_site1_to_siteN(Gamma,Lambda,index1+3*N_wav+1+(iters-2)*(2*N_wav+1)+(j-1)*2+1,index1+3*N_wav+1+(iters-1)*(2*N_wav+1)+2*(j-1),maxSchmidtrank, Eps);
            [Gamma,Lambda] = LocalUnitaryEvolution(Gamma,Lambda,index1+3*N_wav+1+(iters-1)*(2*N_wav+1)+2*(j-1),3,Ug,maxSchmidtrank, Eps);
            [Gamma,Lambda]=Shift_site1_to_siteN(Gamma,Lambda,index1+3*N_wav+1+(iters-1)*(2*N_wav+1)+2*(j-1),index1+3*N_wav+1+(iters-2)*(2*N_wav+1)+(j-1)*2+1,maxSchmidtrank, Eps);
        end
       
    end
     
    % kappa interactions
    
    for j=N_wav-1:-1:1
        [Gamma,Lambda]=Shift_site1_to_siteN(Gamma,Lambda,index1+(j-1)*3+3,index1+3*N_wav+1-(N_wav-j+1),maxSchmidtrank, Eps);
    end
    [Gamma,Lambda] = LocalUnitaryEvolution(Gamma,Lambda,index1+3*N_wav+1-N_wav,N_wav+1,Ukappa{1},maxSchmidtrank, Eps);
    for j=1:N_wav-1
        [Gamma,Lambda]=Shift_site1_to_siteN(Gamma,Lambda,index1+3*N_wav+1-(N_wav-j+1),index1+(j-1)*3+3,maxSchmidtrank, Eps);
    end
    
    for iters=1:N_sys-1
        for j=N_wav-1:-1:1
            [Gamma,Lambda]=Shift_site1_to_siteN(Gamma,Lambda,index1+3*N_wav+1+(j-1)*2+2+(iters-1)*(2*N_wav+1),index1+3*N_wav+1+2*N_wav+1-(N_wav-j+1)+(iters-1)*(2*N_wav+1),maxSchmidtrank, Eps);
        end
        [Gamma,Lambda] = LocalUnitaryEvolution(Gamma,Lambda,index1+3*N_wav+1+2*N_wav+1-N_wav+(iters-1)*(2*N_wav+1),N_wav+1,Ukappa{1+iters},maxSchmidtrank, Eps);
        for j=1:N_wav-1
            [Gamma,Lambda]=Shift_site1_to_siteN(Gamma,Lambda,index1+3*N_wav+1+2*N_wav+1-(N_wav-j+1)+(iters-1)*(2*N_wav+1),index1+3*N_wav+1+(j-1)*2+2+(iters-1)*(2*N_wav+1),maxSchmidtrank, Eps);
        end
    end

    % gamma interactions
    
    [Gamma,Lambda]=Shift_site1_to_siteN(Gamma,Lambda,indexN+1,index1+(N_wav-1)*3+2,maxSchmidtrank, Eps);
    indexN=indexN+1;
    [Gamma,Lambda] = LocalUnitaryEvolution(Gamma,Lambda,index1+(N_wav-1)*3+1,2,Ugamma,maxSchmidtrank, Eps);
    [Gamma,Lambda]=Shift_site1_to_siteN(Gamma,Lambda,index1+(N_wav-1)*3+2,index1+1,maxSchmidtrank, Eps);
    index1=index1+1;
    
    [Gamma,Lambda]=Shift_site1_to_siteN(Gamma,Lambda,indexN+1,indexN-1,maxSchmidtrank, Eps);
    [Gamma,Lambda] = LocalUnitaryEvolution(Gamma,Lambda,indexN-2,2,Ugamma,maxSchmidtrank, Eps);
    [Gamma,Lambda]=Shift_site1_to_siteN(Gamma,Lambda,indexN-1,index1+1,maxSchmidtrank, Eps);
    indexN=indexN+1;
    index1=index1+1;
    
    
    
     % Display progression
     count=count+1;
     if showprogress && count>Nt/20
         waitbar(k/Nt,prog_bar,sprintf('%3.2f%%',100*k/Nt));
         count=0;
     end

end


if showprogress
    toc
    close(prog_bar)
end