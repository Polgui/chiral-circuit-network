function [Gamma,Lambda,ErrorTracker] = LocalUnitaryEvolution(Gamma,Lambda,site1,Nsites,Unitary,maxSchmidtrank, Eps)

%% Merge everything

temp=size(Unitary);
dimensions=temp(1:length(temp)/2);

Gamma_tot=Gamma{site1};
for l=1:Nsites-1
Gamma_tot=Lambda_multiplication(Gamma_tot,Lambda{site1+(l-1)},l+1);
Gamma_tot=tensor_contraction(Gamma_tot,Gamma{site1+l},l+1,1);
end

Gamma_tot=permute(Gamma_tot,[1 Nsites+1 2:Nsites Nsites+2]);

BondlengthL=size(Gamma_tot,1); % left MPS dimension of Gamma_feedback 
BondlengthR=size(Gamma_tot,2); % right MPS dimension of Gamma_feedback 

if site1 > 1
    Gamma_tot=Lambda_multiplication(Gamma_tot,Lambda{site1-1},1);
end

Gamma_tot=Lambda_multiplication(Gamma_tot,Lambda{site1+Nsites-1},2);


%% Apply the unitary

Gamma_tot=tensor_contraction(Gamma_tot,Unitary,(3:length(size(Gamma_tot))),((length(size(Unitary))/2 +1):length(size(Unitary))));
Gamma_tot=permute(Gamma_tot,[1 (3:length(size(Gamma_tot))) 2]);

%% Bring everything back in the original form
ErrorTracker=ones(1,2);
for l=1:Nsites-1
    dimB=dimensions(l);
    dimrest=prod(dimensions(l+1:end));
    Gamma_tot=reshape(Gamma_tot,[BondlengthL*dimB,dimrest*BondlengthR]); % reshape into a matrix

    % Apply the SVD
    [U,L1,V, projection_norm]=MySVD(Gamma_tot,maxSchmidtrank,Eps);
    ErrorTracker(1)=ErrorTracker(1)*projection_norm;
    ErrorTracker(2)=min(ErrorTracker(2),projection_norm);
    
    Bondlength=length(L1);
    Lambda{site1+(l-1)}=L1;
    % New state
    AL_feedback=permute(reshape(U,[BondlengthL,dimB,Bondlength]),[1,3,2]);
    
    if site1+(l-2)>0
        AL_feedback=Lambda_multiplication(AL_feedback,1./Lambda{site1+(l-2)},1);
    end
    Gamma{site1+(l-1)}=AL_feedback;
    % All the rest
    V=V(:,1:Bondlength);
    
    Gamma_tot=diag(L1)*V';
    BondlengthL=Bondlength;

end

% New state for the last one
Gamma{site1+Nsites-1}=Lambda_multiplication(permute(reshape(V',[BondlengthL,dimensions(end),BondlengthR]),[1,3,2]),1./Lambda{site1+Nsites-1},2);

