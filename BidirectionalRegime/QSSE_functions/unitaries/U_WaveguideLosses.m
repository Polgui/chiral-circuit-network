function U_loss=U_WaveguideLosses(dimB,theta)
% U_loss=U_WaveguideLosses(dimB,theta) 
% Returns the unitary tensor of dimension (dimB dimB dimB dimB)
% corresponding to the interaction between a waveguide time-bin and a loss
% channel time-bin. The loss is implemented as a beam-splitter coupling the
% two bins, with theta the beam-splitter angle. theta=0 corresponds to no
% loss, theta=pi/2 to 100% losses

a1=[0,1,0;0,0,0;0,0,0];
a2=[0,0,1;0,0,0;0,0,0];

H=KroneckerProduct(a1',a1)+KroneckerProduct(a1,a1')+KroneckerProduct(a2',a2)+KroneckerProduct(a2,a2');
U_loss=reshape(expm(1i*theta*H),[dimB,dimB,dimB,dimB]);

end

