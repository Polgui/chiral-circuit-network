function rho = N_sites_density_matrix(Gamma,Lambda,site1,N)
if length(Gamma)<site1+N-1
    error 'Gamma too small'
end

if N<2
    error 'N too small'
end

tens=Lambda_multiplication(Gamma{site1},Lambda{site1},2);
if site1>1
    tens=Lambda_multiplication(tens,Lambda{site1-1},1);
end

tens=permute(tensor_contraction(tens,conj(tens),1,1),[1,3,2,4]);

for iters=1:N-2
   
    newtens=Lambda_multiplication(Gamma{site1+iters},Lambda{site1+iters},2);
    
    tens=tensor_contraction(conj(newtens),tens,1,2);
    tens=tensor_contraction(newtens,tens,1,3);
    
    tens=permute(tens,[1,3,2,4,5:2+2*(iters+1)]);
    
end

newtens=Lambda_multiplication(Gamma{site1+N-1},Lambda{site1+N-1},2);
    
tens=tensor_contraction(conj(newtens),tens,1,2);
tens=tensor_contraction(newtens,tens,[1 2],[3 1]);

tens=permute(tens,[1:2:2*N,2:2:2*N]);

rho=reshape(tens,2^N,2^N);

end

