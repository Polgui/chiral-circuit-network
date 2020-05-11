function [Gamma,Lambda,output]=DoubleProjectiveMeasurement(Gamma,Lambda,index,indexq1,indexq2,input,maxSchmidtrank,Eps)
% index labels the ancilla, indexq1 the first qubit in the pair and
% indexq2 the second one

rho=single_site_reduced_state_Gamma_efficient(Gamma,Lambda,index);
projE=rho(2,2)/(rho(2,2)+rho(1,1));
r=rand;

[Gamma,Lambda]= sweep(Gamma,Lambda,1,length(Gamma),maxSchmidtrank, Eps);

if r>projE
    output=0;
    Gamma{index}=tensor_contraction(Gamma{index},[1,0;0,0],3,2);
     if input==0
         Gamma{indexq1}=tensor_contraction(Gamma{indexq1},[-1,0;0,1],3,2);
         Gamma{indexq2}=tensor_contraction(Gamma{indexq2},[-1,0;0,1],3,2);
     elseif input==1
         Gamma{indexq1}=tensor_contraction(Gamma{indexq1},1j*[-1,0;0,1],3,2);
     else
         error('wrong intput')
     end
elseif r<=projE
    output=1;
    Gamma{index}=tensor_contraction(Gamma{index},[0,1;0,0],3,2);
    if input==0
         Gamma{indexq2}=tensor_contraction(Gamma{indexq2},1j*[-1,0;0,1],3,2);
    end
else
    error('wrong output')
end   

[Gamma,Lambda]= sweep(Gamma,Lambda,1,length(Gamma),maxSchmidtrank, Eps);

end
