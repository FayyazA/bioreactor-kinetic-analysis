function data_vec = generate_synthetic(K,A)
Inpfunc=zeros(3,100);
Mest_dyn = zeros(3,100);
%A = [0,0,0,2.854500877,21.32582161,491427684.1];
TR=3;
for k = 2:100
       %gamma pdf approximation of input function , func = @(a,x)( a(3)*x.^(a(1)-1) .* exp(-x/a(2)) / (a(2)^a(1) * gamma(a(1))));
       Inpfunc(:,k) =[ ( A(6)*(k*TR).^(A(4)-1) .* exp(-(k*TR)/A(5)) / (A(5)^A(4) * gamma(A(4))));...
                        0;...
                        0]; 
       
       flip=30 ; % low flip angle (degree)
        
       Mest_dyn(:,k) = expm(K*TR)*Mest_dyn(:,k-1)*cos(flip*pi/180) +Inpfunc(:,k-1);
end
data_vec = Mest_dyn
end
