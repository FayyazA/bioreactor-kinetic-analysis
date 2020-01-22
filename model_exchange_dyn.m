function Mest_dyn = model_exchange_dyn(x)
% Inital conditions
    Mest_dyn(:,1) = S_dyn(:,1);
    K  = [-C(1)-x(2)-x(3),x(5),0;... %Pyr=P(-R1P-kpl-Flow)+Lac_in(Klp) Replace flow by 0.5 and put KMCT1 in entry(1,1)
          x(2),-x(4)-x(5)-x(6),x(9);...   %Lac_in=P(kpl)+Lac_in(-R1PL-klp-K(MCT4))+Lac_out(K(MCT1))
          0,x(6), -x(7)-x(8)-x(9)] ;   %Lac_ex=Lac_in(K(MCT4))+Lac_out(-R1Lac_out-Flow-K(MCT1))   
   
      Inpfunc(:,1) = [0;0;0];

    for k = 2:length(S_dyn)
       %gamma pdf approximation of input function , func = @(a,x)( a(3)*x.^(a(1)-1) .* exp(-x/a(2)) / (a(2)^a(1) * gamma(a(1))));
       Inpfunc(:,k) =[ ( A(6)*(k*TR).^(A(4)-1) .* exp(-(k*TR)/A(5)) / (A(5)^A(4) * gamma(A(4))));...
                        0;...
                        0]; 
       
       flip=30 ; % low flip angle (degree)
        
       Mest_dyn(:,k) = expm(K*TR)*Mest_dyn(:,k-1)*cos(flip*pi/180) +Inpfunc(:,k-1);
    end
end
