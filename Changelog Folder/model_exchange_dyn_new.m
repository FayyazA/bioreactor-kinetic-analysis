function Mest_dyn = model_exchange_dyn_new(S_pyr,S_lac_in,S_lac_out,R_T1_pyr,Kpl,F,R_T1_lacin,k_inex,R_T1_lacex,k,theta,mag)
% Inital conditions
    S_dyn = S_pyr';
    S_lacin = S_lac_in';
    S_lacex = S_lac_out';
    TR = 3;
    k_constant = k;
    Mest_dyn = zeros(3,length(S_dyn));
    Mest_dyn(:,1) = S_dyn(:,1);
    K  = [-R_T1_pyr-Kpl-F,0,0;... %Pyr=P(-R1P-kpl-Flow)+Lac_in(Klp)
          Kpl,-R_T1_lacin-0-k_inex,0;...   %Lac_in=P(kpl)+Lac_in(-R1PL-klp-K(MCT4))+Lac_out(K(MCT1))
          0,k_inex, -R_T1_lacex-0-F] ;   %Lac_ex=Lac_in(K(MCT4))+Lac_out(-R1Lac_out-Flow-K(MCT1))   
   
      Inpfunc(:,1) = [0;0;0];

    for k = 2:length(S_dyn)
       %gamma pdf approximation of input function , func = @(a,x)( a(3)*x.^(a(1)-1) .* exp(-x/a(2)) / (a(2)^a(1) * gamma(a(1))));
       Inpfunc(:,k) = [( mag*(k*TR).^(k_constant-1) .* exp(-(k*TR)/theta) / (theta^k_constant * gamma(k_constant)));...
           0;
           0]
       
       flip=30 ; % low flip angle (degree)
        
       Mest_dyn(:,k) = expm(K*TR)*Mest_dyn(:,k-1)*cos(flip*pi/180) +Inpfunc(:,k-1);
    end
time = 1:length(S_dyn);
plot(time,Mest_dyn(1,:),'r-',time,S_pyr,'k*--')
print(gcf,'-dtiff','-r300',strcat("data3Pyr_startingpoint",".tif"))
plot(time,Mest_dyn(2,:),'r-',time,S_lac_in,'k*--')
print(gcf,'-dtiff','-r300',strcat("data3Lacin_startingpoint",".tif"))
plot(time,Mest_dyn(3,:),'r-',time,S_lac_out,'k*--')
print(gcf,'-dtiff','-r300',strcat("data3Lacex_startingpoint",".tif"))
end