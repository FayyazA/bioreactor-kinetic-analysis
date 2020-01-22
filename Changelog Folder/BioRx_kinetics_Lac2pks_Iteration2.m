function [E,Sfit_dyn,result_vector] = BioRx_kinetics_Lac2pks_Iteration2(S_dyn,A,C,filename)
%BIORX_KINETICS_LAC2PKS_ITERATION2 Summary of this function goes here
%   Detailed explanation goes here
% Written by Renuka for 5mm Bioreactor on September 9 2013
%
% S_dyn is a matrix of size # metabolites(pyr, lac_in, Lac_ex)  by time
% A is the output of function BioRx_kinetics_Pyrfit
% A=[1/T1_Pyr, Kpl,FlowPyr 
%    alpha beta and constant values of gmma fit of input function]
% C is the output of function BioRx_Pyr_Lacin
% C=[1/T1Pyr, Kpl, Flow_pyr
%   1/T1Lin, K(MCT4),Klp]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           DYNAMIC fit                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format shortG
TR = 3; %sec 
%a=[10,7.391,5.078e+09]; % initial values of gamma pdf fit for input function
F = 0.5; %flow of bioreactor 0.5ml/min
Inpfunc=zeros(3,length(S_dyn)); % Input function initialization
S_dyn=cell2mat(num2cell(S_dyn));
Nt = length(S_dyn);

%Intial conditions
% FROM ANALYSIS ON 2/3/15 FOR empty bioreactor data from 6/19/14
% T1Pyr=50.67 +- 0.49 s, T1Urea = 45.3 +- 0.27 sec,T1Lac=31.49 +- 0.28 s
% (from Wilson et al (2010) T1 @ 11.7T is Pyr: 48.3 s +- 1.2 s, Urea:43 +- 1 s)  
% (from Keshari & WIlson (2014) T1 @ 14.1 T for Lac: 32.7 s , DH @ 11.7t =34.8 � 0.1)  
% Flow=0.081 to 0.27 
% InputFunction- X(4)=2.4 - 3.7 (for Pyr), x(5)=17.3-21.4
%From UOK262 baseline dataset, X(4)=3.89+-0.4 , x(5)=15.7 +-1.6
%From empty beads - T1_Pyr = 48.61+-.4, T1_Lacex=36.7+-.05

X0 = [C(1), C(2), C(3),... % 1/T1Pyr, Kpl, Flow
    C(4), C(5),... %1/T1Lin, Klp, K(MCT4)
    1/36.7]; % 1/T1Lout
    
function Mest_dyn = model_exchange_dyn(x)
% Inital conditions
    Mest_dyn(:,1) = S_dyn(:,1);
    K  = [-x(1)-x(2)-x(3),0,0;... %Pyr=P(-R1P-kpl-Flow)+Lac_in(Klp)
          x(2),-x(4)-x(5),0;...   %Lac_in=P(kpl)+Lac_in(-R1PL-K(MCT4))
          0,x(5), -x(6)-x(3)] ;   %Lac_ex=Lac_in(K(MCT4))+Lac_out(-R1Lac_out-Flow)   
   
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
%%
function res_dyn = g_dyn(x)
        res_dyn = model_exchange_dyn(x) - S_dyn; 
        res_dyn = res_dyn(:);
end

% This should be optimized
opts = optimset('MaxIter',10000,'MaxFunEvals', 1e30,'TolX',1e-6,'TolFun', min(abs(S_dyn(:,end)))/1e14,'FinDiffType','central'); %Changed Tolfun term /1e12 to /1e14, did nothing

lb = [C(1)-.001, C(2)-.02,C(3)-0.02,...%.0001, C(3)*0.9,... % 1/T1Pyr, Kpl, Flow_pyr
      C(4)-0.001,C(5)-0.5,...%C(4)-.001, 0.001,... %1/T1Lin, K(MCT4)
      1/45];%1/45,]; 1/40,]; % 1/T1Lout,Flow_lacout, K(MCT1) Changed from 38.5 to 39
 ub = [C(1)+.001, C(2)+0.02,C(3)+0.02,...%0.08, C(3)*1.1,... % 1/T1Pyr, Kpl, Flow_pyr
      C(4)+0.001,C(5)+0.5,...%C(4)+.001, 2,... %1/T1Lin, K(MCT4)
      1/25];%1/25,]; %1/32,]; % 1/T1Lout,Flow_lacout, K(MCT1) Changed from 35.7 to 35
    
[X,resnorm,residual,exitflag,output,lambda,jacobian]  = lsqnonlin(@g_dyn, X0,lb,ub,opts);
Sfit_dyn = model_exchange_dyn(X);

exitflag
output
T1P=(1/X(1))
T1Lin=(1/X(4))
T1Lex=(1/X(6))
Kpl=(X(2))
Flow_pyr=X(3)
K_MCT4=X(5)
Inp_func=[A(4) A(5) A(6)]
E=[1/X(1), X(2), X(3),... % 1/T1Pyr, Kpl, Flow_pyr
    1/X(4), X(5),... %1/T1Lin, K(MCT4)
    1/X(6),... % 1/T1Lout,
    A(4), A(5), A(6)]; %Input Function

for i=1:3
    rmse(i)=gfit(S_dyn(i,:),Sfit_dyn(i,:),'4');     %NORMALIZED!
    R(i)=gfit(S_dyn(i,:),Sfit_dyn(i,:),'7');    
end

Rsqrd_pyr_total=R(1)
Rsqrd_lac_in_magn=R(2)
Rsqrd_lac_ex_magn=R(3)
rmse_pyr = rmse(1)
rmse_lac_in = rmse(2)
rmse_lac_ex = rmse(3)

time=(1:size(Sfit_dyn,2))*TR;  %sec

result_vector = [T1P,T1Lin,T1Lex,Kpl,Flow_pyr,K_MCT4,A(4),A(5),A(6),Rsqrd_pyr_total,Rsqrd_lac_in_magn,Rsqrd_lac_ex_magn,rmse_pyr,rmse_lac_in,rmse_lac_ex];
csvwrite(strcat(filename,".csv"),result_vector);
figure
subplot(3,1,1)
plot(time,(S_dyn(1,:)),'k*--',time,(Sfit_dyn(1,:)),'r-'); 
    legend('pyr magn data','pyr magn fit');
subplot(3,1,2)
plot(time,(S_dyn(2,:)),'k*--',time,(Sfit_dyn(2,:)),'r-'); 
    legend('lac_in magn data','lac(in) magn fit');
subplot(3,1,3)
plot(time,(S_dyn(3,:)),'k*--',time,(Sfit_dyn(3,:)),'r-'); 
    legend('lac_ex magn data','lac(ex) magn fit');
print(gcf,'-dtiff','-r300',strcat(filename,".tif"))



end

