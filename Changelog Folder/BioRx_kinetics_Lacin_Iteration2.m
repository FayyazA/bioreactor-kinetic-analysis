function [X,Sfit_dyn,result_vector] = BioRx_kinetics_Lacin_Iteration2(S_dyn,A,filename)
%BIORX_KINETICS_LACIN ITERATION2 Summary of this function goes here
%   Detailed explanation goes here
% Written by Renuka for 5mm Bioreactor on December 18 2014
%
% S_dyn is a matrix of size # metabolites(pyr, lac_in)  by time
% A is the output of function BioRx_kinetics_Pyrfit
% A=[1/T1_Pyr, Kpl,FlowPyr 
%    alpha beta and constant values of gmma fit of input function]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           DYNAMIC fit                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format shortG
TR = 3; %sec 
Inpfunc=zeros(size(S_dyn)); % Input function initialization
S_dyn=cell2mat(num2cell(S_dyn));
Nt = length(S_dyn);

%Intial conditions
% FROM ANALYSIS ON 2/3/15 FOR empty bioreactor data from 6/19/14
% T1Pyr=50.67 +- 0.49 s, T1Urea = 45.3 +- 0.27 sec,T1Lac=31.49 +- 0.28 s
% (from Wilson et al (2010) T1 @ 11.7T is Pyr: 48.3 s +- 1.2 s, Urea:43 +- 1 s)  
% (from Keshari & WIlson (2014) T1 @ 14.1 T for Lac: 32.7 s )  
% Flow=0.081 to 0.27 
% InputFunction- X(4)=2.4 - 3.7 (for Pyr), x(5)=17.3-21.4
%From UOK262 baseline dataset, X(4)=3.89+-0.4 , x(5)=15.7 +-1.6
%From empty beads - T1_Pyr = 48.61+-.4, T1_Lacex=36.7+-.05

X0 = [A(1),A(2),A(3),...    %1/T1Pyr, Kpl, Flow
    1/25, .04]; %1/T1Lin, K(MCT4)
    
function Mest_dyn = model_exchange_dyn(x)
% Inital conditions
    Mest_dyn(:,1) = S_dyn(:,1);

   
       K  = [-x(1)-x(2)-x(3),0;... %Pyr=P(-R1P-kpl-Flow)
          x(2),-x(4)-x(5)];   %Lac_in=P(kpl)+Lac_in(-R1PL-K(MCT4))   
      
    Inpfunc(:,1) = [0;0];

    for k = 2:length(S_dyn)
       %gamma pdf approximation of input function , func = @(a,x)( a(3)*x.^(a(1)-1) .* exp(-x/a(2)) / (a(2)^a(1) * gamma(a(1))));
       Inpfunc(:,k) =[ ( A(6)*(k*TR).^(A(4)-1) .* exp(-(k*TR)/A(5)) / (A(5)^A(4) * gamma(A(4))));...
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
%opts = optimset('MaxIter',800,'MaxFunEvals', 1e30,'TolX',1e-12,'TolFun', min(abs(S_dyn(:,end)))/1e12,'FinDiffType','central');
opts = optimset('MaxIter',10000,'MaxFunEvals', 1e30,'TolX',1e-6,'TolFun', min(abs(S_dyn(:,end)))/1e12,'FinDiffType','central');

lb = [A(1)-.0001, 0.0001, .01,...                                      %0,A(3)-.05,...                    
      1/35, .01];%1/35, .01];
 ub = [A(1)+.0001, 0.08,.3,...                                                     %A(2)+.1,A(3)+.05,...                   
     1/10,2];%1/10, 2];


[X,resnorm,residual,exitflag,output,lambda,jacobian]  = lsqnonlin(@g_dyn, X0,lb,ub, opts);
Sfit_dyn = model_exchange_dyn(X);

exitflag
output
T1Pyr=1/X(1)
T1Lin=(1/X(4))
Kpl=(X(2))
Flow=X(3)
K_MCT4=X(5)
C=X;

time=(1:size(Sfit_dyn,2))*TR;  %sec

for i=1:2
    rmse(i)=gfit(S_dyn(i,:),Sfit_dyn(i,:),'4');     %NORMALIZED!
    R(i)=gfit(S_dyn(i,:),Sfit_dyn(i,:),'7');    
end

Rsqrd_pyr_total=R(1)
Rsqrd_lac_in_magn=R(2)
rmse_pyr = rmse(1)
rmse_lac_in = rmse(2)

result_vector = [T1Pyr,T1Lin,Kpl,Flow,K_MCT4,Rsqrd_pyr_total,Rsqrd_lac_in_magn,rmse_pyr,rmse_lac_in];
csvwrite(strcat(filename,".csv"),result_vector);
figure
subplot(2,1,1)
plot(time,(S_dyn(1,:)),'k*--',time,(Sfit_dyn(1,:)),'r-'); 
    legend('pyr magn data','pyr magn fit');
subplot(2,1,2)
plot(time,(S_dyn(2,:)),'k*--',time,(Sfit_dyn(2,:)),'r-'); 
    legend('lac_in magn data','lac_in magn fit');
print(gcf,'-dtiff','-r300',strcat(filename,".tif"))


end

