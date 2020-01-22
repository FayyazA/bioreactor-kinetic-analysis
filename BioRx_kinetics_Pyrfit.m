function [X,Sfit_dyn,result_vector]=BioRx_kinetics_Pyrfit(S_dyn,filename)
% Written by Renuka for 5mm Bioreactor on September 9 2013
%
% S_dyn is a matrix of size # metabolites(pyr)  by time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           DYNAMIC fit                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


format shortG
TR = 3; %sec 
F = 0.5; %flow of bioreactor 0.5ml/min
Inpfunc=zeros(1,length(S_dyn)); % Input function initialization of matrix with one row
S_dyn=cell2mat(num2cell(S_dyn));
Nt = length(S_dyn);

%Intial conditions
% FROM ANALYSIS ON 2/3/15 FOR empty bioreactor data from 6/19/14
% T1Pyr=50.67 +- 0.49 s, T1Urea = 45.3 +- 0.27 sec,T1Lac=31.49 +- 0.28 s
% (from Wilson et al (2010) T1 @ 11.7T is Pyr: 48.3 s +- 1.2 s, Urea:43 +- 1 s) (T1_pyr @14.1T =43.8)  
% (from Keshari & WIlson (2014) T1 @ 14.1 T for Lac: 32.7 s )  
% Flow=0.081 to 0.27 
% InputFunction- X(4)=2.4 - 3.7 (for Pyr), x(5)=17.3-21.4
%From UOK262 baseline dataset, X(4)=3.89+-0.4 , x(5)=15.7 +-1.6
%From empty beads - T1_Pyr = 48.61+-.4, T1_Lacex=36.7+-.05
X0 = [0,0.01,0.17,...    %1/T1P, Kpl,FlowPyr
    (2/skewness(S_dyn))^2,19,3e8];  % alpha beta a  nd constant values of gmma fit of input function
function Mest_dyn = model_exchange_dyn(x)                          %model_exchange_dyn(x) returns matrix Mest_dyn
% Inital conditions
    Mest_dyn(:,1) = S_dyn(:,1);                                    %First Column of Mest_dyn set to First Column of Pyruvate Concentration
    K  = [-1/50.67-x(2)-x(3)] ; %Pyr=P(-R1P-kpl-Flow)                  Setting Rate Constant of dpyr = k * pyr
    Inpfunc(1) = [0];                                              %Sets first column element as 0   

    for k = 2:length(S_dyn)
       %gamma pdf approximation of input function , func = @(a,x)( a(3)*x.^(a(1)-1) .* exp(-x/a(2)) / (a(2)^a(1) * gamma(a(1))));
       Inpfunc(:,k) =( x(6)*(k*TR).^(x(4)-1) .* exp(-(k*TR)/x(5)) / (x(5)^x(4) * gamma(x(4))));
       
       flip=30 ; % low flip angle (degree)
        
       Mest_dyn(:,k) = expm(K*TR)*Mest_dyn(:,k-1)*cos(flip*pi/180) +Inpfunc(:,k-1);
    end
end
%%
function res_dyn = g_dyn(x)                                               %Makes function called g_dyn that outputs the residuals vector
        res_dyn = model_exchange_dyn(x) - S_dyn;                          %residuals are set to Input Function - Pyruvate Concentration
        res_dyn = res_dyn(:);                                             %Make a copy of residuals and reassign?
end

% This should be optimized
opts = optimset('MaxIter',2000,'MaxFunEvals', 1e32,'TolX',1e-9,'TolFun', min(abs(S_dyn(:,end)))/1e18,'FinDiffType','central','Display','iter');

lb = [-10^-20,0.0001,0.01,...           %0.0001,...       %[1/50, .00001,0.005,...%.0001, .01,...
      X0(4)-2, 11, 1e3];
ub = [10^-20,0.5,0.2,...           %0.5,...           %[1/46,.1,0.3,...
      X0(4)+2, 22, 1e10];

[X,resnorm,residual,exitflag,output,lambda,jacobian]  = lsqnonlin(@g_dyn, X0,lb,ub, opts);
Sfit_dyn = model_exchange_dyn(X);

exitflag
output
Flow=X(3)
Kpl=X(2)
Inp_func=num2str(X(4:6))
A=X;
time=(1:size(Sfit_dyn,2))*TR;  %sec

rmse=gfit(S_dyn,Sfit_dyn,'4');     %NORMALIZED!
R=gfit(S_dyn,Sfit_dyn,'7');    
Rsqrd_pyr_total=R
rmse_pyr = rmse

result_vector = [Flow,Kpl,X(4),X(5),X(6),rmse_pyr,Rsqrd_pyr_total]
csvwrite(strcat(filename,".csv"),result_vector)
figure
plot(time,S_dyn,'k*--',time,Sfit_dyn,'r-'); 
    legend('pyr data','pyr fit');
    print(gcf,'-dtiff','-r300',strcat(filename,".tif"))
end
    