function [X,Sfit_dyn,results]=BioRx_kinetics_Pyrfit_starting_out_recover_old(S_dyn,filename)
% Written by Renuka for 5mm Bioreactor on September 9 2013
%
% S_dyn is a matrix of size # metabolites(pyr)  by time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           DYNAMIC fit                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%S_dyn = [31977.33,102686.37,208046.14,390548.41,707077.38,1113558.88,1588942,2064040.63,2501565,2891711.5,3154499.5,3310999,3390937.25,3479132,3530190,3501621.5,3501125,3543876.25,3597290.5,3654972.25,3610495.75,3496091.75,3325011.75,3057918.5,2840486,2659053.5,2520255,2387903,2267083.5,2135649.25,1983060.5,1829845.75,1656387.88,1487121.13,1315597.13,1148550.88,996494.13,855521.5,729917.31,621954.56,532742.19,455488.81,389883.69,333007.72,285029.28,247492.3,212415.23,184089.25,160941.84,139125.25,121596.55,107235.84,93839.16,81108.7,72274.75,65253.3,58635.03,51053.24,44795.58,40772.49,35621.19,31576.71,28247.61,25437.99,22846.2,20221.9,18298.46,17559.77,16087.91,14316.51,13643.99,13434.57,12192.03,11818.06,10932.75,10632.71,11025.09,9492.87,8494.15,6665.22,7011.01,5815.17,0,0,0,0,0,0,0,0,0]


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
X0 = [1/48,0.01,0.15,...    %1/T1P, Kpl,FlowPyr
    (2/skewness(S_dyn(1,:)))^2,19,5e8];  % alpha beta and constant values of gmma fit of input function

function Mest_dyn = model_exchange_dyn(x)                          %model_exchange_dyn(x) returns matrix Mest_dyn
% Inital conditions
    Mest_dyn(:,1) = S_dyn(:,1);                                    %First Column of Mest_dyn set to First Column of Pyruvate Concentration
    K  = [-x(1)-x(2)-x(3)] ; %Pyr=P(-R1P-kpl-Flow)                  Setting Rate Constant of dpyr = k * pyr
    Inpfunc(1) = [0];                                              %Sets first column element as 0   

    for k = 2:length(S_dyn)
       %gamma pdf approximation of input function , func = @(a,x)( a(3)*x.^(a(1)-1) .* exp(-x/a(2)) / (a(2)^a(1) * gamma(a(1))));
       Inpfunc(:,k) =( x(6)*(k*TR).^(x(4)-1) .* exp(-(k*TR)/x(5)) / (x(5)^x(4) * gamma(x(4))));
       
       flip=30 ; % low flip angle (degree)
        
       Mest_dyn(:,k) = (K^-1)*Inpfunc(:,k-1) + expm(K*TR) * (Mest_dyn(:,k-1)*cos(flip*pi/180) - (K^-1)*Inpfunc(:,k-1));%Mest_dyn(:,k) = expm(K*TR)*Mest_dyn(:,k-1)*cos(flip*pi/180) +Inpfunc(:,k-1);
    end
end
%%
function res_dyn = g_dyn(x)                                               %Makes function called g_dyn that outputs the residuals vector
        res_dyn = model_exchange_dyn(x) - S_dyn;                          %residuals are set to Input Function - Pyruvate Concentration
        res_dyn = res_dyn(:);                                             %Make a copy of residuals and reassign?
end

% This should be optimized MaxIter final version = 5000, use 2000 for
% starting point analysis
opts = optimset('MaxIter',2000,'MaxFunEvals', 1e32,'TolX',1e-9,'TolFun', min(abs(S_dyn(:,end)))/1e15,'FinDiffType','central','Display','iter');

lb = [1/51, .0001, 1E-8,...
      X0(4)-2, 8, 1e3];
 ub = [1/47,.1,10.0,...
      X0(4)+2, 22, 1e10];

[X,resnorm,residual,exitflag,output,lambda,jacobian]  = lsqnonlin(@g_dyn, X0, lb, ub, opts);
Sfit_dyn = model_exchange_dyn(X);

exitflag
output
T1=num2str(1/X(1))
Flow=X(3)
Kpl=X(2)
Inp_func=num2str(X(4:6))
A=X;
time=(1:size(Sfit_dyn,2))*TR;  %sec

rmse=gfit(S_dyn,Sfit_dyn,'4');     %NORMALIZED!
R=gfit(S_dyn,Sfit_dyn,'7');    
Rsqrd_pyr_total=R
rmse_pyr = rmse

results = [X(1),X(2),X(3),Rsqrd_pyr_total,rmse];

figure
plot(time,S_dyn,'k*--',time,Sfit_dyn,'r-'); 
    legend('pyr data','pyr fit');
print(gcf,'-dtiff','-r300',strcat(filename,"FINAL",".tif"))
end
    