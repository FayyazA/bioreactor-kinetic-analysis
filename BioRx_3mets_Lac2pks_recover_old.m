function [E,Sfit_dyn]=BioRx_3mets_Lac2pks_recover_old(S_dyn,A,C,filename)
% Written by Renuka for 5mm Bioreactor on September 9 2013
%
% S_dyn is a matrix of size # metabolites(pyr, lac_in, Lac_ex)  by time
% A is the output of function BioRx_kinetics_Pyrfit
% A=[1/T1_Pyr, Kpl,FlowPyr 
%    alpha beta and constant values of gmma fit of input function]
% C is the output of function BioRx_Pyr_Lacin
% C=[1/T1Pyr, Kpl, Flow_pyr
%   1/T1Lin, Klp, K(MCT4)]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           DYNAMIC fit                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A = [0.020795808470777,0.047056897441457,0.177125718698039,3.841488640410969,17.960324782396610,5.175951486740151e+08]
%C = [0.020823757168569,0.015955048366388,0.208200244287854,0.045454545454523,2.731220277598615e-12,0.899999999999978]
%S_dyn = [0.01,2792.8,2058.75,19615.59,246385.42,811140.94,1809996.25,3394638,5304828,7336888,8955964,9412852,9370997,9465106,9877435,10494098,11081358,11347210,11290817,10931846,10517499,10306634,10139215,9892399,9494313,9053457,8624561,8073201,7623360,7263704,6936801,6611531,6253984,5873861,5479388.5,5094062,4736785,4409646,4092976.75,3759060,3388128.25,2989936.75,2586248,2209711,1871465.38,1572218.88,1323082.75,1110363.88,930999.38,782904.38,663012.38,563034.13,470105.66,397804.03,341377.19,289465.66,247925.2,206993.78,176032.78,149131,128371.73,108680.34,93985.05,78660.84,74291.62,60434.08,50603.66,41865.52,43139.45,31213.32,23687.73,17038.85,19287.12,15366.49,18372.98,11816.2,15437.3,16437.81,7108.88,8523.41,6441.93,6665.73,4451.49,7139.16,4600.94,3040.3,4058.45,5320.67,4047.55,3473.43,1801.65,3934.43,2773.34,533.03,1609.97,615.41;1.03,11.31,8644.16,0.23,9.03,2545.05,6995.86,13457.82,33614.95,49138.81,77454.23,89274.55,96225.16,97130.09,118560.2,108590.04,103446.23,107998.2,100842.95,97913.92,97363.42,98734.05,87443.72,72589.48,79298.07,72622.41,74245.41,64334.83,67011.81,57573.32,44654.24,46360.18,42693.35,47603.75,34515.66,41437.32,31445.29,32338.33,34172.66,25033.88,21333.98,25904.62,17404.02,14183.56,16717.28,14960.71,10271.08,10870.1,10624.51,8580.38,3793.83,10715.06,4940.39,5121.25,4556.65,713.99,1471.06,726.02,7412.34,2446.78,7224.25,5342.23,4045.21,12.38,3389.9,2514.92,427.21,2074.08,1491.07,1884.18,838.21,1123.89,344,0,1799.51,58.69,2432.65,2144.27,785.28,2300.79,1082.84,38.56,1894.85,1476.5,4044,22.26,64.08,0,2654.86,340.42,441.83,623.71,1714.4,21.64,3589.2,563.61;1871.62,193.12,203.75,301.17,40.56,1991.14,3131.16,11959.22,13928.36,35319.4,45855.14,68883.63,89189.8,110312.09,107267.59,127533.92,139614.86,137560.44,140044.98,139390.83,135064.41,121108.81,126250.37,121796.93,113820.52,108710.34,98564.91,91879.09,85379.74,80787.05,84898.23,75108.41,71605.59,56815.6,58009.42,47768.3,50514.14,43559.21,38924.68,35102,32201.36,26099.77,28005.01,30274.63,22950.76,23148.64,18140.15,17037.76,16864.68,21530,16713.38,9134.84,12613.81,11669.9,8966.81,10439.54,15067.85,10774.6,4144.96,7591.84,1556.85,5840.14,5388.17,7760.62,4452.08,463.51,6079.91,2629.23,2860.69,1323.17,3148.8,2221.56,413.46,5551.66,1805.11,74.3,3397.38,3399.56,37.07,142.76,26.33,17.78,730.6,1163.19,606.05,716.06,1585.59,1466.48,1261.09,2076.47,519.3,740.64,1366.83,2999.99,1982.95,377.35]

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

X0 = [C(1), C(2), C(3),... % 1/T1Pyr, Kpl, Flow_pyr
    C(4), C(5), 0.03,... %1/T1Lin, Klp, K(MCT4)
    1/36.7, C(3)*4, 0]; % 1/T1Lout,Flow_lacout, K(MCT1)
    
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
%%
function res_dyn = g_dyn(x)
        res_dyn = model_exchange_dyn(x) - S_dyn; 
        res_dyn = res_dyn(:);
end

% This should be optimized
opts = optimset('MaxIter',1000,'MaxFunEvals', 1e30,'TolX',1e-6,'TolFun', min(abs(S_dyn(:,end)))/1e12,'FinDiffType','central');

lb = [1/51, .0001, 1E-8,... % 1/T1Pyr, Kpl, Flow_pyr Flow_pyr used to be C(3)*0.9 1/T1P was C(1)-.001
      C(4)-.001, 0, 0.001,... %1/T1Lin, Klp, K(MCT4)
      1/37.7, 1E-8, 0]; % 1/T1Lout,Flow_lacout, K(MCT1)
 ub = [1/47, 0.08, 10.0,... % 1/T1Pyr, Kpl, Flow_pyr Flow_pyr used to be C(3)*1.1 1/T1P was C(1)-.001
      C(4)+.001, 0.1, 10,... %1/T1Lin, Klp, K(MCT4)
      1/35.7, 10.0, 0.9 ]; % 1/T1Lout,Flow_lacout, K(MCT1)
    
[X,resnorm,residual,exitflag,output,lambda,jacobian]  = lsqnonlin(@g_dyn, X0, lb, ub, opts);
Sfit_dyn = model_exchange_dyn(X);

exitflag
output
T1P=(1/C(1))
T1Lin=(1/X(4))
T1Lex=(1/X(7))
Kpl=(X(2))
Klp=X(5)
Flow_pyr=X(3)
Flow_lac=X(8)
K_MCT4=X(6)
K_MCT1=X(9)
Inp_func=[A(4) A(5) A(6)]
E=[1/C(1), X(2), X(3),... % 1/T1Pyr, Kpl, Flow_pyr
    1/X(4), X(5), X(6),... %1/T1Lin, Klp, K(MCT4)
    1/X(7), X(8), X(9),... % 1/T1Lout,Flow_lacout, K(MCT1)
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
E = [E(1),E(2),E(3),E(4),E(5),E(6),E(7),E(8),E(9),E(10),E(11),E(12),Rsqrd_pyr_total,Rsqrd_lac_in_magn,Rsqrd_lac_ex_magn,rmse_pyr,rmse_lac_in,rmse_lac_ex];
time=(1:size(Sfit_dyn,2))*TR;  %sec

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
print(gcf,'-dtiff','-r300',strcat(filename,"FINAL",".tif"))

end
    