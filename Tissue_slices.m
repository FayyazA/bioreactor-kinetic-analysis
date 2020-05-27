function [result_vector,Sfit_dyn]=Tissue_slices(S_dyn,A,filename)
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
%S_dyn = [31977.33,102686.37,208046.14,390548.41,707077.38,1113558.88,1588942,2064040.63,2501565,2891711.5,3154499.5,3310999,3390937.25,3479132,3530190,3501621.5,3501125,3543876.25,3597290.5,3654972.25,3610495.75,3496091.75,3325011.75,3057918.5,2840486,2659053.5,2520255,2387903,2267083.5,2135649.25,1983060.5,1829845.75,1656387.88,1487121.13,1315597.13,1148550.88,996494.13,855521.5,729917.31,621954.56,532742.19,455488.81,389883.69,333007.72,285029.28,247492.3,212415.23,184089.25,160941.84,139125.25,121596.55,107235.84,93839.16,81108.7,72274.75,65253.3,58635.03,51053.24,44795.58,40772.49,35621.19,31576.71,28247.61,25437.99,22846.2,20221.9,18298.46,17559.77,16087.91,14316.51,13643.99,13434.57,12192.03,11818.06,10932.75,10632.71,11025.09,9492.87,8494.15,6665.22,7011.01,5815.17,0,0,0,0,0,0,0,0,0;0,0,0,0,0,5512.3,7519,11608,14146,24747,26946,31632,37045,38626,36755,39879,31721,35703,28279,35078,33609,23965,44795,26795,24349,28262,23707,22501,15271,16449,19056,17342,12485,11600,11498,10425,15027,10530,7272,5473,7924,5273,4877,4831,3174,3828,2237.23,3166,3493,2436,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
%A = [0.021276595744659,0.099999999999978,0.299999999999978,3.575810365588349,17.620365805005417,2.100285239949584e+08]
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

X0 = [A(1),A(2),A(3),...    %1/T1Pyr, Kpl, Flow_pyr
    0, 1/20, A(3)*4]; %klp, 1/T1L, Flow_lac NOTE: OLD STARTING POINT FOR T1L WAS 1/20, WORKED BETTER
    
function Mest_dyn = model_exchange_dyn(x)
% Inital conditions
    Mest_dyn(:,1) = S_dyn(:,1);

   
       K  = [-x(1)-x(2)-x(3),x(4);... %Pyr=P(-R1P-kpl-Flow)+Lac_in(Klp)
          x(2),-x(4)-x(5)-x(6)];   %Lac=P(kpl)+Lac(-klp-1/T1L-Flow_Lac)   
      
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
opts = optimset('MaxIter',1000,'MaxFunEvals', 1e30,'TolX',1e-6,'TolFun', min(abs(S_dyn(:,end)))/1e12,'FinDiffType','central');

lb = [1/51,.0001, 10^-8,... %A(1)-.0001
      10^-20, 1/38,0];
 ub = [1/47, 0.08,25,... %A(1)+.0001
     .1, 1/10,1];


[X,resnorm,residual,exitflag,output,lambda,jacobian]  = lsqnonlin(@g_dyn, X0, lb, ub, opts);
Sfit_dyn = model_exchange_dyn(X);

exitflag
output
T1Pyr=1/X(1)
T1Lin=(1/X(5))
Kpl=(X(2))
Klp=X(4)
Flow_pyr=X(3)
Flow_lac=X(6)
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

result_vector = [X(1),X(2),X(3),X(4),X(5),X(6),R(1),R(2),rmse(1),rmse(2)];
figure
subplot(2,1,1)
plot(time,(S_dyn(1,:)),'k*--',time,(Sfit_dyn(1,:)),'r-'); 
    legend('pyr magn data','pyr magn fit');
subplot(2,1,2)
plot(time,(S_dyn(2,:)),'k*--',time,(Sfit_dyn(2,:)),'r-'); 
    legend('lac_in magn data','lac_in magn fit');
print(gcf,'-dtiff','-r300',strcat(filename,"FINAL",".tif"))
figure
plot(time, S_dyn(1,:)./10^2, 'b*--', time, Sfit_dyn(1,:)./10^2, 'b-', time, Inpfunc(1,:)./10^2, 'k-', time, S_dyn(2,:), 'g*--', time, Sfit_dyn(2,:), 'g-');
legend('Pyruvate Data/10^2','Pyruvate Fit/10^2', 'Input Function / 10^2', 'Lactate Data', 'Lactate Fit');
print(gcf,'-dtiff','-r300',strcat(filename,"FINAL_MANUSCRIPT",".tif"))
end
    