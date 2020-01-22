function [R1fit KPXfit  Sfit Xfit] = fit3exchange_vfa_sameT1(S, TR, flips, R1est, KPXest)

Nt = size(S,2);

Nflips = size(flips,2)/Nt;

Mzscale = zeros(3,Nt); Sscale = zeros(3,Nt);
for t= 1:Nt
    Iflips = [1:Nflips] + (t-1)*Nflips;
    Mzscale(:,t) = prod(cos(flips(:,Iflips)),2);
    for n = 1:Nflips
        Sscale(:,t) = Sscale(:,t) + sin(flips(:,Iflips(n))) .* prod(cos(flips(:,Iflips(1:n-1))),2);
    end
end

Mz_m0 = S(:,1) ./ Sscale(:,1);

X0 = [R1est KPXest Mz_m0.']; % R1P R1L R1A KPL KPA

    function Sest = model_exchange(x)
        
        A = [-x(1)-x(4)-x(5) 0 0;
            +x(4) -x(1) 0;
            +x(5) 0 -x(1)];
        
        Mz_m(1:3,1) = x(6:8) ;
%        Mz_m(1:2,1) = Mz_m0;
        for It = 1:Nt*Nflips
             Mz_p(:,It) = Mz_m(:,It) .* cos(flips(:,It));
             Mxy(:,It) = Mz_m(:,It) .* sin(flips(:,It));
             Mz_m(:,It+1) = expm(A * TR/Nflips) * Mz_p(:,It);
        end
        Sest = squeeze(sum(reshape(Mxy,3,Nflips,Nt),2));
%         for It = 1:Nt
%             Mz_p(:,It) = Mz_m(:,It) .* Mzscale(:,It);
%             Sest(:,It) = Mz_m(:,It) .* Sscale(:,It);
%             Mz_m(:,It+1) = expm(A * TR) * Mz_p(:,It);
%             
%         end
    end

    function res = g(x)
        res = model_exchange(x) - S;
        res = res(:);
    end

% some refinement possible here
opts = optimset('Display', 'off');%, 'FinDiffRelStep', [1 1 1 .1 .1]*sqrt(eps)); 
% FinDiffRelStep - use for normalization? or to preferentially fit certain parameters
% ie scale step size especially  for MzM0.  Or else scale data
% FinDiffType
%'MaxFunEvals', 1e8, 'TolFun', min(abs(S(:,end)))/1e14,'TolX', 1e-14);

%X = fminsearch(@g, X0, opts);
lb = [1/65, 1/65, 1/65, 0 0 0 0 0];
ub = [1/10 1/10 1/10 1 1 Inf Inf Inf ];

Xfit = lsqnonlin(@g, X0, lb, ub, opts);
R1fit = [Xfit(1),Xfit(1),Xfit(1)]; KPXfit = Xfit(4:5);
Minitfit = Xfit(6:8);
Sfit = model_exchange(Xfit);

end