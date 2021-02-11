function [Ybar,Cbar,Kbar,Lbar] = calc_ss(alpha, sigma, phi, Gbar, Abar, taubarL, taubarI, beta, delta)
syms Y K L C
f1 = Y == Abar*K^alpha * L^(1-alpha);
f2 = Y == C + delta*K +Gbar*Y;
f3 = L^phi * C^sigma == (1-taubarL)*Abar*(1-alpha)*K^(alpha)*L^(-alpha);
f4 = (1+taubarI) == beta *( Abar * alpha * K^(alpha - 1) * L^(1-alpha) +(1-delta)*(1+taubarI));

soln = solve([f1 f2 f3 f4], [Y K L C]);
% there could be multiple soln. we want only valid ones
validY = zeros(1,length(soln.Y));
validL = validY;
validC = validL;
for j=1:length(soln.Y)
    if (isreal(eval(soln.Y(j,1)))&& eval(soln.Y(j,1))>0)
        validY(j) = 1;
    end
    if (isreal(eval(soln.L(j,1)))&& eval(soln.L(j,1))>0)
        validL(j) = 1;
    end
    if (isreal(eval(soln.C(j,1)))&& eval(soln.C(j,1))>0)
        validC(j) = 1;
    end
end
ind = find((validY>0).*(validL>0) .*(validC>0)); % everything works in these spots
Ybar = eval(soln.Y(ind(1),1));
Cbar = eval(soln.C(ind(1),1));
Kbar = eval(soln.K(ind(1),1));
Lbar = eval(soln.L(ind(1),1));
end

