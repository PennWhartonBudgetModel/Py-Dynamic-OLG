function Vtilde = valfunc(x,kgrid,bgrid,cap_inc,cap_gain,eff_wage,beq,EV,pref_params, avg_deduc, coefs, limit, X, mpci, rpci, tau_ss, v_ss_max,t1,i1,nb,nk)




% home-made linear interpolation
bpr = (1/t1).*(bgrid(i1)*(t1-1) + eff_wage*x(2));
bpr = max(bpr,bgrid(1));
kpr = max(x(1),kgrid(1));

% b_lower = find(bgrid(1:end-1)<=bpr,1,'last');
% k_lower = find(kgrid(1:end-1)<=kpr,1,'last');
% 
% V_lower = EV(k_lower,b_lower);
% V_higher = EV(k_lower+1,b_lower+1);

% m_b = (V_higher-V_lower)/(bgrid(b_lower+1)-bgrid(b_lower));
% m_k = (V_higher-V_lower)/(kgrid(k_lower+1)-kgrid(k_lower));

% V_star = V_lower + m_b*(bpr - bgrid(b_lower)) + m_k*(x(1) - kgrid(k_lower));

% alternative method
b_lower = find(bgrid(1:nb-1)<=bpr,1,'last');
k_lower = find(kgrid(1:nk-1)<=kpr,1,'last');

b1 = bgrid(b_lower);
b2 = bgrid(b_lower+1);

k1 = kgrid(k_lower);
k2 = kgrid(k_lower+1);

V11 = EV(k_lower,b_lower);
V12 = EV(k_lower,b_lower+1);
V21 = EV(k_lower+1,b_lower);
V22 = EV(k_lower+1,b_lower+1);

denom = (k2-k1)*(b2-b1);

Val1 = V11*(k2-kpr)*(b2-bpr);
Val2 = V21*(kpr-k1)*(b2-bpr);
Val3 = V12*(k2-kpr)*(bpr-b1);
Val4 = V22*(kpr-k1)*(bpr-b1);

V_star = (Val1+Val2+Val3+Val4)/denom;






% pref_params = [sigma4,chi3];
% tax_coefs = [avg_deduc; coefs; limit; X,mpci,rpci,tau_ss,v_ss_max];

income = cap_inc + eff_wage*x(2);

fincome = (rpci/mpci)*(cap_gain + eff_wage*x(2));%  + soc_sec*ge(t1,Tr);

deduc = avg_deduc + coefs(1).*fincome + coefs(2).*fincome.^2 + coefs(3).*fincome.^3 + coefs(4).*fincome.^4;

ftax = (mpci/rpci)*limit*((max(fincome-deduc,0)) - ((max(fincome-deduc,0)).^(-X(1)) + (X(2))).^(-1/X(1)));

sstax = tau_ss*min(eff_wage*x(2),v_ss_max);

cons1 = income - ftax - sstax - x(1) + beq;
% size(cap_inc)
% size(eff_wage)
% size(income)
% size(ftax)
% size(sstax)
% size(cons1)



if (x(1)>=kgrid(1))&&(x(2)>=0)&&(x(2)<=1)&&(cons1>0)
    Vtilde = (1/(1-pref_params(1))).*((cons1.^pref_params(2)).*((1-x(2)).^(1-pref_params(2)))).^(1-pref_params(1)) + V_star;
%     Vtilde = (1/(1-pref_params(1))).*((cons1.^pref_params(2)).*((1-x(2)).^(1-pref_params(2)))).^(1-pref_params(1)) + interp2(bgrid,kgrid,EV,bpr,x(1),'linear',-10e10);
else
    Vtilde = -100000000*abs(x(2)*x(1)) - 10e10;
end

Vtilde = -1*Vtilde;
    


                        