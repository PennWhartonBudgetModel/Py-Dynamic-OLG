function Vtilde = valfuncold(x,kgrid,eff_inc,EV,pref_params)

% pref_params = [sigma4,chi3];
% tax_coefs = [avg_deduc; coefs; limit; X,mpci,rpci,tau_ss,v_ss_max];

cons1 = eff_inc - x;


% home-made linear interpolation
kpr = max(x(1),kgrid(1));
k_lower = find(kgrid(1:end-1)<=kpr,1,'last');
w = (kpr-kgrid(k_lower))/(kgrid(k_lower+1)-kgrid(k_lower));
V_star = (1-w)*EV(k_lower) + w*EV(k_lower+1);


if (cons1>0)&&(x>=0)
    Vtilde = (1/(1-pref_params(1))).*((cons1.^pref_params(2)).*(1^(1-pref_params(2)))).^(1-pref_params(1)) + V_star ; % - eta(t1+Tr-1).*WORK;
%     Vtilde = (1/(1-pref_params(1))).*((cons1.^pref_params(2)).*(1^(1-pref_params(2)))).^(1-pref_params(1)) + interp1(kgrid,EV,x,'linear','extrap') ; % - eta(t1+Tr-1).*WORK;
else
    Vtilde = -10e10*abs(x) - 10e10;
end

Vtilde = -1*Vtilde;

% 
% income = cap_inc + eff_wage*x(2);
% 
% fincome = (tax_coefs(6,1)/tax_coefs(5,1))*(cap_gain + eff_wage*x(2));%  + soc_sec*ge(t1,Tr);
% 
% deduc = tax_coefs(1,1) + tax_coefs(2,1).*fincome + tax_coefs(2,2).*fincome.^2 + tax_coefs(2,3).*fincome.^3 + tax_coefs(2,4).*fincome.^4;
% 
% ftax = (tax_coefs(5,1)/tax_coefs(6,1))*tax_coefs(3,1).*((max(fincome-deduc,0)) - ((max(fincome-deduc,0)).^(-tax_coefs(4,1)) + (tax_coefs(4,2))).^(-1/tax_coefs(4,1)));
% 
% sstax = tax_coefs(7,1).*min(eff_wage*x(2),tax_coefs(8,1));
% 
% cons1 = income - ftax - sstax - x(1) + beq;
% 
% 
% 
% if (x(1)>=0)&&(x(2)>=0)&&(x(2)<=1)&&(cons>0)
%     Vtilde = (1/(1-pref_params(1))).*((cons1.^pref_params(2)).*((1-x(2)).^(1-pref_params(2)))).^(1-pref_params(1)) + interp2(bgrid,kgrid,EV,'linear','extrap');
% else
%     Vtilde = 100000000*abs(x(2)*x(1));
% end
% 
% Vtilde = -1*Vtilde;
    


                        