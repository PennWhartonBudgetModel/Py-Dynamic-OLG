function [y, yprob, epsgrid] = cooper(ny,const,lambda,sigep)
%#codegen

% This function creates a Markov matrix that approximates an AR(1) process
% ny = gridsize
% const = constant on ar1 process
% lambda = coefficient on ar1 process
% sigep = standard deviation of error terms
% y(t) = mu*(1-lambda) + lambdat(t-1) + eps
% recovering mu

mu = const/(1-lambda);
sigy = sigep/(sqrt(1-lambda^2));
epsgrid = zeros(ny+1,1);    % creates grid for endpoints

for i1 = 1:ny+1    
    epsgrid(i1) = sigy*norminv((i1-1)/ny,0,1) + mu;
end

y = zeros(ny,1);    % grid of conditional means
for i1 = 1:ny    
    e1 = (epsgrid(i1)-mu)/sigy;
    e2 = (epsgrid(i1+1)-mu)/sigy;
    y(i1) = ny*sigy*(normpdf(e1,0,1) - normpdf(e2,0,1)) + mu;
end

yprob = zeros(ny,ny);
for i1 = 1:ny
    for i2 = 1:ny
        ei1 = epsgrid(i1);
        ei2 = epsgrid(i1+1);
        ej1 = epsgrid(i2);
        ej2 = epsgrid(i2+1);
        coop2([],ej1,ej2,sigy,sigep,mu,lambda);
        yprob(i1,i2) = (ny/(sqrt(2*pi*(sigy^2))))*quadgk(@coop2,ei1,ei2);        
    end    
end

end

function c2 =  coop2(x,ej1_,ej2_,sigy_,sigep_,mu_,lambda_)

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent initialized
persistent ej1
persistent ej2
persistent sigy
persistent sigep
persistent mu
persistent lambda

% Initialize parameters
if isempty(initialized)
    ej1    = 0;
    ej2    = 0;
    sigy   = 0;
    sigep  = 0;
    mu     = 0;
    lambda = 0;
    
    initialized = true;
    
end

% Set parameters if provided
if (nargin > 1)
    ej1    = ej1_;
    ej2    = ej2_;
    sigy   = sigy_;
    sigep  = sigep_;
    mu     = mu_;
    lambda = lambda_;
    
    c2 = [];
    return
end


nx = length(x);

c2 = zeros(1,nx);

for i1 = 1:nx

    term1 = exp(-((x(i1)-mu)^2)/(2*sigy^2));
    term2 = normcdf((ej2-mu*(1-lambda)-lambda*x(i1))/sigep,0,1);
    term3 = normcdf((ej1-mu*(1-lambda)-lambda*x(i1))/sigep,0,1);
    c2(i1) = term1*(term2-term3);
    
end

end

