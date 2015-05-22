
% Matlab code for
% ``Cooperation in a Society with Differential Treatment of Immigrants''
% Paolo Pin, and Brian Rogers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program evaluates the stable equilibrium values of p_i, p_d and q for given
% parameters a, b, delta, mu, nu and s (where s stands for the support (a-s,a+s) of the
% uniform distribution Phi for individuals a's)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Since we know that there is always a solution with [q,p_i,p_d]=[0,0,0]
% and that there is an unstable equilibrium, we must take care to be around
% the stable equilibrium - that is why we first compute as q_0 the q that would
% result without heterogeneity, and then at every step t we compute q_{t+1}
% "near" q_t (with parameter 'anchor')
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The output of the program is a vector of 3 elements - The input is vector of 6 elements: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [q,pi,pd] = equilibrium(a,b,delta,mu,nu,s)


precision = 10^(-6);
anchor=.99;
step_q=.01;

old_q=0;
old_pi=0;
old_pd=0;

% find the rightmost q, in a grid, that solves uc=ud for the model
% without heterogeneity, for immigrants
for q=1-step_q:-step_q:0
    D1= (q-1)*delta*( delta *(mu-1)-mu)-1;
    D2= q + delta^2-1;
    D3= (q-1) * ( delta + mu - delta*mu)^2 + 1;

    uc= ( (b+1)*q/D1 - b*(q-2) + q^2 * (delta^2 * (mu-1) - delta*mu - 1 ) / (  D2*D3 ) - q*(delta+1) / ( D2 * (delta*mu-delta-1) ) ) / (delta-1);
    udi= (a+1) * q * ( ( 1 - ( delta + mu - delta*mu)^2 )/D3 + 1 ) / ( delta*(nu-1) + 1 );
    if uc>udi
        break
    end
end

control=0;

while control==0
    % compute uc, ud/(a+1) and udi/(a+1)
    D1= (q-1)*delta*( delta *(mu-1)-mu)-1;
    D2= q + delta^2-1;
    D3= (q-1) * ( delta + mu - delta*mu)^2 + 1;
    uc= ( (b+1)*q/D1 - b*(q-2) + q^2 * (delta^2 * (mu-1) - delta*mu - 1 ) / (  D2*D3 ) - q*(delta+1) / ( D2 * (delta*mu-delta-1) ) ) / (delta-1);
    ud= q * ( ( 1 - ( delta + mu - delta*mu)^2 )/D3 + 1 ) / ( -delta + 1 );
    udi= q * ( ( 1 - ( delta + mu - delta*mu)^2 )/D3 + 1 ) / ( delta*(nu-1) + 1 );
    
    % compute pd and pi
    ad=uc/ud-1;
    pd=normcdf(ad,a,s);
    
    ai=uc/udi-1;
    pi=normcdf(ai,a,s);


    % compute q
    q= anchor*q + (1-anchor)*(1 - (1-delta)*(1-mu)*(1-pi) / ( (1-mu+mu*pd)*(1-delta+delta*nu*pi) ) );
    
    % check precision
    if max( abs( [q,pi,pd]-[old_q,old_pi,old_pd] ) ) < precision
        control=1;
    else
        
        old_q=q;
        old_pi=pi;
        old_pd=pd; 
    end
end

return
