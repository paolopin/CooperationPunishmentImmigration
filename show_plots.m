% ==============================================================
% subprogram show_results
% ==============================================================

a=1; b=.1; delta=.8; mu=.2;


% find the rightmost q, in a grid, that solves uc=ud for the model
% without heterogeneity and nu=0;
nu=0;
step_q=.01;
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
benchmark_q=q

choice_s=[.01 .05 .2 .5]
iterations=100;
step_nu=.001;

for j=1:4
    s=choice_s(j);
    for i=1:iterations
        nu=0+i*step_nu;
        [q,pi,pd]=equilibrium(a,b,delta,mu,nu,s);
        %[q,pi,pd]=nequilibrium(a,b,delta,mu,nu,s);
        Q(i)=q;
        PI(i)=pi;
        PD(i)=pd;
    end

    subplot(2,2,j)
    nu=step_nu:step_nu:iterations*step_nu;
    plot(nu,Q,'o')
    hold on
    plot(nu,PI,'og')
    plot(nu,PD,'or')
    plot([0],[benchmark_q],'kd')
    xlabel('\nu','FontSize', 16)
    legend('q','p_i','p_d','q_{no heterogeneity}','FontSize', 16 ,2) 
    title([ 's=',num2str(s) ],'FontSize', 16 ) 
    %title([ '\epsilon=',num2str(s) ],'FontSize', 16 ) 
end