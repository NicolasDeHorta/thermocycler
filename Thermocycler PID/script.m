close;clear;clc;

pkg load io
pkg load signal
pkg load control

datos=xlsread('termociclador.xlsx');
t=datos(501:end,1);
V=datos(501:end,2);
T=datos(501:end,3);
t=t.-5;
disp('_____________________')
#___________________________ORDEN 1 SIN DELAY__________________________________
function y=modelo1(x,t,T)
  k=x(1); tau=x(2);
  Du=1;
  
##  for i=1:length(t)
##    Tcalc(i)=T(1)+k*Du*(1-exp(-t(i)/tau));
##  endfor

Tcalc=T(1).+k*Du*(1.-exp(-t./tau));
  
  y=sum((Tcalc-T).^2);
endfunction

Du=1;                  %hallamos valores semilla
ko=(T(end)-T(1))/Du;
Tseed=0.632*ko*Du+T(1);
[i]=find(T<=Tseed);
tauo=t(i(1));

Xo1=[ko tauo];
[X1 fval1]=fminsearch(@(x) modelo1(x,t,T),Xo1);

disp('Ajuste Modelo Gv*Gp primer orden sin tiempo muerto')
k1=X1(1)
tau1=X1(2)
error_cuad1=fval1
disp('_____________________')


#___________________________ORDEN 2 SIN DELAY SUBAMORTIGUADO__________________________________


function y=modelo3(x,t,T)
  k=x(1); tau1=x(2); tau2=x(3);
  Du=1;
  dumping=(tau1+tau2)/(2*sqrt(tau1*tau2));
  tau=sqrt(tau1*tau2);
  beta=sqrt(1-dumping^2)/tau;
  phi=atan(sqrt(1-dumping^2)/dumping);
    Tcalc=T(1).+k.*Du*(1.-(1/sqrt(1-dumping^2)).*exp(-(dumping.*t./tau)).*sin(beta.*t.+phi));
  
  y=sum((Tcalc-T).^2);
endfunction


Xo3=[ko tauo tauo+0.01]; %aprox tau parecidos pero no iguales porque da error
[X3 fval3]=fminsearch(@(x) modelo3(x,t,T),Xo3);

disp('Ajuste Modelo Gv*Gp segundo orden sin tiempo muerto subamortiguado')
k3=X3(1)
tau13=X3(2)
tau23=X3(3)
error_cuad3=fval3
disp('_____________________')

#___________________________ORDEN 2 SIN DELAY AMORTIGUADO DINAMICA DENOM__________________________________

function [y,T]=modelo4(x,t,T)
  k=x(1); tau1=x(2); tau2=x(3); taun=x(4);
  Du=1;

    Tcalc=T(1).+k*Du*(1.+((taun-tau1)/(tau1-tau2))*exp(-t/tau1).+((taun-tau2)/(tau2-tau1))*exp(-t/tau2));
  
  y=sum((Tcalc-T).^2);
  
endfunction


Xo4=[ko tauo tauo+0.01 tauo/2]; %aprox tau parecidos pero no iguales porque da error
[X3 fval3]=fminsearch(@(x) modelo3(x,t,T),Xo4);

disp('Ajuste Modelo Gv*Gp segundo orden sin tiempo muerto subamortiguado')
k4=X3(1)
tau14=X3(2)
tau24=X3(3)
taun=X3(4)
error_cuad4=fval3
##
##Tcalc=T(1).+k4*Du*(1.+(taun-tau14)/(tau14-tau24)*exp(-t/tau14).+(taun-tau24)/(tau24-tau14)*exp(-t/tau24));
##plot(t,Tcalc)
disp('_____________________')


