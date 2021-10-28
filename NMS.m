function J = NMS(x, filtx11,force,time,Samp_freq, phi_0) %x=[Fom, Lom, A, MA]
%% defining calibration parameters

Fom =x(1);
Lom=x(2);
A=x(3);
MA=x(4);
gamma1=x(5);
gamma2=x(6);
%% params Sens analys
% Fom=12;
% Lom=6.3;
% A=-3;
% MA=1.5;
% gamma1=0.5;
%gamma2=0.5;
%% activation dynamics parameters: delay, gamma1, gamma2
delay = 2;
e=filtx11;
d=ceil(delay*Samp_freq/1000);

alpha=1+gamma1+gamma2;
beta1=gamma1+gamma2;
beta2=gamma1*gamma2;
u1=zeros(1,time);
for t=d+1:1:time
    u1(t)=alpha*e(t-d)-beta1*u1(t-1)-beta2*u1(t-2);
end

a=(exp(A.*u1)-1)/(exp(A)-1);

%% Contraction Dynamics parameters:
%Parameters for Flexor Carpi Radialis (FCR)
lst=24.4;
Lmt=lst+1.2*Lom; %cm 
%---------------
lm=ones(1,time);

lt=ones(1,time);
phi_t=zeros(1,time);
fal=zeros(1,time);
d1=0.56;
k=-1/(d1*d1);
aa=1.5;
bb=8;
cc=0.0866;
%------------
fpl=zeros(1,time);
fv=zeros(1,time);
vm=zeros(1,time);
Ft=zeros(1,time);
Fm=zeros(1,time);
Fmt=zeros(1,time);
Fl=zeros(1,time);
eps_t=zeros(1,time);
b=0.325;
a1=380*0.98;
vom=8*Lom;%(b*Fom)/a1;
theta=zeros(1,time);
lm(1:22)=1.001*Lom;
vm(20)=0;
%---------------
for i=d:1:time
    cur_a=a(i-1);
    cur_lm=lm(i-1);
    cur_lm=cur_lm/Lom;
    
    theta(i)=Lom*sin(phi_0)/cur_lm;
    
    phi_t(i)=asin(theta(i));
    lt(i)=Lmt-cur_lm*cos(phi_t(i));
    eps_t(i)=(lt(i)-lst)/lst;
    if(eps_t(i)<0)
        Ft(i)=0;
        
    elseif(eps_t(i)<0.0127 && eps_t(i)>0)
        Ft(i)=1480.3*eps_t(i)*eps_t(i);
        
    else
        Ft(i)=37.5*eps_t(i)-0.2375;
        
    end
    fal(i)=k*cur_lm*cur_lm-2*k*cur_lm+k+1;
    %-------------------------------------
    fpl(i)=exp(10*(cur_lm-1)-5);
    
    Fl(i)=fpl(i)+fal(i);

    fv(i)=aa/(1+exp(bb*vm(i)/vom-cc));
    
    Fm(i)=fal(i)*fv(i)*cur_a*Fom+fpl(i)*Fom;
    Fmt(i)=Fm(i)*cos(phi_t(i));
      
    vm(i)=(b*(Fl(i)-Fm(i)))/(Fm(i)+a1);

    lm(i)=lm(i-1)-trapz(vm(i-1:i))/Samp_freq;
  end
%% RMSE calculation
index_Force=force*40-100;
Mpred=Fmt*MA;
Mmes=index_Force*1.4;
Mmes=Mmes(1:time);
plot(abs(Mmes))
hold on
plot(abs(Mpred))
legend('Measured Moment','Estimated moment')
xlabel('Sample Number')
ylabel('Moment (N.cm)')
title('Measured vs. Estimated Moment')
error=Mmes-Mpred.';
J=sqrt(sum(error.^2)/length(error));

end

