a1=0.9281;a2=0.09573;a3=-0.5053;a4=0.3718;
a5=-2.023;a6=0.7585;a7=0.6299;a8=0.5168;
a9=0.7168;a10=-0.9815;a11=0.2909;a12=-0.3141;
b1=-0.05578;b2=0.01112;b3=-0.1735;b4=1.675;
b5=0.2603;b6=1.106;b7=-1.073;b8=-1.519;b9=1.965;
c1=1.609;c2=0.1138;c3=0.6729;c4=2.357;c5=1.057;
c6=-4.444;c7=0.01727;c8=-0.01354;c9=-0.01224;c10=8.128;
c11=0.2007;c12=-1.594;
d1=0.006822;d2=1.014;d3=1.746;d4=2.378;d5=0.9434;
d6=-1.243;d7=35.79;d8=7.5;d9=62.94;d10=-6.93;
xi=50/7;eta=0.25;delta=0.27;lambda=0.81;psi=0.1297;phi=1.028;
fp=0.3;
A=a1*exp(a2/xi)*exp(a3*eta*xi*exp(a4/xi))*(1+a5*psi^a6+a7*psi^a8*(eta*xi)^a9)*(1+a10*lambda^a11*phi^a12);
B=xi^b1*eta^(b2*xi^b3/eta)*(b4*phi^b5*(b6*phi^b7)^lambda+b8*psi^b9);  
C=exp(c1/xi)*exp(c2*xi/((1+c3*xi)*eta))*(1+c4*psi^c5+c6*psi*exp(-eta)+c7*psi*exp(-xi))*(c8+c9*phi^c10+c11*lambda^c12*phi);
D=d1*xi^d2*eta^d3*(1+d4*psi^d5+d6*eta*xi*psi)*(1+d7*lambda^d8+d9*phi^d10);
delta_A=0.2262;
fdelta_A=(C*delta_A)^(B*delta_A^D);
epsilonTCrit=A*fdelta_A/(1+fdelta_A);
epsilonTCrit0=1.5*epsilonTCrit;
t0=15.9;t=20;
if (abs(t-t0)>0.00001)
    epsilonTCrit=(t0/t)^(0.8096*(1+1.503*psi^1.229))*epsilonTCrit;
end 
epsilonTCritTest1=epsilonTCrit;
epsilonTCritTest2=epsilonTCrit0+5*fp*(epsilonTCrit-epsilonTCrit0)/3;


