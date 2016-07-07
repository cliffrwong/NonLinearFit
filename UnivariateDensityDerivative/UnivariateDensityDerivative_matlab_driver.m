d=1;
N=100;
M=100;
R=1;
X=R*rand(d,N);
Y=R*rand(d,M);
h=0.1;
r=6;
 
[D]=UnivariateDensityDerivative(N,M,X,Y,h,r);

clear functions