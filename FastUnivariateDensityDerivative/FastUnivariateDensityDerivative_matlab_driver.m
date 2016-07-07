d=1;
N=10;
M=10;
R=1;
X=R*rand(d,N);
Y=R*rand(d,M);
h=0.01;
r=0;
epsilon=1e-6;

[D_fast]=FastUnivariateDensityDerivative(N,M,X,Y,h,r,epsilon)
disp('Fast Done');

[D_direct]=UnivariateDensityDerivative(N,M,X,Y,h,r)
disp('Direct Done');

Q=1/(sqrt(2*pi)*(power(h,(r+1))))

error=max(abs(D_direct-D_fast)/Q)

clear functions