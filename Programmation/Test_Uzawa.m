A=(1/2).*[ 5 0 ; 0 3];
b=[-3 ; 4];

eps=10^(-5);
kmax=10000;
pho=0.1;

d_eq=0;
C_eq=0;
Lambda=[0 ; 0];

d_in=[1 ; 1];
C_in=[ 2 , -1 ; -1 , 2 ];
Mu=[0 ; 0];

[U,Lambda,Mu,k]=Uzawa(A,b,d_eq,C_eq,C_in,d_in,pho,Mu,Lambda,eps,kmax);