%% Test 1: 
% Solution exact: U_exact = [-7/23, 8/23]

A = (1/2).*[ 5, 0 ; 0, 3 ];
b = [ -3 ; 4 ];

eps=10^(-5);
kmax=10000;
rho1 = 0.1;
rho2 = 0.1;

C_in = [2,-1;-1,2];
d_in = [1; 1];
Lambda = [0;0];

C_eq = 0;
d_eq = 0;
Mu = [0;0];

[U,Lambda,Mu,k] = ArrowHurwicz(A,b,C_eq,d_eq,C_in,d_in,rho1,rho2,Mu,Lambda,eps,kmax);

%% Test 2:
%Solution exact: U_exact = [1, -1]

A = (1/2)*eye(2);
b = [ 1 ; -1 ];

eps=10^(-5);
kmax=10000;
rho1 = 0.1;
rho2 = 0.1;

C_in = [1,2;0,1];
d_in = [0; 0];
Lambda = [0;0];

C_eq = 0;
d_eq = 0;
Mu = [0;0];

[U,Lambda,Mu,k] = ArrowHurwicz(A,b,C_eq,d_eq,C_in,d_in,rho1,rho2,Mu,Lambda,eps,kmax);


%% Test 3:
%Solution exact pour N = 2: U_exact = [ 1 ; -1 ]
%Solution exact pour N = 3: U_exact = [ 1 ; -1 ; 0]
%Solution exact pour N = 4: U_exact = [ 1 ; -6/5 ; 3/5 ; -1]
%Solution exact pour N = 5: U_exact = [ 1 ; -1.2 ; 0.6 ; -1 ; 0 ]

N = 4;

A = (1/2)*eye(N);
b = ((-ones(1,N)).^(2:(N+1)))';

eps=10^(-5);
kmax=10000;
rho1 = 0.1;
rho2 = 0.1;

C_in = eye(N) + diag(2*ones(1,N-1),1);
d_in = zeros(N,1);
Lambda = zeros(N,1);

C_eq = 0;
d_eq = 0;
Mu = zeros(N,1);

[U,Lambda,Mu,k] = ArrowHurwicz(A,b,C_eq,d_eq,C_in,d_in,rho1,rho2,Mu,Lambda,eps,kmax);
