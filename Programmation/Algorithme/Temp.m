%EXEMPLES D'UTILISATION DE fmincon

%% Test 1:
% Solution exact: U_exact = [-7/23, 8/23]

A = (1/2).*[ 5, 0 ; 0, 3 ];
b = [ -3 ; 4 ];

f = @(x)  x'*A*x - b'*x;


x0 = [0;0];
C_in = [2,-1;-1,2];
d_in = [1; 1];

C_eq = [];
d_eq = [];

[x_sol,optimal_val,~,~,multiplicateur] = fmincon(f,x0,C_in,d_in,C_eq,d_eq);
x_exact = [ -7/23; 8/23 ];
disp('Solution : ');
disp([x_sol,x_exact]);
disp('Multiplicateurs (inégalité): ');
disp(multiplicateur.ineqlin);

%% Test 2:
% Solution exact: U_exact = [-7/23, 8/23]

A = 2*eye(4);
b = [1;-2;3;-1];

f = @(x)  x'*A*x - b'*x;


x0 = [0;0;0;0];
C_in = [];
d_in = [];

C_eq = [1,1,-1,-1;1,-1,1,-1];
d_eq = [1; 5];

[x_sol,optimal_val,~,~,multiplicateur] = fmincon(f,x0,C_in,d_in,C_eq,d_eq);
x_exact = [ 3/2 ; -7/8 ; 9/8 ; -3/2];
disp('Solution : ');
disp([x_sol,x_exact]);
disp('Multiplicateurs (inégalité): ');
disp(multiplicateur.eqlin);

