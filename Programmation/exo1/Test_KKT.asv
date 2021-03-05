function [OK] = Test_KKT(A,b,C_in,d_in,mu,C_eq,d_eq,lambda,u)
    %Contraintes d'inégalite:
    epsilon1 = 10^-2; %Tolerance
    condition1 = C_in*u - d_in <= epsilon1;
    bool1 = prod(condition1);

    %Multiplicateur positifs:
    epsilon2 = 10^-2; %Tolerance
    condition2 = mu >= -epsilon2;
    bool2 = prod(condition2);

    %Contraintes inégalite * multiplicateur null:
    epsilon3 = 10^-2; %Tolerance
    condition3 = abs(lambda.*(C_in*u - d_in)) <= epsilon3;
    bool3 = prod(condition3);

    %Contraintes egalite:
    epsilon4 = 10^-2; %Tolerance
    condition4 = abs(C_eq*u - d_eq) <= epsilon4;
    bool4 = prod(condition4);
    
    %Gradient lagrangien null:
    epsilon5 = 10^-2; %Tolerance
    condition5 = abs(A*u - b + C_in'*mu + C_eq'*lambda) <= epsilon5;
    bool5 = prod(condition5);

    OK = bool1 && bool2 && bool3 && bool4 && bool5;

end
