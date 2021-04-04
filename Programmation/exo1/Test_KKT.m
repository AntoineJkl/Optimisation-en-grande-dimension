function [OK,bool1,bool2,bool3,bool4,bool5] = Test_KKT(A,b,C_in,d_in,mu,C_eq,d_eq,lambda,u,tolerance)
    %lambda -> contraintes egalites
    %mu -> contraintes inegalites
    %Contraintes d'inégalite:
    condition1 = C_in*u - d_in <= tolerance;
    bool1 = prod(condition1);

    %Multiplicateur positifs:
    condition2 = mu >= -tolerance;
    bool2 = prod(condition2);

    %Contraintes inégalite * multiplicateur null:
    condition3 = abs(mu.*(C_in*u - d_in)) <= tolerance;
    bool3 = prod(condition3);

    %Contraintes egalite:
    condition4 = abs(C_eq*u - d_eq) <= tolerance;
    bool4 = prod(condition4);
    
    %Gradient lagrangien null:
    condition5 = abs(A*u - b + C_in'*mu + C_eq'*lambda) <= tolerance;
    bool5 = prod(condition5);

    OK = bool1 && bool2 && bool3 && bool4 && bool5;

end
