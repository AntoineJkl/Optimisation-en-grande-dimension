function [A,C] = CreationInstance(N,min_A,max_A,min_C,max_C)

%Creation de A
r_A=min_A+(max_A-min_A).*rand(N+1,2,2);

A=zeros(N+1,2,2);
for i=1:(N+1)
    deb=2*i-1;
    fin=2*i;
    A(i,:,:)=diag(r_A(deb:fin));
end

%Creation de C
r_C=min_C+(max_C-min_C).*rand(N,2);

C=zeros(N,2);
for i=1:(N)
    deb=2*i-1;
    fin=2*i;
    C(i,:)=[r_C(deb) ; r_C(deb)+r_C(fin)];
end

end

