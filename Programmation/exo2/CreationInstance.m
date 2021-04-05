function [A,C] = CreationInstance(N,min_Ai,max_Ai,min_AN,max_AN,min_C,max_C)

%Creation de A
A=zeros(N+1,2,2);
for i=1:(N)
    r=min_Ai+(max_Ai-min_Ai).*rand(1,2);
    A(i,:,:)=diag(r);
end
r=min_AN+(max_AN-min_AN).*rand(1,2);
A(N+1,:,:)=diag(r);

%Creation de C
r_C=min_C+(max_C-min_C).*rand(N,2);

C=zeros(N,2);
for i=1:(N)
    deb=2*i-1;
    fin=2*i;
    C(i,:)=[r_C(deb) ; r_C(deb)+r_C(fin)];
end

end

