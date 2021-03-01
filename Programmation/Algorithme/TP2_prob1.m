%% Première question

T=10;
N=5000;
x=0.5;
YN=[N*x];
for k = 1:N*T+1
    p=YN(k)/N;
    YN=[YN, binornd(N,p)];
end

XN=YN/N;
dt=(0:N*T+1);
dtN=dt/N;

plot(dt,YN,'b')
xlabel('Temps t')
ylabel('Processus Y^N_k')
title('Processus Y^N_k en fonction du temps')
xlim([0,N*T])
ylim([0,N])

figure()
plot(dtN,XN,'r')
xlabel('Temps t')
ylabel('Processus X^N_k')
title('Processus X^N_k en fonction du temps')
ylim([0,1])
xlim([0,T])

%% Comparaisons des x0
T=10;
N=5000;
% On fait varier x
i=1;
YN=zeros(5,N*T);

for x = [0.1,0.2,0.3,0.4,0.5]
    YN(i,1)=[N*x];
    for k = 1:N*T+1
        p=YN(i,k)/N;
        YN(i,k+1)=binornd(N,p);
    end
    i=i+1;
end

XN=YN/N;
dt=(0:N*T+1);
dtN=dt/N;

figure()
plot(dtN,XN(1,:),'r',dtN,XN(2,:),'b',dtN,XN(3,:),'g',dtN,XN(4,:),'m',dtN,XN(5,:),'y')
xlabel('Temps t')
ylabel('Processus X^N_k')
title('Processus X^N_k en fonction du temps')
legend('x0=0.1','x0=0.2','x0=0.3','x0=0.4','x0=0.5')
ylim([0,1])
xlim([0,T]) %le T peut être remplacé par 5 pour mieux voir dans 90% du temps


%% Deuxième Question
T=10;
N=5000;
x=0.4;
h=2^-8;

n=floor(T/h);
dt = 0:h:T;

XN=[x];

%Calcul des incréments Browniens
B = randn(1,n)*sqrt(h);

%Formule d'Euler
for k = 1:n
    temp=XN(k)+((XN(k)>0)&(XN(k)<1))*sqrt(XN(k)*(1-XN(k)))*B(k);
    XN=[XN, temp*((temp>0)&(temp<1))+1*(temp>=1)];
end

plot(dt,XN,'r')
xlabel('Temps t')
ylabel('Processus X^N_k')
title('Processus X^N_k en fonction du temps')
ylim([0,1])
xlim([0,T])


%% Comparaison des x
T=10;
N=5000;
h=2^-8;
n=floor(T/h);
dt = 0:h:T;
B = randn(5,n)*sqrt(h);
i=1;
XN=zeros(5,n);

for x = [0.1,0.2,0.3,0.4,0.5]
    XN(i,1)=x;
    for k = 1:n
        temp=XN(i,k)+((XN(i,k)>0)&(XN(i,k)<1))*sqrt(XN(i,k)*(1-XN(i,k)))*B(i,k);
        XN(i,k+1)=temp*((temp>0)&(temp<1))+1*(temp>=1);
    end
    i=i+1;
end
plot(dt,XN(1,:),'r',dt,XN(2,:),'b',dt,XN(3,:),'g',dt,XN(4,:),'m',dt,XN(5,:),'y')
xlabel('Temps t')
ylabel('Processus X^N_k')
title('Processus X^N_k en fonction du temps')
legend('x0=0.1','x0=0.2','x0=0.3','x0=0.4','x0=0.5')
ylim([0,1])
xlim([0,5])



%% Troisième question
T=10;
N=5000;
x=0.4;
h=2^-8;
M=10^4;

n=floor(T/h);

tauN=zeros(M,1);

XN=ones(M,n)*x;
    
%Calcul des incréments Browniens
B = randn([M,n])*sqrt(h);

%Formule d'Euler
for k = 1:n
    temp=XN(:,k)+((XN(:,k)>0)&(XN(:,k)<1)).*sqrt(XN(:,k).*(1-XN(:,k))).*B(:,k);
    XN(:,k+1)=temp.*((temp>0)&(temp<1))+1.*(temp>=1);
    tauN=(tauN==0).*k.*(XN(:,k+1)==0 | XN(:,k+1)==1)+tauN;
end

mean(tauN)*h
-2*(x*log(x)+(1-x)*log(1-x))


%% On fait un joli graphique en fonction du x choisi (attention avec 10^5)
T=10;
N=5000;
h=2^-8;
M=10^4;
n=floor(T/h);

meansTau=[];
dx=(0:0.1:1);
for x = dx
    tauN=zeros(M,1);
    XN=ones(M,n)*x;
    B = randn([M,n])*sqrt(h);

    for k = 1:n
        temp=XN(:,k)+((XN(:,k)>0)&(XN(:,k)<1)).*sqrt(XN(:,k).*(1-XN(:,k))).*B(:,k);
        XN(:,k+1)=temp.*((temp>0)&(temp<1))+1.*(temp>=1);
        tauN=(tauN==0).*k.*(XN(:,k+1)==0 | XN(:,k+1)==1)+tauN;
    end

    meansTau=[meansTau,mean(tauN)*h];
end

plot(meansTau,dx,'rx')
hold on
dx2=(0:0.001:1);
plot(-2*(dx2.*log(dx2)+(1-dx2).*log(1-dx2)),dx2)
title('Temps moyen d atteinte pour M=10^4')
ylabel('x_0')


%% Quatrième question
T=10;
N=500; %Déjà un peu long, si on augmente M faudra passer à N=100
x=0.3;
M=10^3;

tauN=zeros(M,1);
YN=ones(M,N*T)*N*x;
for k = 1:N*T+1
    p=YN(:,k)/N;
    YN(:,k+1)=binornd(N,p);
    tauN=(tauN==0).*k.*(YN(:,k+1)==0 | YN(:,k+1)==N)/N+tauN;
end

mean(tauN)

%% Cinquième question
% On mesure avec tic/toc le temps d'execution des deux derniers programmes
T=10;
N=500;
x=0.3;
h=2^-9;
M=10^3;
n=floor(T/h);


%
t1=tic;
tauN=zeros(M,1);
XN=ones(M,n)*x;
B = randn([M,n])*sqrt(h);

for k = 1:n
    temp=XN(:,k)+((XN(:,k)>0)&(XN(:,k)<1)).*sqrt(XN(:,k).*(1-XN(:,k))).*B(:,k);
    XN(:,k+1)=temp.*((temp>0)&(temp<1))+1.*(temp>=1);
    tauN=(tauN==0).*k.*(XN(:,k+1)==0 | XN(:,k+1)==1)+tauN;
end
MHT1=mean(tauN)*h
tf1=toc(t1)


%
t2=tic;
tauN=zeros(M,1);
YN=ones(M,N*T)*N*x;
for k = 1:N*T+1
    p=YN(:,k)/N;
    YN(:,k+1)=binornd(N,p);
    tauN=(tauN==0).*k.*(YN(:,k+1)==0 | YN(:,k+1)==N)/N+tauN;
end
MHT2=mean(tauN)
tf2=toc(t2)


% Moins de 0.5s vs plus de 40s, y'a pas photo, pour un même M, Wright
% Fischer c'est caca, alors que l'estimation de diffusion c'est le bien

