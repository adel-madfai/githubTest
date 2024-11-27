clc
clear all
close all
format long

%% par entré
alpha = 1;
N=5;
Nt=5;
L=5;
Tfinal=5;
Tleft=1;
Tright=2;
Tmid=2;

%% par dér
deltax=L/N;
deltat=Tfinal/Nt;

%% vect
x=linspace(0,L,N+1);
t=linspace(0,Tfinal,Nt+1);
T=Tmid*ones(N,1);

%% Coefficient du schéma
A=-((alpha*deltat)/(2*deltax*deltax));
C=-((alpha*deltat)/(2*deltax*deltax));
B=1+((alpha*deltat)/(deltax*deltax));

%% Cond limites
T(1)=Tleft;
T(N+1)=Tright;

%%
%matriceM:
M = diag(B * ones(N, 1)) + diag(A * ones(N-1, 1), -1) + diag(C * ones(N-1, 1), 1)
% Adjust the first and last rows for the boundary conditions
% (As seen in Equation (21) structure)

% First row includes the left boundary term (boundary conditions)
M(1, 1) = B; 
M(1, 2) = C;

% Last row includes the right boundary term (boundary conditions)
M(N-1, N-1) = A;
M(N-1, N-1) = B;
%% calcul de D:
for n=1:Nt
D=zeros(N,1);
    for j=2:N
    D(j-1)=T(j)+((alpha*deltat)/(2*deltax*deltax))*(T(j-1)-2*T(j)+T(j+1));
    end
    D(1)=D(1)-A*Tleft;
    D(N-1)=D(N-1)-C*Tright;
    
    T_final=M\D;
    
    %cond aux limites de T
    T(2:N)=T_final;

    T(1)=Tleft;
    T(N+1)=Tright;
    
    
end

T;

plot(x,T), xlabel('distance'), ylabel('Temperature')