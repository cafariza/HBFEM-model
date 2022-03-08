function Result=Newmark_Linear(M_bc,C_damp,K_bc,F_nm,fs)
%Input
% M:mass matrix (n*n) 
% C:damping matrix (n*n) 
% K:stiffness matrix (n*n) 
% f:external force matrix(n,N) 
% fs: sampling frequency 
% where n is the number of degrees of freedom, N is the length of data points of dynamic force
% Output: 
% Result: is a structure consist of 
% Result.Displacement: Displacement (n*N) 
% Result.Velocity: Velocity (n*N) 
% Result.Acceleration: Acceleration (n*N) 
%-----------------------------------------
[n, N]=size(F_nm);
dt=1/fs;   %Sampling rate
beta=1/4;  %Newmark parameters
gamma=1/2; %Newmark parameters
x(:,1)=zeros(n,1);                     %initial conditions
v(:,1)=zeros(n,1);                     %initial conditions
a(:,1)=zeros(n,1);
%a(:,1)=mldivide(M_bc,(F_nm(:,1)-C_damp*v(:,1)-K_bc*x(:,1)));%initial conditions

a1=1/(beta*dt^2)*M_bc+gamma/(beta*dt)*C_damp; %OK
a2=1/(beta*dt)*M_bc+(gamma/beta-1)*C_damp;
a3=(1/(2*beta)-1)*M_bc+dt*(gamma/(2*beta)-1)*C_damp;
K_hat=K_bc+a1; %OK
for i=1:1:N-1     
F_hat=F_nm(:,i+1)+a1*x(:,i)+a2*v(:,i)+a3*a(:,i);

x(:,i+1)=K_hat\F_hat;

v(:,i+1)=gamma/(beta*dt)*(x(:,i+1)-x(:,i))+(1-gamma/beta)*v(:,i)+dt*(1-gamma/(2*beta))*a(:,i);
a(:,i+1)=1/(beta*dt^2)*(x(:,i+1)-x(:,i))-1/(beta*dt)*v(:,i)-(1/(2*beta)-1)*a(:,i);
end
Result.Displacement=x;
Result.Velocity=v;
Result.Acceleration=a;
end