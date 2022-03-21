% This program solves the <<Equivalent Beam>> problem based on the
% homogeneization method using Finite Element Method (FEM). 
% Example given in (Franco et al., 2022) 
%Units : mm, s, MN

%Do you need the following files and sub-routines to run this code:
%*mainvar.mat           Input data to fill the elementary matrices 
%                       ( Ki, Kg, K, lm, H, Ec, I, A, \Lambda etc..)
%*Matrices.mat          Elementary stiffness and mass matrices Ke qnd Me of
%                       the HBFEM model
%*amatrix.m             Assembly procedure function
%*Newmark_Linear        Function: Direct solution by Newmark procedure
%*EVN_dyn_Japan2_12000_EV1.csv Detailed model results used here for
%                       comparison purposes
%*m_elemu and k_elemu  functions within amatrix
%*IBRH161708020202.R2.sac.asc  Ground motion accelerogram in cm/s^2

clc, clear, close all
%% INPUT DATA AND EXTERNAL DATA FROM MF_LocalScale.m and KM_MatrixSolver.m

load('mainvar.mat') %Contents Ki, Kg, K, lm, \Lambda and other values obtained analytically.
k=K;
load('Matrices.mat')
%N_elem=input('How many elements do you want to use (Type 0 if you want to use the default number)?  ');
N_elem=0;
switch N_elem
    case 0
    N_elem=N; % Number of Finite Elements by default equal to number of stories
end
tic
%% Matrix assembly procedure 
Mele=M;
Kele=Kred;
ndg=3;
node=zeros(N_elem+1,2); % nodes
for i=1:N_elem+1
   node(i,1)=i; node(i,2)=H/N_elem*(i-1);
end
I=Kgb/Ew;
[K_bc, M_bc]=amatrix(Kele,Mele, node, ndg, H,Ki,Kgb,k, Lam, I);

%% EIGENVALUE PROBLEM : FUNDAMENTAL FRECUENCIES AND MODE SHAPES COMPUTATION
% --------------------------------------------------------------------------
% Calculate eigenvalues 
% --------------------------------------------------------------------------
%  
  ei=eig(K_bc,M_bc); % eigenvalues
% sorted natural angular frequencies [rad/s] 
  ef=sort(real(sqrt(ei))); 
% sorted natural frequencies [Hz]
  f_fem=ef/(2*pi);

 
% --------------------------------------------------------------------------
%  Mode shapes 
% --------------------------------------------------------------------------

[V,D]=eig(K_bc,M_bc); % eigenvalues and eigenvectors

[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs= zeros(size(V)+3);
Vs(4:end,1:end-3) = V(:,ind);

eigval = D(1,1);
eigvec = V(:,1);
error=K_bc*eigvec - eigval*M_bc*eigvec;

v1 = Vs(1:(ndg):end,1:end) ;    % Collecting only displacement degree's of Freedom
%
Vec = zeros(size(v1)) ;
for i = 1:size(v1,1)
    for j=1:size(v1,2)
        v(i,j)=v1(i,j)-v1(1,j);
    end
end

for i = 1:size(v1,2)
   ma= abs(max(v(:,i)));
   mi= abs(min(v(:,i)));
    Vec(:,i) = v(:,i)./(max(ma,mi)) ;
end

%Location of analysed nodes
x_f=H/L;x_0=0;
dx=H/(L*N_elem);    %Element length of the Discretized beam according to N_elem
x_efem=x_0:dx:x_f;       %Position vector of the Discretized beam according to N_elem
x_efem1=x_0:dx*L:x_f*L;       %Position vector of the Discretized beam according to N_elem
dx1=lw/L;           %Element length : cell size of the HM beam and detailed stucture
x_hm=[x_0:dx1:x_f];   %Position vector of the HM beam and detailed structure



%% TIME HISTORY ANALYSIS
% -----------------------------------------------------------------------
%  EARTHQUAKE HORIZONTAL LOAD
% -----------------------------------------------------------------------

grav=9810; %mm/s2
macc = load('IBRH161708020202.R2.sac.asc'); %units acceleration cm/s^2

% read acceleration history
TimAcc = macc(1:end,1);

% scale acceleration values and store them in Loading data structure
Loading.AccHst.Value = macc(1:end,2)*10/grav; 

% put corresponding time values in Loading data structure
Loading.AccHst.Time = TimAcc;
% time step increment
Deltat = TimAcc(3,1)-TimAcc(2,1);

fs=1/Deltat;                                               % Sampling Frequency

%*-----------------------------------------------------------------------
%*  RAYLEIGH DAMPING
%*-----------------------------------------------------------------------
omega1=omega(1); omega2=omega(2);
Xi=0.02;
a0=2*Xi*omega1*omega2/(omega1+omega2); a1=2*Xi/(omega1+omega2);
C_damp=a0*M_bc+a1*K_bc;


%*-----------------------------------------------------------------------
%*  Newmark  method
%*-----------------------------------------------------------------------

%Convert to equivalent nodal loads
[ns,ms]=size(M_bc);
F_nm=zeros(length(TimAcc),ns);
M_vec=M_bc*ones(ns,1);
F_nm=-M_vec*Loading.AccHst.Value'*1000;

tspan=[TimAcc(1,1) TimAcc(end,1)];
n_nm=tspan(2)/Deltat;


% Constant Average Newmark
Result3=Newmark_Linear(M_bc,C_damp,K_bc,F_nm,fs); 
toc


%% PLOTTING SECTION
% --------------------------------------------------------------------------
% Plot first three FEM mode shapes based on the displacement, U
% --------------------------------------------------------------------------

C = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]}; % Cell array of colors.
%Eigen analysis
figure
for i=1:3 
plot(Vec(1:end,i),x_efem,'color',C{i},'marker','o','LineWidth',2); hold on
end
legend({'Mode 1';'Mode 2';'Mode 3'})
xlabel ('u/u_{max}','FontName','Arial','FontSize',12);
ylabel('x/L','FontName','Arial','FontSize',12);
toc

% plot ground displacement history
figure  
plot(Loading.AccHst.Time,Loading.AccHst.Value,'-','Color','k','LineWidth',2);
set (gca,'FontSize',14,'FontName','Arial');
xlabel ('Time (sec)','FontName','Arial','FontSize',12);
ylabel('Ground Aceleration (g)','FontName','Arial','FontSize',12);

% plot TOP displacement history

figure 
plot(Loading.AccHst.Time,Result3.Displacement(28,:),'-','Color','k','LineWidth',2);
set (gca,'FontSize',12,'FontName','Arial');
xlabel ('Time (sec)','FontName','Arial','FontSize',12);
ylabel('Top roof displacement (mm)','FontName','Arial','FontSize',12);


%% Comparison with detailed structural model 
% --------------------------------------------------------------------------
% HBFEM model is compared with the detailed model
% --------------------------------------------------------------------------

% read acceleration history- Detailed model
TimDisp_num = load('EVN_dyn_Japan2_12000_EV1.csv');
Response.DispHst.Value = TimDisp_num(:,2);
Response.DispHst.Time = TimDisp_num(:,1);
Deltat_num = TimDisp_num(3,1)-TimDisp_num(2,1);

figure
plot(Loading.AccHst.Time,Result3.Displacement(28,:),'-',  'Color', [.2 .2 .2],'LineWidth',2);hold on
plot(Response.DispHst.Time,Response.DispHst.Value, '-',  'Color',[.7 .7 .7] ,'LineWidth',2);hold on

hold off;
set (gca,'FontSize',12,'FontName','Arial');
xlabel ('Time (sec)','FontName','Arial','FontSize',12);
ylabel('Top roof displacement (mm)','FontName','Arial','FontSize',12);
legend({'HBFEM';'Detailed FEM model'},  'location', 'NorthEast')

 grid on
 box on
 set(gca,'position',[0.09 0.10 0.88 0.88],'units','normalized')
 orient landscape

 
%*-----------------------------------------------------------------------
%*  FFT
%*-----------------------------------------------------------------------
%Ground Acceleration spectrum 
Fs = 1/Deltat;                                      % Sampling frequency                    
T = 1/Fs;                                       % Sampling period       
L = length(Loading.AccHst.Time);             % Length of signal
t = (0:L-1)*T;                                   % Time vector
 spec_dispAcc=fft((Loading.AccHst.Value));
 P2 = abs(spec_dispAcc/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f1 = Fs*(0:(L/2))/L;

figure
plot(f1,P1,'k') 
title('Single-Sided Amplitude Spectrum of X(t)-Ground motion')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%Response Acceleration spectrum: Detailed model
Dt=Deltat_num;
Fs = 1/Dt;                                      % Sampling frequency                    
T = 1/Fs;                                       % Sampling period       
L = length(Response.DispHst.Time);             % Length of signal
t = (0:L-1)*T;                                   % Time vector
spec_dispAcc=fft((Response.DispHst.Value));
P2 = abs(spec_dispAcc/L);
P3 = P2(1:L/2+1);
P3(2:end-1) = 2*P3(2:end-1);
f3 = Fs*(0:(L/2))/L;



%Response spectrum : HBFEM model 
Fs = 1/Deltat;                                          % Sampling frequency                    
T = 1/Fs;                                               % Sampling period       
L_g = length(Result3.Displacement(28,:));          % Length of signal
t = (0:L_g-1)*T;        % Time vector
 spec_disp=fft(Result3.Displacement(28,:));
 P2 = abs(spec_disp/L_g);
P4 = P2(1:L_g/2+1);
P4=P4';
P4(2:end-1) = 2*P4(2:end-1);
f4 = Fs*(0:(L_g/2))/L_g;



figure
plot(f4,P4, '-',  'Color', [.2 .2 .2]) ;hold on 
plot(f3,P3, '-.',  'Color', [.7 .7 .7])
xlabel('f (Hz)')
ylabel('|P(f)|')
legend({'HBFEM';'Detailed model'},  'location', 'NorthEast')
% specify axis dimensions for "ideal" plot
% set(gca,'GridLineStyle','--','GridColor',[0 0 0],'LineWidth',2)
grid on
box on
orient landscape
