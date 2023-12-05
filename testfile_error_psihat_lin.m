% Code file for Figure 3.2

clear, clc, close all
fprintf('Started %s\n', datestr(datetime('now')))

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

%% Setup

% Switch flag for saving results to txt-file
save_results = 0;

% Set bandwidth
a = 7;
N = 2^a;

% Set parameters to study
n = 2.^(0:15); % truncation parameter
lambda = [0.5;1;2]; % oversampling parameter
sigma = 1+lambda; % auxiliary parameter

% Initialization of error vectors
const = zeros(length(n),length(lambda)); 
err = zeros(length(n),length(lambda)); 

%% Error constant

% Computation of error constant for the linear frequency window function
for k = 1:length(lambda) 
    const(:,k) = real(sqrt(2*sigma(k)*N/3).*2*(1+lambda(k))./(pi^2*lambda(k)).*(n-sigma(k)*N).^(-3/2));
end%for

%% Reconstruction error

% Initialization of a fine grid for evaluation of reconstruction error
S = 1e5;
s = (-S:S)';
t = s/S;

% Initialization of vectors
Rm = zeros(length(t),1);
err_tilde = zeros(length(t),1);

% Loop for computation of the error
% Set function evaluations
f = @my_sinc;
ft = N^(1/2)*prod(f(N*pi,t),2);
% ft = (4*N/5)^(1/2)*prod((f(N*pi,t)+f(N*pi,(t-1))/2),2); % different test function

for i2 = 1:length(n) 
    % Set truncation parameter
    T = n(i2);
    j = (-T:T)'; % Corresponding index set
        
    for k = 1:length(sigma)
        % Set oversampling
        L = sigma(k)*N;

        % Set function evaluations
        fj = N^(1/2)*prod(f(N*pi,j/L),2);
%         fj = (4*N/5)^(1/2)*prod((f(N*pi,j/L)+f(N*pi,(j/L-1))/2),2); % different test function

        % Setup
        for i3 = 1:length(t)
            x = t(i3,:)-j/L;
            psi = (N+L)/2.*f((N+L)/2*pi,x).*f((L-N)/2*pi,x);

            % Evaluation of regularized WKS sums
            Rm(i3) = psi.'/L*fj;
        end%for

        % Computation of reconstruction errors
        err(i2,k) = norm(Rm-ft,inf);
    end%for
    
    fprintf(['T=',num2str(n(i2)),' done %s\n'], datestr(datetime('now')))
end%for

%% Visualization 

figure(1); loglog(n,const(:,1),'--',n,err(:,1),[192,192],[1e+03,1e-14],'-.',n,const(:,2),'--',n,err(:,2),[256,256],[1e+03,1e-14],'-.',n,const(:,3),'--',n,err(:,3),[384,384],[1e+03,1e-14],'-.'); 
xlabel('$T$'); legend('',['$\lambda=$ ',num2str(lambda(1))],'','',['$\lambda=$ ',num2str(lambda(2))],'','',['$\lambda=$ ',num2str(lambda(3))],''); title(['$N=$ ',num2str(N)])
title({'Figure 3.2: Maximum approximation error (3.13) (solid) and error constant (3.9) (dashed)', 'using the linear frequency window $\psi_{\mathrm{lin}}$ from (3.4) in (3.8) for the function', ' $f(t) = \sqrt{N} \,\mathrm{sinc}(N \pi t)$ with $N=128$, $T = 2^c$, $c \in \{0,\dots,15\}$, and $\lambda\in\{0.5,1,2\}$.'})
colororder(["#FF007F";"#FF007F";"#FF007F";"#D95319";"#D95319";"#D95319";"#008080";"#008080";"#008080"])
ylim([5e-13,2e+2])

%% Generate tables for tikz

if (save_results == 1)
fileID = fopen('error_psihat_lin.txt','w');
format = '%d %1.4e \n';
fprintf(fileID,'---------------------------------------------------------------\n');
fprintf(fileID,'---------------------------------------------------------------\n');
fprintf(fileID,['\n\n Reconstruction error for different lambda with N= ',num2str(N),'\n']);
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% Error constant
fprintf(fileID,'Error constant');
for k = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(k)),'\n\n']);
    matrix = [n.',const(:,k)];
    fprintf(fileID,format,matrix.');
end%for
% Reconstruction error 
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
fprintf(fileID,'Error');
for k = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(k)),'\n\n']);
    matrix = [n.',err(:,k)];
    fprintf(fileID,format,matrix.');
end%for
fclose(fileID);
end%if
fprintf('\n Finished %s\n', datestr(datetime('now')))

%% Function definitions

% Definition of the sinc function
function y = my_sinc(N,x)
    y = (sin(N*x)./(N*x));
    y(x==0) = 1; 
end