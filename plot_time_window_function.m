% Code file for Figure 4.1

%% Setup

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

% Switch flag for saving results to txt-file
save_results = 0;

%% Initialization of parameters

N = 80; % bandwidth
m = 5; % truncation parameter
lambda = 1; % oversampling parameter
L = (lambda+1)*N; % oversampling
S = 1e3; % number of points

%% Plot time window functions

% Evaluation of sinc function
R = 2*m/L; % range of points
x = (-R:2*R/S:R-2*R/S)'; % set evaluation points
f_sinc = my_sinc(L*pi,x); % sinc

% Evaluation of sinh-type regularized sinc function
beta = m*pi*lambda./(1+lambda); % set shape parameter
sh = sinh(beta*sqrt(1-(L*x/m).^2))/sinh(beta); sh = real(sh); % sinh-type window
f_sh = f_sinc.*sh; % sinh-type regularized sinc function

% Evaluation of the continuous Kaiser-Bessel regularized sinc function
beta = m*pi*lambda./(1+lambda); % set shape parameter
ckb = (besseli(0,beta*sqrt(1-(L*x/m).^2))-1)/(besseli(0,beta)-1); ckb = real(ckb); % cKB window
f_ckb = f_sinc.*ckb; % cKB regularized sinc function

% Comparison of the regularized sinc functions
figure(1); plot(x,f_sinc,'-k',x,f_sh,x,f_ckb); 
legend('$\mathrm{sinc}(L \pi x)$','$\xi_\mathrm{sinh}(x)$','$\xi_\mathrm{cKB}(x)$');
title('Comparison of the regularized $\mathrm{sinc}$ functions $\xi$');
xlabel('$t$'); xlim([-R,R]); xticks([-R,-R/2,0,R/2,R]); xticklabels({'-$\frac{2m}{L}$','$-\frac{m}{L}$','0','$\frac{m}{L}$','$\frac{2m}{L}$'})
yticks([0 0.5 1])

%% Plot the corresponding functions in frequency domain

% Evaluation of the characteristic function
k = (-L:2*L/S:L)'; % Set evaluation points
chi = zeros(length(k),1); chi(abs(k)<L/2) = 1;

% Evaluation of sinh-type window function
beta = m*pi*lambda./(1+lambda); % set shape parameter
phihat_sinh = zeros(length(k),1);
rad = @(u) sqrt(u.^2/L^2-(lambda)^2/((2+2*lambda)^2));
func = @(u) besselj(1,2*pi*m*rad(u))./rad(u);
for i = 1:length(k)
    phihat_sinh(i) = integral(func,k(i)-L/2,k(i)+L/2);
end%for
phihat_sinh = beta.*phihat_sinh./(2*L*sinh(beta)); %phihat_sh = phihat_sh/L;

% Evaluation of the continuous Kaiser-Bessel window function
beta = m*pi*lambda./(1+lambda); % set shape parameter
phihat_cKB = zeros(length(k),1);
subst = @(u) 2*pi*m/(beta*L)*u;
func = @(u) sinh(beta*sqrt(1-subst(u).^2))./(beta*sqrt(1-subst(u).^2))-my_sinc(beta,subst(u));
for i = 1:length(k)
    phihat_cKB(i) = integral(func,k(i)-L/2,k(i)+L/2);
end%for
phihat_cKB = 2*m./((besseli(0,beta)-1)*L).*phihat_cKB; %phihat_ckb = phihat_ckb/L;

% Comparison of the functions in frequency domain
figure(2); plot(k,chi,'k--',k,phihat_sinh,k,phihat_cKB); 
legend('$\frac 1L\mathbf 1_{[-\frac L2,\frac L2]}(v)$','$\hat\xi_{\mathrm{sinh}}(v)$','$\hat\xi_{\mathrm{cKB}}(v)$'); 
title('Comparison of the function $\hat\xi$ in frequency domain');
xlim([-L,L]); xticks([-L,-L/2,0,L/2,L]); xticklabels({'$-L$','$-\frac{L}{2}$','$0$','$\frac{L}{2}$','$L$'}); 
ylim([0,1]); yticks([0,1/(2),1]); yticklabels({'$0$','$\frac{1}{2L}$','$\frac 1L$'});

%% Generate the paper plots

% Visualization for the sinh-type window function
figure(3); subplot(1,2,1); plot(k,chi,'k--',k,phihat_sinh)
legend('$\frac 1L\,\mathbf 1_{[-\frac L2,\frac L2]}(v)$','$\hat\xi_{\mathrm{sinh}}(v)$','Location','south'); 
title('$\hat\xi_{\mathrm{sinh}}$');
xlabel('$v$'); xlim([-L,L]); xticks([-L,-L/2,0,L/2,L]); xticklabels({'$-L$','$-\frac{L}{2}$','$0$','$\frac{L}{2}$','$L$'});
ylim([0,1]); yticks([0,1/(2),1]); yticklabels({'$0$','$\frac{1}{2L}$','$\frac 1L$'});
subplot(1,2,2); p = plot(x,f_sinc,'-k',x,sh,'--',x,f_sh); p(2).Color = [0.75 0.75 0.75]; p(3).Color = [0.8500 0.3250 0.0980];
legend('$\mathrm{sinc}(L \pi x)$','$\varphi_{\sinh}(x)$','$\xi_\mathrm{sinh}(x)$');
title('window function (4.1)');
xlim([-R,R]); xticks([-R,-R/2,0,R/2,R]); xticklabels({'-$\frac{2m}{L}$','$-\frac{m}{L}$','0','$\frac{m}{L}$','$\frac{2m}{L}$'})
yticks([0 0.5 1])
sgtitle('The sinh-type regularized $\mathrm{sinc}$ function $\xi_{\sinh}(t) \colon = \mathrm{sinc}(L\pi t)\,\varphi_{\sinh}(t)$ and its Fourier transform $\hat\xi_{\sinh}$.');

% Visualization for the cKB window function
figure(4); subplot(1,2,1); plot(k,chi,'k--',k,phihat_cKB)
legend('$\frac 1L\,\mathbf 1_{[-\frac L2,\frac L2]}(v)$','$\hat\xi_{\mathrm{cKB}}(v)$','Location','south'); 
title('$\hat\xi_{\mathrm{cKB}}$');
xlabel('$v$'); xlim([-L,L]); xticks([-L,-L/2,0,L/2,L]); xticklabels({'$-L$','$-\frac{L}{2}$','$0$','$\frac{L}{2}$','$L$'});
ylim([0,1]); yticks([0,1/(2),1]); yticklabels({'$0$','$\frac{1}{2L}$','$\frac 1L$'}); 
subplot(1,2,2); p = plot(x,f_sinc,'-k',x,ckb,'--',x,f_ckb); p(2).Color = [0.75 0.75 0.75]; p(3).Color = [0.8500 0.3250 0.0980];
legend('$\mathrm{sinc}(L \pi x)$','$\varphi_{\mathrm{cKB}}(x)$','$\xi_{\mathrm{cKB}}(x)$');
title('window function (4.2)');
xlim([-R,R]); xticks([-R,-R/2,0,R/2,R]); xticklabels({'-$\frac{2m}{L}$','$-\frac{m}{L}$','0','$\frac{m}{L}$','$\frac{2m}{L}$'})
yticks([0 0.5 1])
sgtitle('Figure 4.1: The regularized $\mathrm{sinc}$ function $\xi_{\mathrm{cKB}}(t) \colon = \mathrm{sinc}(L\pi t)\,\varphi_{\mathrm{cKB}}(t)$ using the continuous Kaiser-Bessel window function (4.2) and its Fourier transform $\hat\xi_{\mathrm{cKB}}$.');

%% Generate tables for tikz

if (save_results == 1)
fileID = fopen('time_windows.txt','w');
format = '%d %1.4e \n';
fprintf(fileID,['xmin=',num2str(-L),', xmax=',num2str(L),', xtick = {',num2str(-L),',',num2str(-L/2),',0,',num2str(L/2),',',num2str(L),'}, \n\n']);
fprintf(fileID,'chi\n\n');
matrix = [k,chi];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n sinh-type window function\n\n');
matrix = [k,phihat_sinh];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n Continuous Kaiser-Bessel window function\n\n');
matrix = [k,phihat_cKB];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n -----------------------------------------------\n\n');
fprintf(fileID,['xmin=',num2str(-10*R),', xmax=',num2str(10*R),', xtick = {',num2str(-10*R),',',num2str(-10*R/2),',0,',num2str(10*R/2),',',num2str(10*R),'}, \n\n']);
fprintf(fileID,'sinc\n\n');
matrix = [10*x,f_sinc];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n sinh-type window function\n\n');
matrix = [10*x,sh];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n sinh-type regularized sinc function\n\n');
matrix = [10*x,f_sh];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n Continuous Kaiser-Bessel window function\n\n');
matrix = [10*x,ckb];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n cKB regularized sinc function\n\n');
matrix = [10*x,f_ckb];
fprintf(fileID,format,matrix.');
fclose(fileID);
end%if

%% Function definitions

% Definition of the sinc function
function y = my_sinc(N,x)
    y = (sin(N*x)./(N*x));
    y(x==0) = 1; 
end%function