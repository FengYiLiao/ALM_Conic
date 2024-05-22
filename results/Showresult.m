close all;clc; clear;

%description

%parameter setting 
width  = 10;     % Width in inches
height = 3;    % Height in inches
alw    = 0.75;    % AxesLineWidth
fsz    = 12;      % Fontsize
lw     = 2;      % LineWidth
msz    = 8;       % MarkerSize

colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4940 0.1840 0.5560],...
          [0.9290 0.6940 0.1250],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],...
          [1 0 0],[0 0 1],[0 1 0]}; 


datapath = "";
name     = {'G1_n20','n20r5MC','BQP_22'};



range = 1:6;

for idx = 1:3
    load(datapath+name{idx}+"_result");
    subplot(1,3,1);
    semilogy(Out.Dist(range),'-o');
    hold on 
    %
    Feasibility = max([Out.Affinefeasi;Out.DAffinefeasi;Out.Conefeasi;abs(Out.PCost - Out.DCost)/(1+abs(Out.OptimalCost))]);
    subplot(1,3,2);
    semilogy(Feasibility(range),'-o');
    hold on 

    subplot(1,3,3);
    semilogy(Out.DCostgap(range),'-o');
    hold on 

    
end

set(gcf, 'Position', [300 100  width*100, height*100]); %<- Set size

subplot(1,3,1);
xlim([range(1),range(end)]);
xlabel('Iteration','interpreter','latex');
set(gca,'TickLabelInterpreter','latex' ,'FontSize', fsz, 'LineWidth', alw); %<- Set properties
ylabel('$\log_{10}\left(\frac{\mathrm{Dist}(X_k,\Omega_{\mathrm{P}})}{1+\|X^\star\|}\right)$','interpreter','latex','FontSize', 14);
set(gca, 'Position', [0.1 0.3 0.225 0.6]); %<- Set properties
title('Primal distance','interpreter','latex','FontSize', fsz);
hold on 

subplot(1,3,2);
xlim([range(1),range(end)]);
xlabel('Iteration','interpreter','latex');
set(gca,'TickLabelInterpreter','latex' ,'FontSize', fsz, 'LineWidth', alw); %<- Set properties
ylabel('$\log_{10}\left ( \frac{ | \langle b,y_k \rangle -d^\star | }{| d^\star |} \right )$','interpreter','latex','FontSize', 16);
set(gca, 'Position', [0.42 0.3 0.225 0.6]); %<- Set properties
title('Dual cost value gap','interpreter','latex','FontSize', fsz);
hold on 


subplot(1,3,3);
xlim([range(1),range(end)]);
xlabel('Iteration','interpreter','latex');
set(gca,'TickLabelInterpreter','latex' ,'FontSize', fsz, 'LineWidth', alw); %<- Set properties
ylabel('$\log_{10}\left ( e \right )$','interpreter','latex', 'LineWidth', alw,'FontSize', 14)
set(gca, 'Position', [0.73 0.3 0.225 0.6]); %<- Set properties
title('KKT Residuals','interpreter','latex','FontSize', fsz);
hold on 
    

% lg = legend( '$\frac{\mathrm{Dist}(X_k,\Omega_{\mathrm{P}})}{1+\|X^\star\|}$',...
%     '$\left |\frac{\langle b,y_k \rangle -d^\star }{d^\star}\right|$', ...
%     '$\max\left \{\frac{\|\mathcal{A}(X_k)-b\|}{1+\|b\|},\frac{\|X_k - \Pi_{\mathcal{K}}(X_k)\|}{1+\|X_k\|},\frac{\|C-\mathcal{A}^*(y_k)-Z_k\|}{1+\|C\|},\frac{\|Z_k - \Pi_{\mathcal{K}^*}(Z_k)\|}{1+\|Z_k\|} \right\}$',...
%         'interpreter','latex','NumColumns',4,'Position',[0.115,0.03,0.8,0.1],'FontSize', fsz);
lg = legend( 'Max-Cut',...
    'Matrix-Completion', ...
    'Binary quadratic program','interpreter','latex','NumColumns',4,'Position',[0.115,0.03,0.8,0.1],'FontSize', fsz);

legend('boxoff');

%print("Numerical-result",'-depsc','-tiff');