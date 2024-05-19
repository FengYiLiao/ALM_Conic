close all;clc; clear;

%description

%parameter setting 
width  = 9;     % Width in inches
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
    subplot(1,3,idx);
    load(datapath+name{idx}+"_result"); 
    semilogy(Out.Dist(range),'-o');
    hold on 
    Feasibility = max([Out.Affinefeasi;Out.DAffinefeasi;Out.Conefeasi]);
    semilogy(Feasibility(range),'-o');
    semilogy(Out.DCostgap(range),'-o');
    xlabel('Iteration','interpreter','latex');
    xlim([range(1),range(end)]);
    set(gca,'TickLabelInterpreter','latex' ,'FontSize', fsz, 'LineWidth', alw); %<- Set properties
    if idx == 1
        set(gcf, 'Position', [300 100  width*100, height*100]); %<- Set size
        set(gca, 'Position', [0.1 0.3 0.225 0.6]); %<- Set properties
        title('Max-Cut','FontSize', fsz);
    elseif idx == 2
        set(gca, 'Position', [0.405 0.3 0.225 0.6]); %<- Set properties
        title('Matrix-Completion','FontSize', fsz);
    elseif idx == 3
        set(gca, 'Position', [0.71 0.3 0.225 0.6]); %<- Set properties
        title('BQP','FontSize', fsz);
    end
end
lg = legend( '$\frac{\mathrm{Dist}(X_k,\Omega_{\mathrm{P}})}{1+\|X^\star\|}$',...
    '$\left |\frac{\langle b,y_k \rangle -d^\star }{d^\star}\right|$', ...
    '$\max\left \{\frac{\|\mathcal{A}(X_k)-b\|}{1+\|b\|},\frac{\|X_k - \Pi_{\mathcal{K}}(X_k)\|}{1+\|X_k\|},\frac{\|C-\mathcal{A}^*(y_k)-Z_k\|}{1+\|C\|},\frac{\|Z_k - \Pi_{\mathcal{K}^*}(Z_k)\|}{1+\|Z_k\|} \right\}$',...
        'interpreter','latex','NumColumns',4,'Position',[0.115,0.03,0.8,0.1],'FontSize', fsz);

legend('boxoff');

%print("Numerical-result",'-depsc','-tiff');