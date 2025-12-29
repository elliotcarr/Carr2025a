clc, close all, clear

% Store figures in subfolder
mkdir Figures
path_name = './Figures/';

% Test problems
tau = 0.001; problem = '1'; sigma_vals = [0,0.05]; te = 0.1; Ms = 1001;
% tau = 0.0007; problem = '2'; sigma_vals = [0,0.05]; te = 0.1; Ms = 1001;
% tau = 0.0004; problem = '3'; sigma_vals = [0,0.05];te = 0.1; Ms = 1001;
% tau = 0.0001; problem = '4'; sigma_vals = [0,0.05]; te = 0.1; Ms = 1001;
% tau = 0.0015; problem = '5'; sigma_vals = [0,0.05]; te = 0.1; Ms = 1001;
% tau = 0.003; problem = '6'; sigma_vals = [0,0.05]; te = 0.1; Ms = 1001;
% tau = 0.01; problem = '7'; sigma_vals = [0,0.05]; te = 0.4; Ms = 1001;
% tau = 0.03; problem = '8'; sigma_vals = [0,0.05]; te = 0.4; Ms = 1001;

L = 0.002; % length of sample
k = 222; % thermal conductivity
rhoc = 2700*896; % volumetric heat capacity 
Qinf = 7000; % total amount of heat absorbed
alpha = k/rhoc; % thermal diffusivity
beta = 0.001; % heat pulse parameter
T0 = 10; % initial uniform sample temperature
Nx = 1001; % number of spatial nodes
x = linspace(0,L,Nx)'; % spatial locations
qt = @(t) (Qinf/beta^2)*(tau + (1-tau/beta)*t).*exp(-t/beta); % equation (10) in paper when q(t) = Qinf*t*exp(-t/beta)/beta^2
IQ = 2*beta; % integral(1-Q(t)/Qinf,0,Inf)
Tinf = T0 + Qinf/(rhoc*L); % steady state temperature

AbsTol = 1e-10; RelTol = 1e-6; % integration tolerances (default values)
font_size = 24;
% plot_legend = true;
plot_legend = false;
rng default; rng(3)

%% Temperature fields (plots for Figure 2 and Figure 4 of paper)
figure(1)
if isequal(problem,'7') || isequal(problem,'8')
    tplot = [0.0001,0.0005,0.0012,0.003,0.006,0.01,0.1]*4;
else
    tplot = [0.0001,0.0005,0.0012,0.003,0.006,0.01,0.1];
end
colors = turbo(9); colors(1,:) = []; colors(4,:) = [];
Ta = zeros(Nx,length(tplot));
for j = 1:length(tplot)
    for i = 1:length(x)
        Ta(i,j) = T0 + integral(@(u) integrand(u,tplot(j),qt,alpha,tau,k,x(i),L),0,tplot(j),...
            'AbsTol',AbsTol,'RelTol',RelTol);
    end
    plot(x,(Ta(:,j)-T0)/(Tinf-T0),'LineWidth',2,'Color',colors(j,:))
    hold on
    xw = tplot(j)*sqrt(alpha/tau);
    if xw <= L
        plot(xw,0,'.','MarkerSize',36,'Color',colors(j,:));
    end
    drawnow
end
if isequal(problem,'5') || isequal(problem,'6')
    ylim([0,5])
    set(gca,'FontSize',font_size,'Xtick',[0,L/2,L],'XTickLabel',{'$0$','$L/2$','$L$'},'YTick',0:5,...
        'YTickLabel',{'$0$','$1$','$2$','$3$','$4$','$5$'},'TickLabelInterpreter','LaTeX')
elseif isequal(problem,'7')
    ylim([0,8])
    set(gca,'FontSize',font_size,'Xtick',[0,L/2,L],'XTickLabel',{'$0$','$L/2$','$L$'},'YTick',0:2:8,...
        'YTickLabel',{'$0$','$2$','$4$','$6$','$8$'},'TickLabelInterpreter','LaTeX') 
elseif isequal(problem,'8')
    ylim([0,14])
    set(gca,'FontSize',font_size,'Xtick',[0,L/2,L],'XTickLabel',{'$0$','$L/2$','$L$'},'YTick',0:4:12,...
        'YTickLabel',{'$0$','$4$','$8$','$12$'},'TickLabelInterpreter','LaTeX')  
else
    ylim([0,4])
    set(gca,'FontSize',font_size,'Xtick',[0,L/2,L],'XTickLabel',{'$0$','$L/2$','$L$'},'YTick',0:4,...
        'YTickLabel',{'$0$','$1$','$2$','$3$','$4$'},'TickLabelInterpreter','LaTeX')
end
xlabel('$x$','FontSize',font_size+4,'Interpreter','LaTeX');
ylabel('$\widehat{T}(x,t)$','FontSize',font_size+4,'Interpreter','LaTeX');
text(0.6,0.9,['$s_{p} = ',num2str(sqrt(alpha/tau)),'$'],'FontSize',font_size,'Units','Normalized','Interpreter','LaTeX')
text(0.6,0.8,['$t_{p} = ',num2str(L*sqrt(tau/alpha)),'$'],'FontSize',font_size,'Units','Normalized','Interpreter','LaTeX')
print(gcf,[path_name,'flux_temp_distribution',num2str(problem)],'-depsc2')

%% Legend (for Figure 2 and Figure 4 of paper)
if plot_legend
    figure(2)
    set(gcf,'Color','w')
    scsize = get(gcf,'Position');
    set(gcf,'Position',[scsize(1),scsize(2),scsize(3)*3.2,scsize(4)*1.2]);
    labels = cell(size(colors,1),1);
    for j = 1:size(colors,1)
        plot(NaN,NaN,'Color',colors(j,:),'LineWidth',20)
        hold on
        labels{j} = ['{ }$t=',num2str(tplot(j)),'${ }'];
    end
    hold on
    hold off
    set(gca,'FontSize',font_size)
    [h] = legend(labels,'Location','NorthWest');
    set(h,'TextColor','k','Interpreter','LaTeX','Orientation','horizontal','EdgeColor','k')
    h.ItemTokenSize = [40,20];
    h.Position = h.Position + [0.025 0.045 0 0] ;
    h.NumColumns = 7;
    axis off
    set(gca, 'visible', 'off');
    if isequal(problem,'7') || isequal(problem,'8')
        print(gcf,[path_name,'temp_field_legend2'],'-depsc2')
    else
        print(gcf,[path_name,'temp_field_legend'],'-depsc2')
    end
end

%% Synthetic data (perfect signal)
ts = linspace(0,te,Ms);
dt = te/(Ms-1);

tp_indx = find(ts>=L*sqrt(tau/alpha),1,'first');

xa = L;
Tperfect = zeros(1,Ms);
for j = 1:Ms
    fprintf('Generating synthetic data: %g%% complete\n',round(j/Ms*100,1));
    Tperfect(j) = T0 + integral(@(u) integrand(u,ts(j),qt,alpha,tau,k,xa,L),0,ts(j),...
        'AbsTol',AbsTol,'RelTol',RelTol);
end

cnt = 2;
for sigma = sigma_vals
    
    cnt = cnt + 1;
    Tdata = Tperfect + sigma*randn(1,length(ts));
    Tdata(1:tp_indx-1) = T0;

    %% Calculate thermal diffusivity and relaxation time (perfect signal)
    alpha_true = alpha;
    alpha_approx = (L^2/6)/(trapz(ts,(Tinf-Tdata)/(Tinf-T0)) - IQ);

    tp_true = L*sqrt(tau/alpha);
    indx = find(Tdata~=T0,1,'first');
    tp = (ts(indx-1)+ts(indx))/2;
    tp_approx = tp;
    tau_true = tau;
    tau_approx = alpha_approx*(tp_approx/L)^2;

    %% Solve Cattaneo heat transfer model with computed values of alpha and tau
    Tfit = zeros(1,Ms);
    for j = 1:Ms
        fprintf('Generating fitted curve: %g%% complete\n',round(j/Ms*100,1));
        Tfit(j) = T0 + integral(@(u) integrand(u,ts(j),qt,alpha_approx,tau_approx,k,xa,L),0,ts(j),...
            'AbsTol',AbsTol,'RelTol',RelTol);
    end

    %% Back surface temperature history (plots for Figure 3 and Figure 5 of paper)
    figure(cnt)
    plot(ts,(Tdata-T0)/(Tinf-T0),'LineWidth',2,'Color',colors(2,:))
    hold on
    plot(ts,(Tfit-T0)/(Tinf-T0),'k--','LineWidth',2)       
    if isequal(problem,'5') || isequal(problem,'6')
        ylim([0,2.2])
        set(gca,'FontSize',font_size,'Xtick',[0,te],'XTickLabel',{'$0$','$t_{e}$'},'YTick',[0,1,2],...
            'YTickLabel',{'$0$','$1$','$2$'},'TickLabelInterpreter','LaTeX')        
        pos_lab = 0.3;
    elseif isequal(problem,'7')
        ylim([0,6.6])
        set(gca,'FontSize',font_size,'Xtick',[0,te],'XTickLabel',{'$0$','$t_{e}$'},'YTick',[0,2,4,6],...
            'YTickLabel',{'$0$','$2$','$4$','$6$'},'TickLabelInterpreter','LaTeX')
        pos_lab = 0.9;
    elseif isequal(problem,'8')
        ylim([0,14*1.1])
        set(gca,'FontSize',font_size,'Xtick',[0,te],'XTickLabel',{'$0$','$t_{e}$'},'YTick',[0:5:20],...
            'YTickLabel',{'$0$','$5$','$10$','$15$','$20$'},'TickLabelInterpreter','LaTeX')
        pos_lab = 0.9;        
    else
        ylim([0,1.1])
        set(gca,'FontSize',font_size,'Xtick',[0,te],'XTickLabel',{'$0$','$t_{e}$'},'YTick',[0,1],...
            'YTickLabel',{'$0$','$1$'},'TickLabelInterpreter','LaTeX')
        pos_lab = 0.3;
    end
    xlim([0,te])
    xl = xlabel('$t$','FontSize',font_size+4,'Interpreter','LaTeX');
    yl = ylabel('$\widehat{T}(L,t)$','FontSize',font_size+4,'Interpreter','LaTeX');
    text(0.5,pos_lab-0.1,['$\widehat{\tau} = ',num2str(tau_approx/10^(floor(log10(abs(tau_approx)))),'%.4f'),'\times 10^{',...
        num2str(floor(log10(abs(tau_approx))),'%i'),'}$'],'FontSize',font_size,'Units','Normalized','Interpreter','LaTeX')
    text(0.5,pos_lab-0.2,['$\widehat{\alpha} = ',num2str(alpha_approx/10^(floor(log10(abs(alpha_approx)))),'%.4f'),'\times 10^{',...
        num2str(floor(log10(abs(alpha_approx))),'%i'),'}$'],'FontSize',font_size,'Units','Normalized','Interpreter','LaTeX')
    if sigma == 0
        text(0.5,pos_lab,'Noise-free data','FontSize',font_size,'Units','Normalized','Interpreter','LaTeX')
        print(gcf,[path_name,'perfect_rear_temp',num2str(problem)],'-depsc2')
    else
        text(0.5,pos_lab,'Noisy data','FontSize',font_size,'Units','Normalized','Interpreter','LaTeX')
        print(gcf,[path_name,'noisy_rear_temp',num2str(problem)],'-depsc2')
    end

end

%% Legend (for Figure 3 and Figure 5 of paper)
if plot_legend
    figure(5)
    set(gcf,'Color','w')
    scsize = get(gcf,'Position');
    set(gcf,'Position',[scsize(1),scsize(2),scsize(3)*3.2,scsize(4)*1.2]);
    labels = cell(2,1);
    plot(NaN,NaN,'Color',colors(2,:),'LineWidth',20); labels{1} = '{ }Synthetic Data{ }';
    hold on
    plot(NaN,NaN,'k--','LineWidth',20); labels{2} = '{ }Fitted Curve{ }';
    hold off
    set(gca,'FontSize',font_size)
    [h] = legend(labels,'Location','NorthWest');
    set(h,'TextColor','k','Interpreter','LaTeX','Orientation','horizontal','EdgeColor','k')
    h.ItemTokenSize = [40,20];
    h.Position = h.Position + [0.025 0.045 0 0] ;
    h.NumColumns = 7;
    axis off
    set(gca, 'visible', 'off');
    print(gcf,[path_name,'temp_history_legend'],'-depsc2')
end

return

%% Table 1
tau_vals = [0.001,0.0007,0.0004,0.0001];
Ms_vals = 100001;
tau_approx_vals = zeros(length(Ms_vals),length(tau_vals));
alpha_approx_vals = zeros(length(Ms_vals),length(tau_vals));
for jj = 1:length(tau_vals)

    jj

    tau = tau_vals(jj);

    for ii = 1:length(Ms_vals)

        Ms = Ms_vals(ii);

        ts = linspace(0,te,Ms);
        tp_indx = find(ts>=L*sqrt(tau/alpha),1,'first');
        qt = @(t) (Qinf/beta^2)*(tau + (1-tau/beta)*t).*exp(-t/beta);

        xa = L;
        Tperfect = zeros(1,Ms);
        for j = 1:Ms
            fprintf('%i: %g%% complete\n',jj,round(j/Ms*100,1));
            Tperfect(j) = T0 + integral(@(u) integrand(u,ts(j),qt,alpha,tau,k,xa,L),0,ts(j),...
                'AbsTol',AbsTol,'RelTol',RelTol);
        end

        Tdata = Tperfect;
        Tdata(1:tp_indx-1) = T0;

        % Calculate thermal diffusivity and relaxation time (perfect signal)
        alpha_true = alpha;
        alpha_approx_vals(ii,jj) = (L^2/6)/(trapz(ts,(Tinf-Tdata)/(Tinf-T0)) - IQ);

        tp_true = L*sqrt(tau/alpha);
        indx = find(Tdata~=T0,1,'first');
        tp = (ts(indx-1)+ts(indx))/2;
        tp_approx = tp;
        tau_true = tau;
        tau_approx_vals(ii,jj) = alpha_approx_vals(ii,jj)*(tp_approx/L)^2;
    end
end
for ii = 1:length(Ms_vals)
    Ms = Ms_vals(ii);
    dt = te/(Ms-1);
    fprintf('%i & \\num{%1.2e} & \\num{%1.4e} & \\num{%1.4e} & \\num{%1.4e} & \\num{%1.4e}\\\\\n',...
        Ms,dt,tau_approx_vals(ii,1),alpha_approx_vals(ii,1),...
        tau_approx_vals(ii,2),alpha_approx_vals(ii,2));
end
for ii = 1:length(Ms_vals)
    Ms = Ms_vals(ii);
    dt = te/(Ms-1);
    fprintf('%i & \\num{%1.2e} & \\num{%1.4e} & \\num{%1.4e} & \\num{%1.4e} & \\num{%1.4e}\\\\\n',...
        Ms,dt,tau_approx_vals(ii,3),alpha_approx_vals(ii,3),...
        tau_approx_vals(ii,4),alpha_approx_vals(ii,4));
end

%% Sub functions
function I = integrand(u,t,qt,alpha,tau,k,x,L)
% See equations (16) and (17) of manuscript
I = 0;
M1 = floor(1/(2*L)*(sqrt(alpha/tau)*t-x));
M2 = floor(1/(2*L)*(sqrt(alpha/tau)*t+x)-1);
for n = 0:M1
    I = I + heaviside(u-sqrt(tau/alpha)*(x+2*n*L)).*besseli(0,(1/(2*tau))*sqrt(u.^2-tau/alpha*(x+2*n*L)^2));
end
for n = 0:M2
    I = I + heaviside(u-sqrt(tau/alpha)*(2*(n+1)*L-x)).*besseli(0,(1/(2*tau))*sqrt(u.^2-tau/alpha*(2*(n+1)*L-x)^2));
end
I = sqrt(alpha/tau)/k*qt(t-u).*exp(-u/(2*tau)).*I;

end