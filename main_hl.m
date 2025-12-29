clc, close all, clear

% Check if CME method is included in directory
flag1 = exist('matlab_ilt.m');
flag2 = exist('iltcme.json');
if flag1 == 0 || flag2 == 0
    error(['Code uses CME-based numerical inverse Laplace transform method. ',...
        'Download ''matlab_ilt.m'' and ''iltcme.json'' from ',...
        'https://github.com/ghorvath78/iltcme ', ...
        'and place in current directory.']);
end

% Store figures in subfolder
mkdir Figures
path_name = './Figures/';

% Test problem
tau = 0.001; problem = '1'; sigma_vals = [0,0.05]; te = 0.4; N = 1001;

L = 0.002; % length of sample
k = 222; % thermal conductivity
rhoc = 2700*896; % volumetric heat capacity 
Qinf = 7000; % total amount of heat absorbed
alpha = k/rhoc; % thermal diffusivity
beta = 0.001; % heat pulse parameter
T0 = 0; % initial uniform sample temperature
Nx = 1001; % number of spatial nodes
x = linspace(0,L,Nx)'; % spatial locations
qt = @(t) (Qinf/beta^2)*(tau + (1-tau/beta)*t).*exp(-t/beta); % equation (10) in paper when q(t) = Qinf*t*exp(-t/beta)/beta^2
IQ = 2*beta; % integral(1-Q(t)/Qinf,0,Inf)
Tinf = T0 + Qinf/(rhoc*L); % steady state temperature (for no heat losses only)
h0 = 1e4; % heat transfer coefficient at x = 0
hL = 1e5; % heat transfer coefficient at x = L
maxFnEvals = 10000; % parameter for numerical Laplace inversion
method = 'cme'; % method for numerical Laplace inversion
qtL = @(s) (Qinf*(tau/(s + 1/beta) + (beta - tau)/(beta*(s + 1/beta)^2)))/beta^2; % Laplace transform of qt

AbsTol = 1e-10; RelTol = 1e-6; % integration tolerances (default values)
font_size = 24;
% plot_legend = true;
plot_legend = false;
rng default; rng(3)

%% Temperature fields
figure(1)
tplot = [0.0001,0.0005,0.0012,0.003,0.006,0.01,0.1];
colors = turbo(9); colors(1,:) = []; colors(4,:) = [];
Ta = zeros(Nx,length(tplot));
for j = 1:length(tplot)
    for i = 1:length(x)
        m = @(s) sqrt(s*(1+tau*s)/alpha);
        h0t = @(s) h0*(1+tau*s);
        hLt = @(s) hL*(1+tau*s);
        TL = @(s) T0/s + qtL(s)*((hLt(s)+k*m(s))*exp(-m(s)*x(i)) - (hLt(s)-k*m(s))*exp(-m(s)*(2*L-x(i))))/...
            ((h0t(s)+k*m(s))*(hLt(s)+k*m(s)) - (h0t(s)-k*m(s))*(hLt(s)-k*m(s))*exp(-2*m(s)*L));
        Ta(i,j) = matlab_ilt(TL, tplot(j), maxFnEvals, method);
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
else
    ylim([0,4])
    set(gca,'FontSize',font_size,'Xtick',[0,L/2,L],'XTickLabel',{'$0$','$L/2$','$L$'},'YTick',0:4,...
        'YTickLabel',{'$0$','$1$','$2$','$3$','$4$'},'TickLabelInterpreter','LaTeX')
end
xlabel('$x$','FontSize',font_size+4,'Interpreter','LaTeX');
ylabel('$\widehat{T}(x,t)$','FontSize',font_size+4,'Interpreter','LaTeX');
text(0.6,0.9,['$s_{p} = ',num2str(sqrt(alpha/tau)),'$'],'FontSize',font_size,'Units','Normalized','Interpreter','LaTeX')
text(0.6,0.8,['$t_{p} = ',num2str(L*sqrt(tau/alpha)),'$'],'FontSize',font_size,'Units','Normalized','Interpreter','LaTeX')
print(gcf,[path_name,'flux_temp_distribution_hl',num2str(problem)],'-depsc2')

%% Legend
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
    print(gcf,[path_name,'temp_field_legend_hl'],'-depsc2')
end

%% Synthetic data (perfect signal)
ts = linspace(0,te,N);
dt = te/(N-1);

tp_indx = find(ts>=L*sqrt(tau/alpha),1,'first');

xa = L;
Tperfect = zeros(1,N);
for j = 1:N
    fprintf('Generating synthetic data: %g%% complete\n',round(j/N*100,1));
    m = @(s) sqrt(s*(1+tau*s)/alpha);
    h0t = @(s) h0*(1+tau*s);
    hLt = @(s) hL*(1+tau*s);
    TL = @(s) T0/s + qtL(s)*((hLt(s)+k*m(s))*exp(-m(s)*x(i)) - (hLt(s)-k*m(s))*exp(-m(s)*(2*L-x(i))))/...
        ((h0t(s)+k*m(s))*(hLt(s)+k*m(s)) - (h0t(s)-k*m(s))*(hLt(s)-k*m(s))*exp(-2*m(s)*L));
    Tperfect(j) = matlab_ilt(TL, ts(j), maxFnEvals, method);
end
Tperfect(1) = T0;

cnt = 2;
for sigma = sigma_vals
    
    cnt = cnt + 1;
    Tdata = Tperfect + sigma*randn(1,length(ts));
    Tdata(1:tp_indx-1) = T0;

    %% Calculate thermal diffusivity and relaxation time (perfect signal)
    alpha_true = alpha;
    alpha_approx = h0*hL*L*trapz(ts,Tdata-T0)/(rhoc*(Qinf-(h0+hL)*trapz(ts,Tdata-T0)));

    tp_true = L*sqrt(tau/alpha);
    indx = find(Tdata~=T0,1,'first');
    tp = (ts(indx-1)+ts(indx))/2;
    tp_approx = tp;
    tau_true = tau;
    tau_approx = alpha_approx*(tp_approx/L)^2;

    %% Solve Cattaneo heat transfer model with computed values of alpha and tau
    Tfit = zeros(1,N);
    for j = 1:N
        fprintf('Generating fitted curve: %g%% complete\n',round(j/N*100,1));
        m = @(s) sqrt(s*(1+tau_approx*s)/alpha_approx);
        h0t = @(s) h0*(1+tau_approx*s);
        hLt = @(s) hL*(1+tau_approx*s);
        TL = @(s) T0/s + qtL(s)*((hLt(s)+k*m(s))*exp(-m(s)*x(i)) - (hLt(s)-k*m(s))*exp(-m(s)*(2*L-x(i))))/...
            ((h0t(s)+k*m(s))*(hLt(s)+k*m(s)) - (h0t(s)-k*m(s))*(hLt(s)-k*m(s))*exp(-2*m(s)*L));        
        Tfit(j) = matlab_ilt(TL, ts(j), maxFnEvals, method);
    end

    %% Back surface temperature history (plots for Figure 6 of paper)
    figure(cnt)
    plot(ts,Tdata,'LineWidth',2,'Color',colors(2,:))
    hold on
    plot(ts,Tfit,'k--','LineWidth',2)
    plot(ts,zeros(size(ts)),'--','LineWidth',2,'Color',0.6*ones(1,3))
    if isequal(problem,'5') || isequal(problem,'6')
        ylim([0,2])
        set(gca,'FontSize',font_size,'Xtick',[0,te],'XTickLabel',{'$0$','$t_{e}$'},'YTick',[0,1,2],...
            'YTickLabel',{'$0$','$1$','$2$'},'TickLabelInterpreter','LaTeX')        
    else
        ylim([-0.1,1.1])
        set(gca,'FontSize',font_size,'Xtick',[0,te],'XTickLabel',{'$0$','$t_{e}$'},'YTick',[0,1],...
            'YTickLabel',{'$0$','$1$'},'TickLabelInterpreter','LaTeX')
    end
    xlim([0,te])
    xl = xlabel('$t$','FontSize',font_size+4,'Interpreter','LaTeX');
    yl = ylabel('$T(L,t)$','FontSize',font_size+4,'Interpreter','LaTeX');
    text(0.5,0.8,['$\widehat{\tau} = ',num2str(tau_approx/10^(floor(log10(abs(tau_approx)))),'%.4f'),'\times 10^{',...
        num2str(floor(log10(abs(tau_approx))),'%i'),'}$'],'FontSize',font_size,'Units','Normalized','Interpreter','LaTeX')
    text(0.5,0.7,['$\widehat{\alpha} = ',num2str(alpha_approx/10^(floor(log10(abs(alpha_approx)))),'%.4f'),'\times 10^{',...
        num2str(floor(log10(abs(alpha_approx))),'%i'),'}$'],'FontSize',font_size,'Units','Normalized','Interpreter','LaTeX')
    if sigma == 0
        text(0.5,0.9,'Noise-free data','FontSize',font_size,'Units','Normalized','Interpreter','LaTeX')
        print(gcf,[path_name,'perfect_rear_temp_hl',num2str(problem)],'-depsc2')
    else
        text(0.5,0.9,'Noisy data','FontSize',font_size,'Units','Normalized','Interpreter','LaTeX')
        print(gcf,[path_name,'noisy_rear_temp_hl',num2str(problem)],'-depsc2')
    end

end

%% Legend (for Figure 6 of paper)
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
    %print(gcf,[path_name,'temp_history_legend'],'-depsc2')
end