clc
clear variables
close all

%%  Log scale plot or linear scale plot
use_log = true;

%%

nIDs = 6;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); 

figure('Position',[2,296,1915,521])
t = tiledlayout("horizontal",...
    "TileSpacing","tight",...
    "Padding","loose");

%%  Plot settings
avg_x = 10;                                                                
avg_y = 20;                                                               

spread_x = 1.5;                                                            
spread_y = 1.5; 

plot_range = [1e0 20 1e0 40];
plot_range_log = [1e0 100 1e0 100];
levellist = linspace(1e-3,1e-2,500);
mymap = [1,1,1;parula];

%% uncorrelated
nexttile
                                                           
corr = 0.0;                                                             

% Log-normal distribution parameters
mu = log([avg_x, avg_y]);
Sigma = [log(spread_x), corr*log(spread_x)*log(spread_y);
    corr*log(spread_x)*log(spread_y), log(spread_y)];

% Function for graphing
q0 = @(x,y) fcn_LN2D(x,y,mu,Sigma);

colormap(mymap)
fcn_range = plot_range .* [1.05 0.95 1.05 0.95];
fcontour(q0,...
    fcn_range,...
    'Fill','on',...
    'LevelList',levellist)

colorbar
xlabel('x')
ylabel('y')
title('\sigma_{xy}=0')
text(0.025,0.9,charlbl{1+3*use_log},'Units','normalized','FontSize',12)
if use_log
    xlim(plot_range_log(1:2))
    ylim(plot_range_log(3:4))
    set(gca,'YScale','log')
    set(gca,'XScale','log')
end

%% positive correlation
nexttile                                                          
corr     = 0.9;                                                             

% Log-normal distribution parameters
mu = log([avg_x, avg_y]);
Sigma = [log(spread_x), corr*log(spread_x)*log(spread_y);
    corr*log(spread_x)*log(spread_y), log(spread_y)];

% Function for graphing
q0 = @(x,y) fcn_LN2D(x,y,mu,Sigma);


colormap(mymap)
fcn_range = plot_range .* [1.05 0.95 1.05 0.95];
fcontour(q0,...
    fcn_range,...
    'Fill','on',...
    'LevelList',levellist)
colorbar
xlabel('x')
ylabel('y')
title('\sigma_{xy}>0')
text(0.025,0.9,charlbl{2+3*use_log},'Units','normalized','FontSize',12)
if use_log
    xlim(plot_range_log(1:2))
    ylim(plot_range_log(3:4))
    set(gca,'YScale','log')
    set(gca,'XScale','log')
end

%% negative correlation
nexttile                                                          
corr     = -0.9;                                                             

% Log-normal distribution parameters
mu = log([avg_x, avg_y]);
Sigma = [log(spread_x), corr*log(spread_x)*log(spread_y);
    corr*log(spread_x)*log(spread_y), log(spread_y)];

% Function for graphing
q0 = @(x,y) fcn_LN2D(x,y,mu,Sigma);


colormap(mymap)
fcn_range = plot_range .* [1.05 0.95 1.05 0.95];
fcontour(q0,...
    fcn_range,...
    'Fill','on',...
    'LevelList',levellist)
colorbar
xlabel('x')
ylabel('y')
title('\sigma_{xy}<0')
text(0.025,0.9,charlbl{3+3*use_log},'Units','normalized','FontSize',12)
if use_log
    xlim(plot_range_log(1:2))
    ylim(plot_range_log(3:4))
    set(gca,'YScale','log')
    set(gca,'XScale','log')
end

%%  Save plot

fontsize(gcf,scale=2.25)
if use_log
    figname = "figure1_log";
else
    figname = "figure1_linear";
end
print(gcf,figname,"-dsvg")
print(gcf,'-depsc','-image',figname);