clc
clear variables
close all

%%  Definition of qfeed

avg_d = 20;                                                                 % median diameter in nm
avg_l = 100;                                                                % median length in nm

spread_d = 1.1;                                                            % 68% of prob in [avg_d/spread_d, avg_d*spread_d] and 95% in [avg_d/spread_d^2, avg_d*spread_d^2]
spread_l = 1.1;                                                            % 68% of prob in [avg_l/spread_l, avg_l*spread_l] and 95% in [avg_l/spread_l^2, avg_l*spread_l^2]
corr     = 0.3;                                                             % Correlation between parameters for the 2D distribution

% Log-normal distribution parameters
mu = log([avg_d, avg_l]);
Sigma = [log(spread_d), corr*log(spread_d)*log(spread_l);
    corr*log(spread_d)*log(spread_l), log(spread_l)];

% Function for graphing
fcn_Qfeed = @(d,l) fcn_LN2D(d,l,mu,Sigma);

%%  Define the cut

%   Cut through median size
V_mid = fcn_CylinderVolume(exp(mu(1)),exp(mu(2)));
%   Cut */ some threshold
cut_width = 1.25;
V_lo = V_mid / cut_width;
V_hi = V_mid * cut_width;

%   Cut efficiency function
fcn_T = @(d,l) fcn_CutEfficiency(d,l,V_lo,V_hi);

%%  Define qcoarse

%   Normalizing factor
C = fcn_Cintegral(fcn_T, fcn_Qfeed);
%   Coarse PSD
fcn_Qcoarse = @(d,l) fcn_T(d,l) .* fcn_Qfeed(d,l) / C;

%%  Define qfine

%   Normalizing factor
F = fcn_Fintegral(fcn_T, fcn_Qfeed);
%   Fine PSD
fcn_Qfine = @(d,l) (1-fcn_T(d,l)) .* fcn_Qfeed(d,l) / F;

%%  Create Plots

nIDs = 3;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); 

mymap = [1,1,1;parula];

figure('Position',[43,458,1738,420])
t = tiledlayout("horizontal",...
    "TileSpacing","compact",...
    "Padding","loose");

dlims = [5 40];
llims = [5 190];
plotdlims = [0 50];
plotllims = [0 200];

%   Feed PSD
ax1 = axes(t);
box on
ax2 = copyobj(ax1,t);
ax3 = copyobj(ax1,t);
hold on
colormap(ax1,mymap)
fcontour(ax1,fcn_Qfeed,...
    [dlims, llims],...
    'Fill','on',...
    'MeshDensity',200,...
    'LevelList',linspace(1e-4,1e-3))

fcontour(ax2,@fcn_CylinderVolume,...
    [plotdlims, plotllims],...
    'LevelList',[V_lo, V_hi],...
    'LineStyle','-',...
    'LineColor',[0.6350 0.0780 0.1840],...
    'LineWidth',2)
fcontour(ax3,@fcn_CylinderVolume,...
    [plotdlims, plotllims],...
    'LevelList',V_mid,...
    'LineStyle','--',...
    'LineColor','k',...
    'LineWidth',2)

xlim(plotdlims)
ylim(plotllims)
colorbar(ax1)
xlabel("Diameter / nm")
ylabel("Length / nm")
title("Feed")
text(0.025,0.95,charlbl{1},'Units','normalized','FontSize',12)

ax2.UserData = linkprop([ax1,ax2,ax3],...
    {'Position','InnerPosition','DataAspectRatio','xtick','ytick', ...
    'ydir','xdir','xlim','ylim'}); % add more props as needed
ax2.Visible = 'off';
ax3.UserData = linkprop([ax1,ax2,ax3],...
    {'Position','InnerPosition','DataAspectRatio','xtick','ytick', ...
    'ydir','xdir','xlim','ylim'}); % add more props as needed
ax3.Visible = 'off';
hold off

%   Coarse PSD
axC1 = nexttile;
box on
axC2 = copyobj(axC1,t);
axC3 = copyobj(axC1,t);
hold on
colormap(axC1,mymap)
fcontour(axC1,fcn_Qcoarse,...
    [dlims, llims],...
    'Fill','on',...
    'MeshDensity',200,...
    'LevelList',linspace(1e-4,1.5e-3))

fcontour(axC2,@fcn_CylinderVolume,...
    [plotdlims, plotllims],...
    'LevelList',[V_lo, V_hi],...
    'LineStyle','-',...
    'LineColor',[0.6350 0.0780 0.1840],...
    'LineWidth',2)
fcontour(axC3,@fcn_CylinderVolume,...
    [plotdlims, plotllims],...
    'LevelList',V_mid,...
    'LineStyle','--',...
    'LineColor','k',...
    'LineWidth',2)

xlim(plotdlims)
ylim(plotllims)
colorbar(axC1)
xlabel("Diameter / nm")
ylabel("Length / nm")
title("Coarse")
text(0.025,0.95,charlbl{2},'Units','normalized','FontSize',12)

axC2.UserData = linkprop([axC1,axC2,axC3],...
    {'Position','InnerPosition','DataAspectRatio','xtick','ytick', ...
    'ydir','xdir','xlim','ylim'}); % add more props as needed
axC2.Visible = 'off';
axC3.UserData = linkprop([axC1,axC2,axC3],...
    {'Position','InnerPosition','DataAspectRatio','xtick','ytick', ...
    'ydir','xdir','xlim','ylim'}); % add more props as needed
axC3.Visible = 'off';
hold off

%   Fine PSD
axF1 = nexttile;
box on
axF2 = copyobj(axF1,t);
axF3 = copyobj(axF1,t);
hold on
colormap(axF1,mymap)
fcontour(axF1,fcn_Qfine,...
    [dlims, llims],...
    'Fill','on',...
    'MeshDensity',200,...
    'LevelList',linspace(1e-4,2e-3))

fcontour(axF2,@fcn_CylinderVolume,...
    [plotdlims, plotllims],...
    'LevelList',[V_lo, V_hi],...
    'LineStyle','-',...
    'LineColor',[0.6350 0.0780 0.1840],...
    'LineWidth',2)
fcontour(axF3,@fcn_CylinderVolume,...
    [plotdlims, plotllims],...
    'LevelList',V_mid,...
    'LineStyle','--',...
    'LineColor','k',...
    'LineWidth',2)

xlim(plotdlims)
ylim(plotllims)
colorbar(axF1)
xlabel("Diameter / nm")
ylabel("Length / nm")
title("Fine")

text(0.025,0.95,charlbl{3},'Units','normalized','FontSize',12)

axF2.UserData = linkprop([axF1,axF2,axF3],...
    {'Position','InnerPosition','DataAspectRatio','xtick','ytick', ...
    'ydir','xdir','xlim','ylim'}); % add more props as needed
axF2.Visible = 'off';
axF3.UserData = linkprop([axF1,axF2,axF3],...
    {'Position','InnerPosition','DataAspectRatio','xtick','ytick', ...
    'ydir','xdir','xlim','ylim'}); % add more props as needed
axF3.Visible = 'off';
hold off


fontsize(gcf,scale=1.75)

print("figure10","-dsvg")
print(gcf,'-depsc','-image',"figure10");

%%  Helper Functions
function V = fcn_CylinderVolume(d,l)
    V = pi/4 * d.^2 .* l;
end

function T = fcn_CutEfficiency(d,l,v_lower,v_upper)
    V = fcn_CylinderVolume(d,l);
    m = 1/(v_upper-v_lower);
    b = -v_lower/(v_upper-v_lower);
    T = min(1,max(0,m*V+b));
end

function C = fcn_Cintegral(fcn_T,fcn_Qfeed)
    fcn_integrand = @(d,l) fcn_T(d,l).*fcn_Qfeed(d,l);
    C = integral2(fcn_integrand,...
        1e-6,Inf,1e-6,Inf);
end

function F = fcn_Fintegral(fcn_T,fcn_Qfeed)
    fcn_integrand = @(d,l) (1-fcn_T(d,l)).*fcn_Qfeed(d,l);
    F = integral2(fcn_integrand,...
        1e-6,Inf,1e-6,Inf);
end

