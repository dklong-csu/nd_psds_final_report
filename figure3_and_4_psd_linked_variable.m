clc
clear variables
close all

%%  Particle Size Distribution of Linked Variables
%   Theory
%   Rather than the PSD with respect to typical size parameters, sometimes
%   it is more desireable to have a PSD with respect to some combination of
%   the measured sizes. For instance, aspect ratio and specific surface
%   area are meaningful quantities for many particles.
%
%   To compute this particle size distribution, the normal PSD (e.g. q0) is
%   integrated such that the integration bounds remain in the region less
%   than the linked variable for a specified value of the linked variable.
%   This computes the CDF of the linked variable PSD. Taking the derivative
%   the results in the desired PSD.
%
%   This can be done numerically. Here we demonstrate this procedure for
%   the case of q0 being a log-normal distribution, the particles are
%   cylindrical, and we want the PSD with respect to aspect ratio and
%   specific surface area.

%   q0 distribution

avg_d = 20;                                                                 % median diameter in nm
avg_l = 100;                                                                % median length in nm

spread_d = 1.10;                                                            % 68% of prob in [avg_d/spread_d, avg_d*spread_d] and 95% in [avg_d/spread_d^2, avg_d*spread_d^2]
spread_l = 1.10;                                                            % 68% of prob in [avg_l/spread_l, avg_l*spread_l] and 95% in [avg_l/spread_l^2, avg_l*spread_l^2]
corr     = 0.3;                                                             % Correlation between parameters for the 2D distribution

% Log-normal distribution parameters
mu = log([avg_d, avg_l]);
Sigma = [log(spread_d), corr*log(spread_d)*log(spread_l);
    corr*log(spread_d)*log(spread_l), log(spread_l)];

% Function evaluating q0
q0 = @(d,l) fcn_LN2D(d,l,mu,Sigma);

%%  Labels of subfigures

nIDs = 3;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); 


%%  Aspect Ratio -- Computations
%   Aspect ratio is defined as: v = l/d
%   Thus we can relate: l = v*d
%   For a given aspect ratio v, we compute
%   CDF(v) = int_0^Inf int_0^(v*d) q0(l,d)dl dd
%   Then once this is done for all (reasonable) aspect ratios, the
%   derivative can be evaluated numerically

% Which aspect ratios to solve for
v = linspace(1e-6,25,1000);
CDF_v = 0*v;
for ii=1:length(v)
    CDF_v(ii) = integral2(q0,1e-6,Inf,0,@(x) v(ii)*x);
end

PSD_v = diff(CDF_v)./diff(v);


%%  Aspect Ratio -- Figures

text_scale = 1.75;

%   Plot of original distribution
dlims = linspace(1e-6,50);
llims = linspace(1e-6,250);
[D,L] = meshgrid(dlims,llims);
q0_vals = q0(D,L);
asp_ratio = L ./ D;

figure('Position',[105,440,1480,469])
t = tiledlayout("horizontal",...
    "TileSpacing","compact",...
    "Padding","compact");

ax1 = axes(t);
box on
ax2 = copyobj(ax1,t);
hold on
contourf(ax1,D,L,q0_vals, ...
    linspace(1e-4,1e-3),...
    'LineStyle','none')
colormap(ax1,"parula")

[C,h] = contour(ax2,D,L,asp_ratio, ...
    1.5:1.5:6.5,...
    'FaceAlpha',0.2,...
    'ShowText','on',...
    'LabelFormat',"AR = %0.1f",...
    'LineWidth',2,...
    'LabelSpacing',500);
colormap(ax2,"cool")
clabel(C,h,'FontSize',18)

xlabel('Diameter / nm')
ylabel('Length / nm')

xticks(0:10:max(dlims))
yticks(0:50:max(llims))

cb = colorbar(ax1);
xlim([0,max(dlims)])
ylim([0,max(llims)])
text(0.025,0.95,charlbl{1},'Units','normalized','FontSize',12)

ax2.UserData = linkprop([ax1,ax2],...
    {'Position','InnerPosition','DataAspectRatio','xtick','ytick', ...
    'ydir','xdir','xlim','ylim'}); % add more props as needed
ax2.Visible = 'off';
hold off

%   Plot of CDF
nexttile
box on
plot(v,CDF_v,...
    'LineWidth',2)
xlabel('Aspect ratio / -')
ylabel('CDF / -')
text(0.025,0.95,charlbl{2},'Units','normalized','FontSize',12)

%   Plot of PSD
nexttile
box on
plot(v(1:end-1),PSD_v,...
    'LineWidth',2)
xlabel('Aspect ratio / -')
ylabel('PSD / -')
text(0.025,0.95,charlbl{3},'Units','normalized','FontSize',12)
fontsize(gcf,scale=text_scale)

print('figure3','-dsvg')
print(gcf,'-depsc','-image','figure3')

%% Specific Surface Area -- Computations
%   Specific Surface Area is defined as: ssa = SA/V
%   For a cylinder we have
%   SA = pi*d*l + 1/2*pi*d^2
%   V  = pi/4 *d^2*l
%   => ssa = 4/d + 2/l
%   => l = (2d)/(SSA*d-4)

% Which aspect ratios to solve for
% step_size = 0.1;
% ssa = 0:step_size:5;
ssa = linspace(1e-6,1,1000);
CDF_ssa = 0*ssa;
for ii=1:length(ssa)
    CDF_ssa(ii) = integral2(q0, ...
        4/ssa(ii),Inf, ...
        @(x) 2*x./(ssa(ii)*x-4),Inf);
end
CDF_ssa = CDF_ssa;

PSD_ssa = diff(CDF_ssa)./diff(ssa);

%% Specific Surface Area -- Figures

text_scale = 1.75;

%   Plot of original distribution
dlims = linspace(1e-6,50);
llims = linspace(1e-6,250);
[D,L] = meshgrid(dlims,llims);
q0_vals = q0(D,L);
SSA = 4./D + 2./L;

figure('Position',[105,440,1480,469])
t = tiledlayout("horizontal",...
    "TileSpacing","compact",...
    "Padding","compact");

ax1 = axes(t);
box on
ax2 = copyobj(ax1,t);
hold on
contourf(ax1,D,L,q0_vals, ...
    linspace(1e-4,1e-3),...
    'LineStyle','none')
colormap(ax1,"parula")

[C,h] = contour(ax2,D,L,SSA, ...
    0.15:0.1:0.45,...
    'FaceAlpha',0.2,...
    'ShowText','on',...
    'LabelFormat',"SSA = %0.2f",...
    'LineWidth',2,...
    'LabelSpacing',500);
colormap(ax2,"cool")
clabel(C,h,'FontSize',18)

xlabel('Diameter / nm')
ylabel('Length / nm')

xticks(0:10:max(dlims))
yticks(0:50:max(llims))

cb = colorbar(ax1);
xlim([0,max(dlims)])
ylim([0,max(llims)])
text(0.025,0.95,charlbl{1},'Units','normalized','FontSize',12)

ax2.UserData = linkprop([ax1,ax2],...
    {'Position','InnerPosition','DataAspectRatio','xtick','ytick', ...
    'ydir','xdir','xlim','ylim'}); % add more props as needed
ax2.Visible = 'off';
hold off

%   Plot of CDF
nexttile
box on
plot(ssa,CDF_ssa,...
    'LineWidth',2)
xlabel('Specific Surface Area / 1/nm')
ylabel('CDF / -')
text(0.025,0.95,charlbl{2},'Units','normalized','FontSize',12)

%   Plot of PSD
nexttile
box on
plot(ssa(1:end-1),PSD_ssa,...
    'LineWidth',2)
xlabel('Specific Surface Area / 1/nm')
ylabel('PSD / nm')
text(0.025,0.95,charlbl{3},'Units','normalized','FontSize',12)

fontsize(gcf,scale=text_scale)

print('figure4','-dsvg')
print(gcf,'-depsc','-image','figure4')


%%  Helper Functions

%   Probability density function for log-normal
function pdf = fcn_LN2D(D,L,mu,sigma)
    X = log(D);
    Y = log(L);
    x = X(:);
    y = Y(:);
    R = [x,y];
    psd = mvnpdf(R,mu,sigma)./exp(x)./exp(y);
    [r,c] = size(X);
    pdf = reshape(psd,[r,c]);
end
