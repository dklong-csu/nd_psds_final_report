clc
clear variables
close all

%%  Plot TEM data versus normal/lognormal fit

TEM = readmatrix("nanorod_tem_data.txt");
TEM = TEM(:,[2 1]);
lb = min(TEM);
ub = max(TEM);

text_scale = 1.75;

figure('Position',[43,458,1738,420])
tiledlayout("horizontal",...
    "TileSpacing","compact",...
    "Padding","compact");

nIDs = 4;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); 

%%  Normal Distribution

%   Normal distribution fit
mu = mean(TEM);
S = cov(TEM);
MVNpsd = @(d,l) fcn_MVN2D(d,l,mu,S);

%   Create Plot
nexttile
hold on
box on
mymap = [1,1,1;colormap("parula")];
colormap(mymap)
fcontour(MVNpsd,...
    [lb(1)*0.9 ub(1)*1.1 lb(2)*0.9 ub(2)*1.1],...
    'Fill','on',...
    'MeshDensity',200)
scatter(TEM(:,1),TEM(:,2),...
    75,'.',...
    'MarkerFaceAlpha',0.3,...
    'MarkerEdgeAlpha',0.5,...
    'MarkerEdgeColor',[0.6350 0.0780 0.1840],...
    'MarkerFaceColor',[0.6350 0.0780 0.1840])
xlim([lb(1)*0.8 ub(1)*1.15])
ylim([lb(2)*0.8 ub(2)*1.15])
colorbar
xlabel('Diameter / nm')
ylabel('Length / nm')
title("Normal Distribution Fit")
text(0.025,0.95,charlbl{1},'Units','normalized','FontSize',12)
hold off

%   Q-Q Plot
%   Empirical quantiles
dist = zeros(size(TEM,1),1);
N = length(dist);
for ii=1:N
    diff = TEM(ii,:)-mu;
    dist(ii) = diff * (S\diff');
end
qtile = quantile(dist,N);
theory = chi2inv((1:N)/(N+1),2);

nexttile
box on
hold on
plot(theory,theory,'-k','LineWidth',2)
scatter(qtile,theory,36,[0 0.4470 0.7410],"filled")
xlabel("Squared Mahalanobis Distance")
ylabel("\chi^2 Quantile")
title("Normal Distribution Q-Q")
text(0.025,0.95,charlbl{2},'Units','normalized','FontSize',12)
hold off

%%  Lognormal Distribution

%   Lognormal distribution fit
mu = mean(log(TEM));
S = cov(log(TEM));
MVLNpsd = @(d,l) fcn_LN2D(d,l,mu,S);

%   Create Plot
nexttile
hold on
box on
mymap = [1,1,1;colormap("parula")];
colormap(mymap)
% colormap("summer")
% colormap(flipud(colormap("copper")))
fcontour(MVLNpsd,...
    [lb(1)*0.9 ub(1)*1.1 lb(2)*0.9 ub(2)*1.1],...
    'Fill','on',...
    'MeshDensity',200)
scatter(TEM(:,1),TEM(:,2),...
    75,'.',...
    'MarkerFaceAlpha',0.3,...
    'MarkerEdgeAlpha',0.5,...
    'MarkerEdgeColor',[0.6350 0.0780 0.1840],...
    'MarkerFaceColor',[0.6350 0.0780 0.1840])
xlim([lb(1)*0.8 ub(1)*1.15])
ylim([lb(2)*0.8 ub(2)*1.15])
colorbar
xlabel('Diameter / nm')
ylabel('Length / nm')
title("Lognormal Distribution Fit")
text(0.025,0.95,charlbl{3},'Units','normalized','FontSize',12)
hold off

%   Q-Q Plot
%   Empirical quantiles
dist = zeros(size(TEM,1),1);
N = length(dist);
for ii=1:N
    diff = log(TEM(ii,:))-mu;
    dist(ii) = diff * (S\diff');
end
qtile = quantile(dist,N);
theory = chi2inv((1:N)/(N+1),2);
nexttile
box on
hold on
plot(theory,theory,'-k','LineWidth',2)
scatter(qtile,theory,36,[0 0.4470 0.7410],"filled")
xlabel("Squared Mahalanobis Distance")
ylabel("\chi^2 Quantile")
title("Lognormal Distribution Q-Q")
text(0.025,0.95,charlbl{4},'Units','normalized','FontSize',12)
hold off

fontsize(gcf,scale=text_scale)
%%  Print the generated figure
print(gcf,'-dsvg',"figure9")
print(gcf,'-depsc','-image',"figure9");