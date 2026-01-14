
clear; close all; clc

%% Parameters
D     = 1.5;
Hloss = 1.0;
tMax  = 60;
dt    = 0.01;
sMin  = 0.25;

theta0 = 0.10;
s0     = 1.20;
x0     = 0;
y0     = 0;

thetaStar = -atan(D);

hMax  = 1.0;
kList = [0 0.5 1 2 4 8];

rangeAtH = zeros(size(kList));

% Nominal runs for plots/checks
kMain  = 2;
solBase = simulateGlider(theta0,s0,x0,y0,D,0,hMax,thetaStar,Hloss,dt,tMax,sMin);

%% Check (i): k=0 gives exactly baseline (h(t)=0)
kTest = 0;
h_from_base = arrayfun(@(th) sat(-kTest*(th - thetaStar), hMax), solBase.theta);
fprintf("Check (i): max |h(t)| when k=0 = %.3e\n", max(abs(h_from_base)));

solCtrl = simulateGlider(theta0,s0,x0,y0,D,kMain,hMax,thetaStar,Hloss,dt,tMax,sMin);

%% Check (iii): saturation respected, i.e. |h(t)| <= hMax
hCtrl = arrayfun(@(th) sat(-kMain*(th - thetaStar), hMax), solCtrl.theta);
maxAbsH = max(abs(hCtrl));
fprintf("Check (iii): max |h(t)| for k=%g is %.6f (hMax=%.6f)\n", kMain, maxAbsH, hMax);

nSat = sum(abs(hCtrl) >= hMax - 1e-12);
fprintf("Check (iii): saturation hits = %d out of %d time steps\n", nSat, numel(hCtrl));

assert(maxAbsH <= hMax + 1e-12, "Saturation violated: |h(t)| exceeded hMax");

%% Sweep k 
for i = 1:numel(kList)
    sol = simulateGlider(theta0,s0,x0,y0,D,kList(i),hMax,thetaStar,Hloss,dt,tMax,sMin);
    rangeAtH(i) = sol.x(end);
end

%% Robustness check (vary theta0)
theta0_list = [0.10 0.60 1.00];

kPos = kList(kList > 0);
Apos = rangeAtH(kList > 0);
[~, idxBest] = max(Apos);
k_best = kPos(idxBest);

fprintf("\nRobustness check (vary theta0, keep s0=%.2f): baseline k=0 vs controlled k=k_best=%.3g\n", s0, k_best);
fprintf("theta0     A_baseline(k=0)     A_control(k=k_best)     stop_reason_baseline\n");

for th0 = theta0_list
    solB = simulateGlider(th0,s0,x0,y0,D,0,hMax,thetaStar,Hloss,dt,tMax,sMin);
    solC = simulateGlider(th0,s0,x0,y0,D,k_best,hMax,thetaStar,Hloss,dt,tMax,sMin);

    A_B = solB.x(end);
    A_C = solC.x(end);

    if solB.y(end) <= -Hloss
        reason = "hit y=-H";
    elseif solB.s(end) <= sMin
        reason = "stalled (s<=sMin)";
    else
        reason = "other stop";
    end

    fprintf("%5.2f         %10.4f           %10.4f        %s\n", th0, A_B, A_C, reason);
end
fprintf("\n");

%% Check (ii): halve timestep dt -> dt/2 
dtFine = dt/2;

rangeAtH_dt       = rangeAtH;
rangeAtH_dt_over2 = zeros(size(kList));

for i = 1:numel(kList)
    sol2 = simulateGlider(theta0,s0,x0,y0,D,kList(i),hMax,thetaStar,Hloss,dtFine,tMax,sMin);
    rangeAtH_dt_over2(i) = sol2.x(end);
end

absDiff = abs(rangeAtH_dt_over2 - rangeAtH_dt);
relDiff = absDiff ./ max(1e-12, abs(rangeAtH_dt));

fprintf("\nCheck (ii): timestep refinement (dt vs dt/2) for A(H)=x(t_H)\n");

kcol = kList(:);
A_dt = rangeAtH_dt(:);
A_dt2 = rangeAtH_dt_over2(:);
dAbs = absDiff(:);
dRel = relDiff(:);

assert(numel(kcol)==numel(A_dt) && numel(A_dt)==numel(A_dt2) && numel(A_dt2)==numel(dAbs) && numel(dAbs)==numel(dRel), ...
    "Check (ii) failed: vector lengths don't match.");

T = table(kcol, A_dt, A_dt2, dAbs, dRel);
T.Properties.VariableNames = {'k','A_dt','A_dt_over2','AbsDiff','RelDiff'};
disp(T);

fprintf("Max AbsDiff = %.6g\n", max(absDiff));
fprintf("Max RelDiff = %.6g\n\n", max(relDiff));

%% Save figures
outDir = "images_individual";
if ~exist(outDir, "dir")
    mkdir(outDir);
end

% Figure 1: flight path (x vs y)
figure(1)
plot(solBase.x, solBase.y, 'LineWidth', 2); hold on
plot(solCtrl.x, solCtrl.y, 'LineWidth', 2);
yline(-Hloss, '--');
xlabel('\bf x'); ylabel('\bf y');
legend("Baseline (k=0)", "Autopilot (k=" + kMain + ")", "y=-H", "Location", "best");
grid on
set(gca,'FontSize',14);
set(findall(gcf,'Type','text'),'FontSize',16);
exportgraphics(gcf, fullfile(outDir,"fig_paths.png"), "Resolution", 300);

% Figure 2: phase plane (theta vs s)
figure(2)
plot(solBase.theta, solBase.s, 'LineWidth', 2); hold on
plot(solCtrl.theta, solCtrl.s, 'LineWidth', 2);
xlabel('\bf \theta'); ylabel('\bf s');
legend("Baseline (k=0)", "Autopilot (k=" + kMain + ")", "Location", "best");
grid on
set(gca,'FontSize',14);
set(findall(gcf,'Type','text'),'FontSize',16);
exportgraphics(gcf, fullfile(outDir,"fig_phaseplane.png"), "Resolution", 300);

% Figure 3: range vs k
figure(3)
plot(kList, rangeAtH, '-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('\bf gain k'); ylabel('\bf range at fixed height loss');
grid on
set(gca,'FontSize',14);
set(findall(gcf,'Type','text'),'FontSize',16);
exportgraphics(gcf, fullfile(outDir,"fig_range_vs_k.png"), "Resolution", 300);

disp("Done. Figures saved to: " + outDir)

%% functions

function sol = simulateGlider(theta0,s0,x0,y0,D,k,hMax,thetaStar,Hloss,dt,tMax,sMin)
    nMax = ceil(tMax/dt) + 1;

    t = zeros(nMax,1);
    theta = zeros(nMax,1);
    s = zeros(nMax,1);
    x = zeros(nMax,1);
    y = zeros(nMax,1);

    theta(1) = theta0;  s(1) = s0;  x(1) = x0;  y(1) = y0;

    for n = 1:nMax-1
        z = [theta(n); s(n); x(n); y(n)];

        k1 = rhs(z, D, k, hMax, thetaStar);
        k2 = rhs(z + 0.5*dt*k1, D, k, hMax, thetaStar);
        k3 = rhs(z + 0.5*dt*k2, D, k, hMax, thetaStar);
        k4 = rhs(z + dt*k3,     D, k, hMax, thetaStar);

        zNext = z + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

        theta(n+1) = zNext(1);
        s(n+1)     = zNext(2);
        x(n+1)     = zNext(3);
        y(n+1)     = zNext(4);
        t(n+1)     = t(n) + dt;

        if y(n+1) <= -Hloss, break; end
        if s(n+1) <= sMin,   break; end
        if ~isfinite(s(n+1)) || ~isfinite(theta(n+1)), break; end
    end

    last = find(t > 0, 1, 'last');
    if isempty(last), last = 1; end

    sol.t = t(1:last);
    sol.theta = theta(1:last);
    sol.s = s(1:last);
    sol.x = x(1:last);
    sol.y = y(1:last);
end

function dz = rhs(z, D, k, hMax, thetaStar)
    theta = z(1);
    s     = z(2);

    thetaDot0 = (s^2 - cos(theta))/s;
    h         = sat(-k*(theta - thetaStar), hMax);
    thetaDot  = thetaDot0 + h;

    sDot = -sin(theta) - D*s^2;

    xDot = s*cos(theta);
    yDot = s*sin(theta);

    dz = [thetaDot; sDot; xDot; yDot];
end

function y = sat(u, hMax)
    y = min(max(u, -hMax), hMax);
end