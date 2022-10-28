clear all;
addpath('./Functions/');

% System: dx/dt = (A+uB)x+bu, y = Cx, dim x = n

n = 2;
A = [0, 1; -1, 0];
B = [0, 1; -1, 0];
b = [1; 0];
C = [0, 1];

% Experiment time
tmax = 50;

% Tolerance of numerical resolution
tol = 10^(-5);
options = odeset('RelTol',tol);

% Initial conditions
x0 = [-10; 0];%[-10; 0];
xhat0 = [-15; 5];%zeros(n, 1);x0;
S = eye(n, n); % dS/dt = -(A+uB)*S - S*(A+uB) - theta(t)*S + C'*C;
Sinv = inv(S);
Res = eye(n, n);
Gram = zeros(n, n); % Gram(t-DeltaT, t)

phase = 0; % 0 = (observation, duration = tobs) / 1 = (stabilization, duration = tstab) / 2 = (stabilization, while Gram(t-DeltaT, t)>gmin)
tobs = 2;
tstab = 3;
DeltaT = 1;
gmin = 0.5*10^(-3);

alpha = 1; % Gain of S during observation
beta = 1; % Gain of S during stabilization
gamma = 10; % Gain of Gram

Switch = [0;1]; % Switching times [time; index]

Z0 = [x0; xhat0; matrix2vect(Sinv); matrix2vect(S); matrix2vect(Res); matrix2vect(Gram); phase];
t = 0;

while t<tmax

    if phase == 0
        if t == 0
            sol = dde23(@(t,z,zpast) observer(t, z, zpast, DeltaT, n, alpha, gamma, A, B, b, C),DeltaT,Z0,[t, min(t+tobs,tmax)],options);
        else
            sol = dde23(@(t,z,zpast) observer(t, z, zpast, DeltaT, n, alpha, gamma, A, B, b, C),DeltaT,sol,[t, min(t+tobs,tmax)],options);
        end
        phase = 1;
        t = sol.x(end);
        it = size(sol.x, 2);
        Switch = [Switch, [t;it]];
        sol.y(end, end) = phase;
    
    else
        if phase == 1
        sol = dde23(@(t,z,zpast) observer(t, z, zpast, DeltaT, n, beta, gamma, A, B, b, C),DeltaT,sol,[t, min(t+tstab,tmax)],options);
        phase = 2;
        t = sol.x(end);
        it = size(sol.x, 2);
        sol.y(end, end) = phase;

        else
        sol = dde23(@(t,z,zpast) observer(t, z, zpast, DeltaT, n, beta, gamma, A, B, b, C),DeltaT,sol,[t, min(t+tstab,tmax)],options);
        t = sol.x(end);
        i = it;
        while i<= size(sol.x,2) && phase ~=0
            Gram = vect2matrix(sol.y((2*n+3*n^2+1):(2*n+4*n^2), i));
            eigmin = min(eig(Gram));
            if eigmin < gmin
                sol.x = sol.x(:, 1:it);
                sol.y = sol.y(:, 1:it);
                sol.yp = sol.yp(:, 1:it);
                sol.discont = [];
                phase = 0;
                sol.y(end, end) = phase;
                Switch = [Switch, [t;it]];
                t = sol.x(end);
            end
            i = i+1;
        end
        end
    end

end


T = (sol.x)'; % Times
Nt = size(T, 1);

Z = (sol.y)';

X = Z(:, 1:n)'; % State
Xhat = Z(:, (n+1):2*n)'; % Observer

Phase = Z(:, end); % Phases

S = zeros(n, n, Nt);
TrS = zeros(1, Nt);
EigS = zeros(n, Nt);
Eigmin = zeros(1, Nt);
Res = zeros(n^2, Nt);
for i = 1:Nt
    S(:, :, i) = vect2matrix(Z(i, (2*n+n^2+1):(2*n+2*n^2)));
    TrS(i) = trace(S(:, :, i));
    EigS(:, i) = eig(S(:, :, i));
    Gram = vect2matrix(Z(i, (2*n+3*n^2+1):(2*n+4*n^2)));
    Eigmin(i) = min(eig(Gram));
    Res(:, i) = Z(i, (2*n+2*n^2+1):(2*n+3*n^2));
end

Eps = zeros(n, 1, Nt); % Error
Eps(:,1,:) = Xhat - X;
SEps = pagemtimes(S, Eps);
EpsSEps = pagemtimes(Eps, 'transpose', SEps, 'none');

Norm = X(1, :).^2 + X(2, :).^2; % |X|^2 + |Xhat|^2
NormEps = zeros(1, Nt);
NormEps(1,:) = (sum(Eps(:,1,:).^2)); % |Eps|^2
NormEpsSEps = zeros(1, Nt);
NormEpsSEps(1,:) = EpsSEps(1,1,:); % Eps'SEps
P1 = polyfit(T(1:floor(end/2)), log10(NormEpsSEps(1:floor(end/2))), 1); % Linear regression

U = cont(Phase, Xhat); % Control


%% Plot
% close all;

figure(1)
clf
plot(X(1, :), X(2, :), 'k-')
hold on
plot(Xhat(1, :), Xhat(2, :), 'k-.')

i = Switch(2, 1);
scatter(X(1, i), X(2, i), 'k', 'o')
i = Switch(2, 2);
scatter(X(1, i), X(2, i), 'k', 'd')
i = 0;
for ts = Switch
    if mod(i, 2) == 0
        if i>0
            scatter(X(1, ts(2)), X(2, ts(2)), 'k', 'o')
        end
        scatter(Xhat(1, ts(2)), Xhat(2, ts(2)), 'k', 'o')
    else
        if i>1
            scatter(X(1, ts(2)), X(2, ts(2)), 'k', 'd')
        end
        scatter(Xhat(1, ts(2)), Xhat(2, ts(2)), 'k', 'd')
    end
    i = i+1;
end

xlabel('$x_1$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
legend({'$(x_1(t), x_2(t))$', '$(\hat x_1(t), \hat x_2(t))$', 'stabilization $\rightarrow$ observation', 'observation $\rightarrow$ stabilization'}, 'Interpreter', 'latex', 'Location', 'southwest')
xlim([-16, 12])
ylim([-15, 12])

figure(2)
clf
plot(T, Eigmin, 'k-')
hold on
plot(T, -1.5*10^(-3) + (Z(:, end)==0)*1*10^(-3), 'k-','LineWidth',2)
ylim([-3*10^(-3), 3*10^(-3)])
xlabel('Time $t$', 'Interpreter', 'latex')
legend('Observability Gramian', 'Mode (observation/stabilization)', 'Interpreter', 'latex', 'Location', 'southeast')


figure(3)
clf
plot(T, Norm, 'k-')
hold on
plot(T, NormEps, 'k-.')
ylim([10^(-15), 10^3])

i = Switch(2, 1);
scatter(T(i), Norm(i), 'k', 'o')
i = Switch(2, 2);
scatter(T(i), Norm(i), 'k', 'd')
i=0;
for ts = Switch
    if mod(i, 2) == 0
        if i>0
            scatter(T(ts(2)), Norm(ts(2)), 'k', 'o')
        end
        scatter(T(ts(2)), NormEps(ts(2)), 'k', 'o')
    else
        if i>1
            scatter(T(ts(2)), Norm(ts(2)), 'k', 'd')
        end
        scatter(T(ts(2)), NormEps(ts(2)), 'k', 'd')
    end
    i = i+1;
end
xlabel('Time $t$', 'Interpreter', 'latex')
legend({'$|x(t)|^2$', '$|\varepsilon(t)|^2$'}, 'Interpreter', 'latex', 'Location', 'southeast')
set(gca,'yscale','log')

% Article format

% define figure properties
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.width      = 8;
opts.height     = 6;
opts.fontType   = 'Times';
opts.fontSize   = 9;

for i = 1:3
fig = figure(i);
% scaling
fig.Units               = 'centimeters';
fig.Position(3)         = opts.width;
fig.Position(4)         = opts.height;
% set text properties
set(fig.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);
% remove unnecessary white space
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
if i==1
saveas(gcf,'Figures/traj.eps')
else
if i==2
saveas(gcf,'Figures/gram.eps')
else
saveas(gcf,'Figures/norm.eps')
end
end
end

%% Other plots

% figure
% clf
% plot(T, U)
% title('u')

% figure
% plot(T, TrS)
% title('Tr S')
% figure
% plot(T, EigS)
% title('Eigs S')

% figure
% clf
% plot(T, X(1, :))
% hold on
% plot(T, Xhat(1, :))
% plot(T, X(2, :))
% plot(T, Xhat(2, :))
% title('x(1), xhat(1)')
% for ts = Switch
%     plot([ts(1), ts(1)], [-10, 10])
% end

% figure
% clf
% plot(T, Z(:, end))
% title('Phase')

% figure
% clf
% plot(T, Res)
% title('Res')

% figure
% clf
% plot(T)
% title('T')

