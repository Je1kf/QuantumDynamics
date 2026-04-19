%% Time evolution of stationary Gaussian Packet

%%
data = load("../src/stationaryGaussian.txt");

%% psi
real = data(:, 1:2:end-1);
imag = data(:, 2:2:end);

psi = real + 1i.*imag;

x = linspace(-5, 5, 160);
nx = 160;
nt = 1000;

%% Fitting function
gauss = @(x, a, c) a*exp(-(x).^2 / (2*(c^2)));
gaussEqn = 'a*exp(-(x)^2 / (2*(c^2)))';

sigmas = zeros(nt, 1);
for k = 1:999
    k
    i = k*nx + 1;
    j = (k+1) * nx;
    f1 = fit(x', abs(psi(i:j, 80)).^2, gaussEqn);
    % plot(x, abs(psi(i:j, nx/2)).^2, ...
    %      x, gauss(x, f1.a, f1.c), '.')
    sigmas(k) = f1.c * sqrt(2);
end

%% Analytical

dt = 0.00125;

sigma = @(T) sqrt(0.2^4 + T.^2)/0.2;

t = dt:dt:1.25;

figure(1)
clf
hold on
scatter(t, sigmas, 20, 'k')
plot(t, sigma(t), "LineWidth", 2)
hold off
grid on
ylim([0, 20])
xlabel("Time a.u")
legend("Numerical", "Analytical", Location="best")
ylabel("$\sigma$", "Interpreter", "latex")
set(gca, "FontName", "Serif", "FontSize", 14)

% Zoom
figure(2)
clf
hold on
scatter(t, sigmas, 40, 'k')
plot(t, sigma(t), "LineWidth", 2)
hold off
grid on
ylim([0, 1])
xlim([0, 0.2])
xlabel("Time a.u")
legend("Numerical", "Analytical", Location="best")
ylabel("$\sigma$", "Interpreter", "latex")
set(gca, "FontName", "Serif", "FontSize", 14)