data = load("../src/TimeEvolution.txt");

%% psi
real = data(:, 1:2:end-1);
imag = data(:, 2:2:end);

psi = real + 1i.*imag;


%% Visualization
k = 560;
nx = size(psi, 2);

i = k*nx + 1;
j = (k+1) * nx;
figure(1)
imagesc(abs(psi((i:j), :)))
colorbar