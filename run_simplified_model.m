parameters.A = 3.6e5;
parameters.n = 3.6;
parameters.ktr = 0.5/3600;

parameters.e = 0.1;
parameters.u = 1e-3;
parameters.g0 = 0;
parameters.g1 = 0.4*3600;

sol = simplified_model(parameters,[0.1:0.1:10]*3600);

% figure(1000)
% set(1000, 'Name', 'Tryptic')
% plot([0.1:0.1:10]*3600, sol(1,:))

figure(1)
hold on
max( (sol(1,:)./(parameters.g0+parameters.g1*[0.1:0.1:10]*3600)) )
plot([0.1:0.1:10]*3600, ...
    (sol(1,:)./(parameters.g0+parameters.g1*[0.1:0.1:10]*3600)) /  ...
     max( (sol(1,:)./(parameters.g0+parameters.g1*[0.1:0.1:10]*3600)) ), ...
    'r')

figure(2)
hold on
max( (sol(1,:)./(parameters.g0+parameters.g1*[0.1:0.1:10]*3600)) )
plot([0.1:0.1:10]*3600, ...
    (sol(1,:)./(parameters.g0+parameters.g1*[0.1:0.1:10]*3600)) /  ...
     max( (sol(1,:)./(parameters.g0+parameters.g1*[0.1:0.1:10]*3600)) ), ...
    'r --')


figure(3)
hold on
max(sol(2,:))
plot([0.1:0.1:10]*3600, sol(2,:)/max(sol(2,:)), 'r --')
