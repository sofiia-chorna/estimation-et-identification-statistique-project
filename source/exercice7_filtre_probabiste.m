% Filtre probabiste

clear; close all; clc;
load('mesurestrajKalm3.mat');

n_of_mesures = size(mesures, 1);
T = 0.1;
g = 9.81;
sigma_a = 0.05;
mu = 0.01;
beta = exp(-mu);
sigma_w = sigma_a * sqrt(1 - exp(-2 * mu));

% initialisation
x0_initial = mesures(1, 1) * cos(mesures(1, 6));
y0_initial = mesures(1, 1) * sin(mesures(1, 6));
X_initial = [x0_initial/2; y0_initial/2; 0; 0; 0; -g; 0; 0];

Pi = zeros(8,8,n_of_mesures);
P0 = 10^6 * eye(8);

F = [1, 0, T, 0, T^2/2, 0, 0, 0;
     0, 1, 0, T, 0, T^2/2, 0, 0;
     0, 0, 1, 0, T, 0, 0, 0;
     0, 0, 0, 1, 0, T, 0, 0;
     0, 0, 0, 0, beta, 0, 0, 0;
     0, 0, 0, 0, 0, beta, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0];

H = [0, 0, 0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 0, 0, 1];

W = T * [0, 0, 0, 0, 0, 0, 0, 0;
         0, 0, 0, 0, 0, 0, 0, 0;
         0, 0, 0, 0, 0, 0, 0, 0;
         0, 0, 0, 0, 0, 0, 0, 0;
         0, 0, 0, 0, sigma_w^2, 0, 0, 0;
         0, 0, 0, 0, 0, sigma_w^2, 0, 0;
         0, 0, 0, 0, 0, 0, 0, 0;
         0, 0, 0, 0, 0, 0, 0, 0];

V = [sigmesurayon^2, 0; 0, simesuang^2];

X_est = zeros(8,n_of_mesures);
variance = zeros(1,n_of_mesures);
erreur = zeros(1,n_of_mesures);

X_est_polar = zeros(200,2);
eqm = zeros(5,1);
best_measures = zeros(200,2);

for i = 1:n_of_mesures
    % Prediction
    if i > 1
        X_est(:,i) = F*X_est(:,i-1);
        Pi(:,:,i) = F*Pi(:,:,i-1)*F' + W;
    else
        X_est(:,1) = F*X_initial;
        Pi(:,:,1) = F*P0*F' + W;
    end

    x = X_est(1,i);
    y = X_est(2,i);

    % Linearization matrice H
    Hk = [x/sqrt(x^2+y^2) y/sqrt(x^2+y^2)
        -y/(x^2+y^2) x/(x^2+y^2)];
    H(1:2,1:2) = Hk;

    X_est_polar(i,:) = [sqrt(x^2+y^2); atan(y/x)];
    
    % Mise à jour mk_x et mk_y
    X_est(7:8,i) = -H*X_est(:,i) + [sqrt(x^2+y^2); atan(y/x)];
    
    % Sélection des mesures les plus probables
    xk = X_est_polar(i,1);
    yk = X_est_polar(i,2);
    
    for j = 1:5
    eqm(j) = sum(sqrt(([xk;yk]-[mesures(i,j);mesures(i,j+5)]).^2));
    end
    [erreur(1,i), j_best] = min(eqm);
    best_measures(i,:) = [mesures(i,j_best) mesures(i,j_best+5)];
    
    % Correction
    Ki = Pi(:,:,i)*H'*inv(H*Pi(:,:,i)*H'+V);
    
    % Mise à jour
    X_est(:,i) = X_est(:,i) + Ki*(best_measures(i,:)'-H*X_est(:,i));
    Pi(:,:,i) = (eye(8)-Ki*H)*Pi(:,:,i);
    variance(1,i) = trace(Pi(:,:,i));
end

tiledlayout(2, 1);
nexttile
plot(1:n_of_mesures, variance);
legend({'Trace de la matrice de variance'});

nexttile
plot(1:n_of_mesures, erreur);
legend({'Erreur quadratique'});

