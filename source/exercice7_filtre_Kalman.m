% Filtre de Kalman

clear; close all; clc;
load('mesurestrajKalm1.mat');

n_of_mesures = length(Tmesu);
n_of_parametres = 6;

var_noise_a = 0.05;
var_noise_distance = sigmesurayon;
var_noise_angle = simesuang;

W = [0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, var_noise_a, 0;
    0, 0, 0, 0, 0, var_noise_a];

H = [1, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0];

X_est = zeros(n_of_parametres, n_of_mesures);
var_P = zeros(1, n_of_mesures);

x_k = zeros(n_of_parametres, 1);
P_k = 1000 * eye(n_of_parametres);

var_P(1, 1) = trace(P_k);

for k = 1:n_of_mesures
    angle = mesures(k, 2);
    D = mesures(k, 1);

    J = [cos(angle), -D * sin(angle); sin(angle), D * cos(angle)];
    V = J * [var_noise_distance, 0; 0, var_noise_angle] * J';

    % Prediction

    T = Tmesu(k);

    F = [1, 0, T, 0, T ^ 2 / 2, 0;
        0, 1, 0, T, 0, T ^ 2 / 2;
        0, 0, 1, 0, T, 0;
        0, 0, 0, 1, 0, T;
        0, 0, 0, 0, 1, 0;
        0, 0, 0, 0, 0, 1];

    x_k = F * x_k;

    P_k = F * P_k * F' + W;

    % Correction

    K_kp = P_k * H' * inv(H * P_k * H' + V);

    pseudo_mesures = [D * cos(angle); D * sin(angle)];
    x_k = x_k + K_kp * (pseudo_mesures - H * x_k);
    P_k = (eye(n_of_parametres) - K_kp * H) * P_k;

    % Sauvegarder le résultat
    X_est(:, k) = x_k;
    var_P(1, k) = trace(P_k);
end

% Vecteur d'erreur
erreur = zeros(1, n_of_mesures);
for k = 1:n_of_mesures
    angle = mesures(k, 2);
    D = mesures(k, 1);
    pseudo_mesures = [D * cos(angle); D * sin(angle)];
    estimated = H * X_est(:,k);
    erreur(k) = sqrt(sum((pseudo_mesures - estimated).^2));
end

tiledlayout(2, 1);
nexttile
plot(1:n_of_mesures, var_P);
legend({'Trace de la matrice de variance'});
nexttile
plot(1:n_of_mesures, erreur);
legend({'Erreur quadratique'});

