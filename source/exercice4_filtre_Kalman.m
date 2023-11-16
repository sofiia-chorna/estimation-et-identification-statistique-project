% Filtre de Kalman étendu

clear; close all; clc;
load('mesurestrajKalm2.mat');

n_of_mesures = length(Tmesu);
n_of_parametres = 6 + 2;

var_noise_a = 100;
var_noise_distance = sigmesurayon;
var_noise_angle = simesuang;

W = [0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, var_noise_a, 0, 0, 0;
    0, 0, 0, 0, 0, var_noise_a, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0 ];

X_est = zeros(n_of_parametres, n_of_mesures);
var_P = zeros(1,n_of_mesures);

% Initialization

x_k = zeros(n_of_parametres, 1);
angle = mesures(1, 1);
D = mesures(1, 2);
x_k(1,1) = D * cos(angle);
x_k(2,1) = D * sin(angle);

P_k = 10000 * eye(n_of_parametres);

var_P(1,1) = trace(P_k);

V = [var_noise_distance, 0; 0, var_noise_angle];

for k = 2:n_of_mesures
    % Prediction

    T = Tmesu(k);

    F = [1, 0, T, 0, T ^ 2 / 2, 0, 0, 0;
        0, 1, 0, T, 0, T ^ 2 / 2, 0, 0;
        0, 0, 1, 0, T, 0, 0, 0;
        0, 0, 0, 1, 0, T, 0, 0;
        0, 0, 0, 0, 1, 0, 0, 0;
        0, 0, 0, 0, 0, 1, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0];

    x_coor = x_k(1,1);
    y_coor = x_k(2,1);

    H = [
        x_coor / sqrt(x_coor^2 + y_coor^2), y_coor / sqrt(x_coor^2 + y_coor^2),0,0,0,0,0,0;
        -y_coor / x_coor^2 + y_coor^2, x_coor / x_coor^2 + y_coor^2,0,0,0,0,0,0
    ];

    B_k = zeros(n_of_parametres, 1);
    M_k = -H * x_k + [
        sqrt(x_coor^2 + y_coor^2);
        atan(y_coor/x_coor)
    ];
    if (x_k(1,1) == 0)
        M_k(2,:) = 0;
    end
    B_k(7,:) = M_k(1,:);
    B_k(8,:) = M_k(2,:);
    
    x_k = F * x_k + B_k;
    P_k = F * P_k * F' + W;

    % Correction

    K_kp = P_k * H' * inv(H * P_k * H' + V);
    
    x_k = x_k + K_kp * (mesures(k, :)' - H * x_k);
    P_k = (eye(n_of_parametres) - K_kp * H) * P_k;

    % Sauvegarder le resultat
    X_est(:,k) = x_k;
    var_P(1,k) = trace(P_k);
end

% Vecteur d'erreur
erreur = zeros(1,n_of_mesures);
for k = 1:n_of_mesures
    angle = mesures(k, 1);
    D = mesures(k, 2);
    pseudo_mesures = [D * cos(angle); D * sin(angle)];
    estimated = H * X_est(:,k);
    erreur(k) = sqrt(sum((pseudo_mesures - estimated).^2));
end

tiledlayout(2,1);
nexttile
plot(1:n_of_mesures,var_P);
legend({'Trace de la matrice de variance'});
nexttile
plot(1:n_of_mesures,erreur);
legend({'Erreur quadratique'});