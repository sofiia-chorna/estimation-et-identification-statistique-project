clear; close all; clc;
load('xx2.mat');

Tmesu = Tmes;
mesures = Mesures;

mesures

var_noise_x = 2;
var_noise_y = 2;

n_of_mesures = length(Tmesu);
n_of_parametres = 6;

H = zeros(n_of_mesures * 2, n_of_parametres); % on multiplie par 2, car le vecteur de mesures contient 2 variables
R = zeros(n_of_mesures * 2, n_of_mesures * 2);

for i = 1:n_of_mesures
    H(2*i-1, :) = [1, Tmesu(i), Tmesu(i) ^ 2 / 2, 0, 0, 0];
    H(2*i, :) = [ 0, 0, 0, 1, Tmesu(i), Tmesu(i) ^ 2 / 2];

    R(2*i-1, 2*i-1) = var_noise_x;
    R(2*i, 2*i) = var_noise_y;
end

disp(H);
disp(R);

normal_matrix = H' / R * H;
disp(normal_matrix);

mesures_vector = mesures';
mesures_vector = mesures_vector(:);

theta = normal_matrix \ H' / R * mesures_vector;
disp(theta);

% Plot result

predictions = zeros(n_of_mesures, 2);
for i = 1:n_of_mesures
    H = [1, Tmesu(i), Tmesu(i) ^ 2 / 2, 0, 0, 0;
         0, 0, 0, 1, Tmesu(i), Tmesu(i) ^ 2 / 2];
    predictions(i,:) = H * theta;
end

plot(mesures(:,1),mesures(:,2),'r+');
hold on
plot(predictions(:,1),predictions(:,2),'bx');
legend({'Donné réel','Donné estimé'});
