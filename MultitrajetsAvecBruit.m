clear;
close all;

%% Impact d’un canal de propagation multitrajets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARTIE 2.2.
% Question 4.

%% Données

M = 2;

nb = 10000; 
nbBits = nb * log2(M);
Fe = 24000;
Rb = 3000;
Tb = 1 / Rb;
Te = 1 / Fe;

Ts = Tb * log2(M);
Ns = Ts / Te;

attenuation = [1, 1 / 2]; retard = [0, Ts];

%% Information binaire à transmettre, mapping, échantillonage

bitsIn = randi([0, 1], 1, nbBits);
symboles = 2 * bitsIn - 1;
echantillons = kron(symboles, [1, zeros(1, Ns - 1)]);

%% Modulation & démodulation

h = ones(1, Ns);
hr = h;
hc = zeros(1, Ns);
hc(floor(retard(1) / Te) + 1) = attenuation(1);
hc(floor(retard(2) / Te) + 1) = attenuation(2);
n0 = length(h);

xe = filter(h, 1, echantillons);

x = filter(hc, 1, xe);

Px = mean(abs(x).^2);
RSBdB = 0 : 10;
RSB = 10.^(RSBdB ./ 10);
TEBavecBruit = zeros(1, length(RSBdB));
zEch = zeros(nbBits / log2(M), length(RSBdB));

for k = 1 : length(RSBdB)
    % Ajout de bruit
    var = Px * Ns / (2 * log2(M) * RSB(k));
    n = sqrt(var) * randn(1, length(x));
    r = x + n;
    
    % Filtre de réception
    hr = h;
    z = filter(hr, 1, r);
    
    % Echantillonnage
    zEch(:, k) = z(n0 : Ns : end);
    
    % Décision et demapping
    bitsOut = (zEch(:, k) > 0).';
    
    % TEB
    TEBavecBruit(k) = sum(bitsIn ~= bitsOut) / nbBits;
end

figure;
TEBtheorique = 0.5 * qfunc(1.5 * sqrt(2 * RSB)) + 0.5 * qfunc(0.5 * sqrt(2 * RSB));
semilogy(RSBdB, TEBavecBruit, "*");
hold on;
semilogy(RSBdB, TEBtheorique);
title("Tracé du TEB pratique (avec filtrage canal, avec bruit)");
xlabel("E_b / N_0");
ylabel("qfunc(E_b / N_0)");
legend("TEB pratique", "TEB théorique (léger décalage)");