clear;
close all;

%% Impact d’un canal de propagation multitrajets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARTIE 2.2.
% Question 2.

%% Données

nbBits = 1000;
Fe = 24000;
Rb = 3000;
Tb = 1 / Rb;
Te = 1 / Fe;

M = 2;

Ts = Tb * log2(M);
Ns = Ts / Te;

%% Information binaire à transmettre, mapping, échantillonage

bitsIn = randi([0, 1], 1, nbBits);
symboles = 2 * bitsIn - 1;
echantillons = kron(symboles, [1 zeros(1, Ns - 1)]);

%% Modulation & démodulation

h = ones(1, Ns);
hr = h;
n0 = length(h);

g = conv(h, hr);
z = filter(g, 1, echantillons);
zEch = z(n0 : Ns : length(symboles) * Ns);

bitsOut = (zEch > 0);
TEB = abs(sum(bitsIn ~= bitsOut)) / nbBits;