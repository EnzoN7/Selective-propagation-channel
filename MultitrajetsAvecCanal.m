clear;
close all;

%% Impact d’un canal de propagation multitrajets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARTIE 2.2.
% Question 3.

%% Données

M = 2;

nb = 5000; 
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
echantillons = kron(symboles, [1 zeros(1, Ns - 1)]);

%% Modulation & démodulation

h = ones(1, Ns);
hr = h;
hc = zeros(1, Ns);
hc(floor(retard(1) / Te) + 1) = attenuation(1);
hc(floor(retard(2) / Te) + 1) = attenuation(2);
n0 = length(h);

xe = filter(h, 1, echantillons);
r = filter(hc, 1, xe);
z = filter(hr, 1, r);
eyediagram(z(4 * Ns + 1 : end), Ns, Ns);

zEch = z(n0 : Ns : length(symboles) * Ns);
scatterplot(zEch);
bitsOut = (zEch > 0);

TEB = abs(sum(bitsIn ~= bitsOut)) / nbBits;