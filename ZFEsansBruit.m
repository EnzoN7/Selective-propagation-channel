clear;
close all;

%% Egalisation ZFE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARTIE 3.2.
% Sans bruit

%% Données

M = 2;

nb = 1000; 
nbBits = nb * log2(M);
Fe = 24000;
Rb = 3000;
Tb = 1 / Rb;
Te = 1 / Fe;

Ts = Tb * log2(M);
Ns = Ts / Te;

attenuation = [1, 1 / 2]; retard = [0, Ts];

%% Information binaire à transmettre, mapping, échantillonage

bitsIn = [1 zeros(1, nbBits - 1)];
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
r = filter(hc, 1, xe);
z = filter(hr, 1, r);
zEch = z(n0 : Ns : end);

% Egalisation
Y0 = [1 zeros(1, nbBits - 1)].';
Z = toeplitz(zEch);
C = Z \ Y0;
y = filter(C, 1, zEch);

bitsOut = (y > 0);

TEB = abs(sum(bitsIn ~= bitsOut)) / nbBits;