function [Xp, D, W, Wp, Phi, R] = PCA_custom(X)

[~, cols] = size(X);

% Z scores of data
X = X - mean(X);                                                            
X = X ./ std(X);

% Correlation matrix
W = X' * X;                                                                 

% Right eigenvectors and diagonal matrix of eigenvalues
[V, D] = eig(W);          

% Eigenvalues in decreasing order. Inertie selon les axes.
D = sum(fliplr(D));                                                         

% Eigenvectors in the same order. Axes principaux.
V = fliplr(V);                                                              

% Calcul des C.P.  pour représentation dans la nouvelle base.
Xp = (inv(V) * X')';       

% Normalisation des C.P. pour projection des variables
Xpn = Xp ./ std(Xp);       

% Matrice de correlation Variables / C.P.
W = (X' * Xpn) ./ size(X, 1); 

% Projection sur les 2 premières composantes.
Wp = (X' * Xp(:, 1:2));   
Wp = [(1:cols)' Wp];

% Rayon de correlation.
R = sqrt(W(:, 1) .^ 2 + W(:, 2) .^ 2); 

% Calcul de la phase de chaque variable
Phi = atan(W(:, 2) ./ W(:, 1));    
