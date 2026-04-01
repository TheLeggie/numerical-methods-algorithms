function [Q, R] = householder(A)
% Projekt 2, Zadanie 09
% Karol Wójcik, 339107
%
% Funkcja wykonuje rozkład QR macierzy A metodą odbić Householdera
%
% WEJŚCIE:
%   A   - macierz kwadratowa, rzeczywista lub zespolona
%
% WYJŚCIE:
%   Q   - macierz ortogonalna uzyskana z kolejnych macierzy Householdera
%   R   - macierz trójkątna górna


[m, n] = size(A);
R = A;
Q = eye(m);

% Pętla po kolumnach
for k = 1:n-1

    % Tworzenie wektora z k-tej kolumny, od przekątnej w dół
    x = R(k:m, k);
    norm_x = norm(x);

    % Tworzenie wektora normalnego v
    % Aby uniknąć błędów numerycznych, dobieramy odpowiednio znak normy,
    % unikniemy odejmowania bliskich liczb
    x1 = x(1);
    if x1 == 0
        a = - norm_x;
    else
        a = -(x1/abs(x1))* norm_x;
    end

    u  = x;
    u(1) = u(1) - a;
    v = u / norm(u); % Ostateczny wektor v

    % Aktualizacja Q i R
    % Macierz Householdera H = I - 2*v*v' - wstawiamy bezpośrednio do Q i R
    % Aby ograniczyć obliczenia, robimy operacje na podmacierzach, reszta
    % elementów macierzy się nie zmienia
    R(k:m, k:n) = R(k:m, k:n) - 2 * (v * (v' * R(k:m, k:n)));
    Q(:, k:m) = Q(:, k:m) - 2 * (Q(:, k:m) * v) * v';

end % function

