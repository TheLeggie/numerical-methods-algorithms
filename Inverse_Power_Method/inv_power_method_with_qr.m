function [lambda, v, errEst] = inv_power_method_with_qr(A, tol, maxIter)
% Projekt 2
% Karol Wójcik
%
% Funkcja znajduje najmniejszą co do wartości bezwzględnej wartość własną
% macierzy kwadratowej A przy użyciu rozkładu QR otrzymanego metodą odbić
% Householdera.
%
% WEJŚCIE:
%   A       - macierz kwadratowa, której najmniejszą wartość własną szukamy
%   tol     - tolerancja dla oszacowania błędu względnego odwróconej
%             metody potęgowej
%   maxIter - maksymalna liczba iteracji w odwróconej metodzie potęgowej
%
% WYJŚCIE:
%   lambda  - przybliżenie najmniejszej wartości własnej
%   v       - wektor własny odpowiadający powyższej wartości własnej
%   errEst  - oszacowanie błędu względnego obliczonego przybliżenia
%             wartości własnej


n = size(A, 1);

% Rozkład QR metodą odbić Householdera
[Q, R] = householder(A);

% Inicjalizacja wektora początkowego - jest losowy
x = rand(n, 1) + 1i * rand(n, 1);

%lambda_prev = 0;
errEst = Inf;

% Pętla odwróconej metody potęgowej
for k =1:maxIter

    x_norm = x / norm(x);

    % Rozwiązanie układu QR * x_new = x_norm
    y = Q' * x_norm;
    x_new = R \ y;

    % Obliczenie przybliżenia wartości własnej
    lambda_inv_curr = x_norm' * x_new;

    % Odwracanie przybliżenia wartości własnej, przybliżenie dla A to
    % odwrotność przybliżenia dla A^-1
    lambda_curr = 1 / lambda_inv_curr;

    v_temp = x_new / norm(x_new);
    % Obliczanie błędu względnego i sprawdzanie czy osągnięto tolerancję
    errEst = norm(A * v_temp - lambda_curr * v_temp ) ...
        / abs(lambda_curr);

    if errEst < tol
        lambda = lambda_curr;
        v = x_new / norm(x_new);
        return;
    end

    x = v_temp;
end

% Ostateczna aktualizacja zwracanych wartości po pętli
lambda = lambda_curr;
v = x / norm(x);

end % function







