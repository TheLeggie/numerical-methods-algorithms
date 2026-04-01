function q = triangular_quadrature(f, x1, y1, x2, y2, x3, y3, n)
% Projekt 1
% Karol Wójcik
% 
% Funkcja oblicza przybliżoną wartość całki podwójnej z funkcji f(x,y)
% po trójkącie o wierzchołkach (x1,y1), (x2,y2), (x3,y3).
% Zastosowano złożoną kwadraturę rzędu 3-go (reguła punktów środkowych
% boków). Trójkąt wejściowy jest dzielony na n^2 przystających
% podtrójkątów.
%
% WEJŚCIE:
%   f          - uchwyt do funkcji podcałkowej @(x,y).
%   x1, y1     - współrzędne pierwszego wierzchołka trójkąta.
%   x2, y2     - współrzędne drugiego wierzchołka trójkąta.
%   x3, y3     - współrzędne trzeciego wierzchołka trójkąta.
%   n          - liczba podziałów boku (generuje n^2 podtrójkątów).
%
% WYJŚCIE:
%   q            - obliczone przybliżenie całki.


% Ustawiamy współrzędne kolejnych wierzchołków w wektory
A = [x1; y1];
B = [x2; y2];
C = [x3; y3];

% Wektory boków
AB = B - A;
AC = C - A;

% Obliczenie pola całego trójkąta (0.5 * moduł z wyznacznika)
detJ = abs(AB(1)*AC(2) - AB(2)*AC(1));
BigArea = 0.5 * detJ;

% Pole pojedynczego małego trójkąta (skalowanie 1/n^2)
SmallArea = BigArea / (n^2);

% Pomocnicze do końcowego wyniku
coeff = SmallArea / 3.0;
sum = 0;

% Podówjna pętla po kolejnych punktach podziału na dwóch bokach trójkąta
for i = 0:(n-1)
    % Ograniczenie j <= n-1-i zapewnia pozostanie wewnątrz trójkąta
    % (nie wychodzimy na niepotrzebną część równoległoboku)
    for j = 0:(n-1-i)
        
        % Tworzymy wierzchołki trójkąta
        % P00 = punkt(i, j)
        % P10 = punkt(i+1, j)
        % P01 = punkt(i, j+1)
        p00 = A + (i/n)*AB + (j/n)*AC;
        p10 = A + ((i+1)/n)*AB + (j/n)*AC;
        p01 = A + (i/n)*AB + ((j+1)/n)*AC;
        
        % Tworzymy punkty będące środkami boków trójkąta
        % Wierzchołki: P00, P10, P01. Zawsze wewnątrz obszaru.
        m1 = 0.5 * (p00 + p10);
        m2 = 0.5 * (p10 + p01);
        m3 = 0.5 * (p01 + p00);
        
        % Dodanie do sumy
        sum = sum + f(m1(1), m1(2)) + f(m2(1), m2(2)) + f(m3(1), m3(2));
        
        % 
        % Gdy nie jesteśmy na "brzegu" podziału, uwzględniamy też drugi 
        % trójkąt, który jest drugą połową równoległoboku stworzonego 
        % z naszych punktów oraz punktu (i+1,j+1)
        if (i + j) < (n - 1)
            % Potrzebujemy czwartego punktu P11 = punkt(i+1, j+1)
            p11 = A + ((i+1)/n)*AB + ((j+1)/n)*AC;
            
            % Wierzchołki odwróconego: P10, P11, P01
            % Tworzymy punkty będące środkami boków trójkąta
            m4 = 0.5 * (p10 + p11);
            m5 = 0.5 * (p11 + p01);
            m6 = 0.5 * (p01 + p10); 
            
            % Dodanie do sumy
            sum = sum + f(m4(1), m4(2)) + f(m5(1), m5(2)) + ... 
            f(m6(1), m6(2));
        end
    end
end
 
q = sum * coeff;
end % function