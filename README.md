syms x
format long 
disp('Metodo de Interpolación de Lagrange');
X = input('Ingrese los valores de x: ');
Y = input('Ingrese los valores de Y: ');

aprox = input('Ingrese el valor a aproximar: ');

n=length(X);
for i=1:n
    numerador=1;
    denominador=1;
    for j=1:n
        if j~=i
            numerador=numerador*(x-X(j));
            denominador=denominador*(X(i)- X(j));
        end
    end
    L(i) = numerador/denominador;
    fprintf('L%d(x)= \n', i-1);
    pretty(L(i))
end
%Constuyendo el polinomio
pol = 0;
for i=1:n
    pol=pol + double(Y(i))*L(i);
end
fprintf('Polinomio de lagrange resultante: \n ');
pretty(vpa(pol,9))
valAprox = subs(pol, aprox);
fprintf('El valor aproximado de la función es: %.15f \n ', double(valAprox));


----------------------------------------------------------------------------

syms x
format long
disp('Metodo de Neville')
X = input('Ingrese los valores de X: ');
Y = input('Ingrese los valores de Y: ');
aprox = input('Ingrese el valor a aproximar: ');

n=length(X);
Q=zeros(n);

Q(:,1)= Y';

for j=2:n
    for i=j:n
        Q(i,j)= ((aprox-X(i-j+1))*(Q(i,j-1))-(aprox - X(i))*(Q(i-1,j-1)))/(X(i)-X(i-j+1));
    end
end

fprintf('La matriz de valores Q: \n')
display(Q)
fprintf('El valor aproximado es: %.15f', double(Q(n,n)))


----------------------------------------------------------------------------

syms x
disp('Metodo de diferencias divididas de Newton');

X = input('Ingrese los valores de x como un vector: ');
Y = input('Ingrese los valores de y como un vector: ');
aprox = input('Ingrese el valor a aproximar: ');

n = length(X);

F = zeros(n); 
F(:, 1) = Y';  

for j = 2:n
    for i = j:n
        F(i, j) = (F(i, j-1) - F(i-1, j-1)) / (X(i) - X(i-j+1));
    end
end

fprintf('Tabla de diferencias divididas:\n');
disp(F);

syms x; 
pol = F(1, 1);
for i = 2:n
    factor = 1; 
    for j = 1:i-1
        factor = factor * (x - X(j)); 
    end
    pol = pol + factor*F(i,i); 
end


% Mostrar el polinomio de Newton
fprintf('El polinomio de Newton es:\n');
pretty(vpa(pol, 9)); 

% Evaluar el polinomio en el punto a aproximar
valorAprox = double(subs(pol, x, aprox)); % Evaluar P en x = aprox
fprintf('El valor aproximado en x = %.4f es: %.9f\n', aprox, valorAprox);


----------------------------------------------------------------------------
syms x
format long

disp('Metodo de Hermite')
X = input('Introduce el arreglo X (valores de x): ');
Y = input('Introduce el arreglo Y (valores de f(x)): ');
DY = input('Introduce el arreglo DY (valores de f''(x)): ');
x_interp = input('Introduce el valor a interpolar: ');

n = length(X);
N = 2 * n;

Q = zeros(N, N);

Z = zeros(1, N);

for i = 1:n
    Z(2*i-1) = X(i); 
    Z(2*i) = X(i);

    Q(2*i-1, 1) = Y(i);  
    Q(2*i, 1) = Y(i);
end

 
for i = 1:N

    if mod(i, 2) == 0
        Q(i, 2) = DY(i/2);
    else      
        if i > 1
            Q(i, 2) = (Q(i, 1) - Q(i-1, 1)) / (Z(i) - Z(i-1));
        end
    end

end

for j = 3:N
    for i = j:N
        Q(i, j) = (Q(i, j-1) - Q(i-1, j-1)) / (Z(i) - Z(i-j+1));
    end
end

disp('Tabla de diferencias divididas:');
disp(Q);

syms x;
pol = Q(1, 1); 
for i = 2:N
    factor = 1; 
    for j = 1:i-1
        factor = factor * (x - Z(j)); 
    end
    pol = pol + factor*Q(i, i); 
end
 
disp('El polinomio de interpolación es:');
disp(pol);

 
h_interp = double(subs(pol, x, x_interp));

fprintf('El valor aproximado en x = %.4f es: %.9f\n', x_interp, h_interp);


----------------------------------------------------------

% Definir nodos
x = [2, 4, 8, 16]; 
y = log2(x);

% Puntos a interpolar
x_eval = [7.1, 12];
y_exact = log2(x_eval);

% Construcción del trazador cúbico natural
yspline = interp1(x, y, x_eval, 'spline');

% Cálculo del error absoluto
error_abs = abs(y_exact - yspline);

% Mostrar resultados
disp(table(x_eval', yspline', y_exact', error_abs', ...
    'VariableNames', {'x', 'Interpolacion', 'ValorExacto', 'ErrorAbsoluto'}));

% Gráfica
x1 = linspace(2, 16, 100); 
y1 = log2(x1);
yspline_fine = interp1(x, y, x1, 'spline');

plot(x1, y1, 'r', x1, yspline_fine, ':g', x, y, '*')
title('Trazador Cúbico para log2(x)')
legend('Función Original log2(x)', 'Trazador Cúbico', 'Nodos')
grid on;


