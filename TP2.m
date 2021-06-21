# Formato para trabajar con toda la precision del programa y poder limitarla luego yo.
format long;
clc
clear

# Carga del archivo

datos = load ("data_00.txt");

#************************************** Modelo a **************************************

#************************************ Modelo lineal ***********************************

function resultado = calcular_producto_interno ( vector1, vector2, decimales )

  if length( vector1 ) != length( vector2 )
    resultado = 0;
    return;
  endif
    
  resultado = 0;
    
  for i = 1:length( vector1 )
    aux = vector1 (i) * vector2 (i);
    resultado = resultado + aux;
  endfor

endfunction
function X = eliminacion_gauss (A, b)
  
  [n,m]=size(A);
  C=[A,b];
  for i=1:(n-1)
    mayor=0;
    fila_m=i;
    for j=i:n
      if mayor<abs(C(j,i))
        mayor=abs(C(j,i));
        fila_m=j;
      endif
    endfor
    if fila_m != i 
      for j=1:(n+1)
        aux=C(i,j);
        C(i,j)=C(fila_m,j);
        C(fila_m,j)=aux;
      endfor
    endif
    for j=(i+1):n
      m(j,i)=C(j,i)/C(i,i);
      for l=i:(n+1)
        C(j,l)= C(j,l) - m(j,i)*C(i,l);
      endfor
    endfor
  endfor
  for i=n:-1:1
    suma=0;
     for j=(i+1):n
      suma = suma + C(i,j)*X(j);
     endfor
    X(i)=(C(i,n+1)-suma)/C(i,i);
  endfor  
endfunction
function error = calcular_error_cuadratico ( y_real, y_estimada )
  error = 0;
  for i = 1:rows (y_real)
    error = ( error + (y_real(i) - y_estimada(i))^2 );
  endfor
  
endfunction
# Lleno las funciones base con sus valores
funciones_base = [];
for i = 1:rows(datos)
  x = datos(i,1);
  funciones_base(1,i) = 1;
  funciones_base(2,i) = x;
  funcion_real(i) = datos(i,2);
endfor

# Realizo los productos internos

matriz_A = [];
matriz_B = [];
for i = 1:rows(funciones_base)
  for j = 1:rows(funciones_base)
    matriz_A(i,j) = calcular_producto_interno(funciones_base(i,:),funciones_base(j,:));
  endfor
  matriz_B(i,1) = calcular_producto_interno(funciones_base(i,:), funcion_real);
endfor

solucion = eliminacion_gauss (matriz_A, matriz_B);
x = datos (:,1);
for i = 1:rows (datos);
  y_estimada(i) = (solucion(2)*x(i)) + (solucion(1));
endfor
error = calcular_error_cuadratico ( datos (:,2), y_estimada );
disp(sprintf("El error de calcular con el modelo lineal es %i", error));

# Grafico
x = 0:0.001:100;
y = (solucion(2)*x) + (solucion(1));
plot (x,y);
hold on;
plot (datos(:,1), datos(:, 2),'b.');
hold off;
axis ([20 60 0 100])
print -djpg modelo_lineal.jpg 

#********************************** Modelo cuadratico *********************************

# Lleno las funciones base con sus valores
funciones_base = [];
for i = 1:rows(datos)
  x = datos(i,1);
  funciones_base(1,i) = 1;
  funciones_base(2,i) = x;
  funciones_base(3,i) = x^2;
  funcion_real(i) = datos(i,2);
endfor

# Realizo los productos internos

matriz_A = [];
matriz_B = [];
for i = 1:rows(funciones_base)
  for j = 1:rows(funciones_base)
    matriz_A(i,j) = calcular_producto_interno(funciones_base(i,:),funciones_base(j,:));
  endfor
  matriz_B(i,1) = calcular_producto_interno(funciones_base(i,:), funcion_real);
endfor

solucion = eliminacion_gauss (matriz_A, matriz_B);
x = datos (:,1);
for i = 1:rows (datos);
  y_estimada(i) = (solucion(3)*(x(i)^2)) + (solucion(2)*x(i)) + (solucion(1));
endfor
error = calcular_error_cuadratico ( datos (:,2), y_estimada );
disp(sprintf("El error de calcular con el modelo polinomial (2) es %i", error));

# Grafico
x = 0:0.001:100;
y = solucion(3)*(power(x,2)) + solucion(2)*x + solucion(1);
plot (x,y);
hold on;
plot (datos(:,1), datos(:, 2),'b.');
hold off;
axis ([20 60 0 100])
print -djpg modelo_cuadratico.jpg 

#********************************** Modelo logistico **********************************

# Lleno las funciones base con sus valores
funciones_base = [];
for i = 1:rows(datos)
  x = datos(i,1);
  funciones_base(1,i) = 1;
  funciones_base(2,i) = -x;
  funcion_real(i) = log ( ( 100 / ( datos(i,2) ) ) - 1); # Cambio la unidad de la sorcion del iodo para evitar numeros negativos
endfor

# Realizo los productos internos

matriz_A = [];
matriz_B = [];
for i = 1:rows(funciones_base)
  for j = 1:rows(funciones_base)
    matriz_A(i,j) = calcular_producto_interno(funciones_base(i,:),funciones_base(j,:));
  endfor
  matriz_B(i,1) = calcular_producto_interno(funciones_base(i,:), funcion_real);
endfor

solucion = eliminacion_gauss (matriz_A, matriz_B);
solucion (1) = exp (solucion(1));
solucion (2) = solucion(2);
a = solucion (1);
b = solucion (2);
x = datos (:,1);
for i = 1:rows (datos);
  y_estimada(i) = power((a * power(e, -b*x(i))) + 1, -1) * 100;
endfor
error = calcular_error_cuadratico ( datos(:,2) , y_estimada );
disp(sprintf("El error de calcular con el modelo logistico es %i", error));

# Grafico
x = 0:0.001:100;
y = power((a * power(e, -b*x)) + 1, -1) * 100;

plot (x,y);
hold on;
plot (datos(:,1), datos(:, 2),'b.');
hold off;
axis ([20 60 0 100])

print -djpg modelo_logistico.jpg 

#************************************** Modelo c **************************************

# Declaracion de funciones
function resultado = redondear (numero, decimales)   
  aux = 1;
  for i = 1:decimales
    aux = aux * 10;
  endfor
  
  resultado = round ( numero * aux ) / aux;
endfunction
function [L, U] = doolitle (matriz, decimales)
  n = rows(matriz);
  L = eye(n);
  U = zeros(n,n);
  
  for i=1:n
    for k = i + 1:n
      multi = redondear( matriz (k,i) / matriz (i,i), decimales) ;
      L (k,i) = multi;
      matriz (k,i) = 0;
      for l = i + 1:n
        matriz (k, l) = redondear ( matriz (k, l) - multi*matriz(i,l), decimales );
      endfor
    endfor
    
  endfor
  U = matriz;
endfunction

function [L, U] = crout ( matriz, decimales )
  n = rows(matriz);
  L = eye(n);
  U = zeros(n,n);
  P = eye (n);
  
  for i=1:n
    L (i, i) = redondear ( matriz (i,i) , decimales);
    for l = i:n
        matriz (i,l) = redondear ( matriz (i,l) / L (i,i) , decimales);
    endfor
    for k = i + 1:n
      multi = redondear ( matriz (k,i) / matriz (i,i) , decimales);
      L (k,i) = multi;
      matriz (k,i) = 0;
      for l = i + 1:n
        matriz (k, l) = redondear ( matriz (k, l) - multi*matriz(i,l) , decimales);
      endfor
    endfor
  endfor
  U = matriz;
endfunction

function resultado = resolver_LU ( L, U, b , decimales)
  n = rows(L);
  for i=1:n
    suma=0;
    for p=1:i-1
      suma=suma+redondear(L(i,p)*z(p),decimales);
    endfor
    z(i)=redondear(redondear((b(i)-suma),decimales)/L(i,i), decimales);
  endfor
  for i=n:-1:1
    suma=0;
    for p=(i+1):n
      suma = suma+redondear(U(i,p)*x(p),decimales);
    endfor
    x(i)=redondear(redondear((z(i)-suma),decimales)/U(i,i), decimales);
  endfor
  resultado = x;
  
endfunction
function resultado = calcular_producto_interno_2 ( vector1, vector2, decimales )

  if length( vector1 ) != length( vector2 )
    resultado = 0;
    return;
  endif
    
  resultado = 0;
    
  for i = 1:length( vector1 )
    aux = redondear(vector1 (i) * vector2 (i), decimales);
    resultado = redondear(resultado + aux, decimales);
  endfor

endfunction

function [ solucion, error ] = resolver_con_doolitle ( datos, precision )
  
  funciones_base = [];
  for i = 1:rows(datos)

    x = datos(i,1);
    funciones_base(1,i) = 1;
    funciones_base(2,i) = redondear(x,precision);
    funciones_base(3,i) = redondear(power(x,2),precision);
    funciones_base(4,i) = redondear(power(x,3),precision);
    funcion_real(i) = redondear(datos(i,2),precision);
  endfor

  # Realizo los productos internos

  matriz_A = [];
  matriz_B = [];
  for i = 1:rows(funciones_base)
    for j = 1:rows(funciones_base)
      matriz_A(i,j) = calcular_producto_interno_2(funciones_base(i,:),funciones_base(j,:), precision);
    endfor
    matriz_B(i,1) = calcular_producto_interno_2(funciones_base(i,:), funcion_real, precision);
  endfor

  # Resuelvo por Doolitle
  [L, U] = doolitle (matriz_A, precision);
  solucion = resolver_LU (L, U, matriz_B, precision);
  
  x = datos (:,1);
  for i = 1:rows (datos);
    y_estimada(i) = (solucion(4)*(x(i)^3) + (solucion(3)*(x(i)^2)) + (solucion(2)*x(i)) + (solucion(1)));
  endfor
  error = calcular_error_cuadratico ( datos (:,2), y_estimada );
  
endfunction
function [ solucion, error ] = resolver_con_crout ( datos, precision )
  funciones_base = [];
  for i = 1:rows(datos)

    x = datos(i,1);
    funciones_base(1,i) = 1;
    funciones_base(2,i) = redondear(x,precision);
    funciones_base(3,i) = redondear(power(x,2),precision);
    funciones_base(4,i) = redondear(power(x,3),precision);
    funcion_real(i) = redondear(datos(i,2),precision);
  endfor

  # Realizo los productos internos

  matriz_A = [];
  matriz_B = [];
  for i = 1:rows(funciones_base)
    for j = 1:rows(funciones_base)
      matriz_A(i,j) = calcular_producto_interno_2(funciones_base(i,:),funciones_base(j,:), precision);
    endfor
    matriz_B(i,1) = calcular_producto_interno_2(funciones_base(i,:), funcion_real, precision);
  endfor

  # Resuelvo por Doolitle
  [L, U] = crout (matriz_A, precision);
  solucion = resolver_LU (L, U, matriz_B, precision);
  
  x = datos (:,1);
  for i = 1:rows (datos);
    y_estimada(i) = (solucion(4)*(x(i)^3) + (solucion(3)*(x(i)^2)) + (solucion(2)*x(i)) + (solucion(1)));
  endfor
  error = calcular_error_cuadratico ( datos (:,2), y_estimada );
  
endfunction
# Resuelvo por Doolitle con 3 decimales
[solucion1, error1] = resolver_con_doolitle (datos, 3);
disp(sprintf("El error de calcular con Doolitle con precision 3 es %i", error1));

# Resuelvo por Crout con 3 decimales
[solucion2, error2] = resolver_con_crout (datos, 3);
disp(sprintf("El error de calcular con Crout con precision 3 es %i", error2));

# Resuelvo por Doolitle con 6 decimales
[solucion3, error3] = resolver_con_doolitle (datos, 6);
disp(sprintf("El error de calcular con Doolitle con precision 6 es %i", error3));

# Resuelvo por Crout con 6 decimales
[solucion4, error4] = resolver_con_crout (datos, 6);
disp(sprintf("El error de calcular con Crout con precision 6 es %i", error4));

% Grafico
x = 0:0.001:100;
y1 = (solucion1(4)*(power(x,3))) + (solucion1(3)*(power(x,2)) + (solucion1(2)*x) + (solucion1(1)));
y2 = (solucion2(4)*(power(x,3))) + (solucion2(3)*(power(x,2)) + (solucion2(2)*x) + (solucion2(1)));
y3 = (solucion3(4)*(power(x,3))) + (solucion3(3)*(power(x,2)) + (solucion3(2)*x) + (solucion3(1)));
y4 = (solucion4(4)*(power(x,3))) + (solucion4(3)*(power(x,2)) + (solucion4(2)*x) + (solucion4(1)));
plot (x,y1);
hold on;
plot (x,y2);
plot (x,y3);
plot (x,y4);
plot (datos(:,1), datos(:, 2),'b.');
hold off;
legend('Doolitle 3 decimales','Crout 3 decimales','Doolitle 6 decimales','Crout 6 decimales')
axis ([20 60 0 100])
print -djpg modelo_polinomial_3.jpg 

# *************************************************************************************************