# Formato para trabajar con toda la precision del programa y poder limitarla luego yo.

format long;

# Carga del archivo

datos = load ("data_00.txt");

# Declaracion de funciones
function resultado = redondear (numero, decimales) 
  str = mat2str(numero, decimales);
  resultado=eval(str);
    
endfunction
function [L, U] = doolitle (A, decimales)
  L=[];
  U=[];
  [n,m] = size (A);
  for k =1:n;
    L(k,k)=1; %la diagonal de la L tiene que ser 1
    suma1 = 0;
    for p=1:k-1 % hallar valores que van en la diagonal de la u 
      suma1=redondear(suma1+redondear(L(k,p)*U(p,k), decimales),decimales);
    endfor
    U(k,k)=redondear((A(k,k)-suma1), decimales);
    for i=k+1:n 
      suma2=0;
      for p=1:(k-1) % hallar valores que van diferentes de  la diagonal de L
        suma2=redondear(suma2+redondear(L(i,p)*U(p,k), decimales),decimales);
      endfor
      L(i,k)=redondear((redondear(A(i,k)-suma2,decimales))/U(k,k),decimales);
    endfor
    for j=k+1:n  % hallar valores que van diferentes de  la diagonal de U
      suma3=0;
      for p=1:k-1
        suma3=redondear(suma3+redondear(L(k,p)*U(p,j), decimales),decimales);
      endfor
    U(k,j)=redondear((A(k,j)-suma3)/L(k,k), decimales);
    endfor
  endfor
endfunction

function [L, U] = crout ( A, decimales )
  
  [n,m]=size(A);
  for k=1:n
    %La instrucción iterativa for permite repetir estamentos a un
    %numero específico de veces  
    U(k,k)=1; %princio del metodo
    suma=0;
    for p=1:k-1 %diagonal
      suma=redondear(suma+redondear(L(k,p)*U(p,k), decimales), decimales);
    endfor
    L(k,k)=redondear((A(k,k)-suma), decimales); 
        
    for i=k+1:n %no hacen parte de la diagonal
      suma=0;
      for r=1:k-1 %matriz u
        suma=redondear(suma+redondear(L(i,r)*U(r,k), decimales), decimales);
      endfor
      L(i,k)=redondear((A(i,k)-suma), decimales); %obtencion de la matriz L
    endfor
    for j=k+1:n
      suma=0;
        for s=1:k-1
          suma=redondear(suma+redondear(L(k,s)*U(s,j), decimales), decimales);
        endfor
      U(k,j)=redondear(redondear((A(k,j)-suma), decimales)/L(k,k), decimales); %obtencion de la matriz U
      endfor
    endfor
endfunction

function resultado = resolver_LU ( L, U, b , decimales)
  n = rows(L);
  for i=1:n %sustitución progresiva
    suma=0;
    for p=1:i-1
      suma=suma+redondear(L(i,p)*z(p),decimales);
    endfor
    z(i)=redondear(redondear((b(i)-suma),decimales)/L(i,i), decimales); %obtencion del vector Z
  endfor
  for i=n:-1:1
    suma=0;
    for p=(i+1):n
      suma = suma+redondear(U(i,p)*x(p),decimales);
    endfor
    x(i)=redondear(redondear((z(i)-suma),decimales)/U(i,i), decimales); % solcion, calculos de las variables
  endfor
  fprintf('\n Matriz L:\n')
  disp(L)
  fprintf('\n Matriz U:\n')
  disp(U)
  fprintf('\n El vector Z:\n')
  disp(z)
  fprintf('\n\nLa solucion de X1 hasta Xn es:\n');
  %a continuacion de utiliza una instruccion for, para mostrar el usuario, 
  %los resultados de una manera mas ordenada
  for i=1:n
    xi=x(1,i);
    fprintf('\nX%g=',i)
    disp(xi);
  endfor
  resultado = x;
  
endfunction
function resultado = calcular_producto_interno ( vector1, vector2, decimales )

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


function solucion = resolver_con_doolitle ( datos, precision )
  
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
      matriz_A(i,j) = calcular_producto_interno(funciones_base(i,:),funciones_base(j,:), precision);
    endfor
    matriz_B(i,1) = calcular_producto_interno(funciones_base(i,:), funcion_real, precision);
    disp(matriz_B(i,1));
  endfor
  disp (matriz_A);

  # Resuelvo por Doolitle
  [L, U] = doolitle (matriz_A, precision);
  solucion = resolver_LU (L, U, matriz_B, precision);
  
endfunction
function solucion = resolver_con_crout ( datos, precision )
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
      matriz_A(i,j) = calcular_producto_interno(funciones_base(i,:),funciones_base(j,:), precision);
    endfor
    matriz_B(i,1) = calcular_producto_interno(funciones_base(i,:), funcion_real, precision);
    disp(matriz_B(i,1));
  endfor
  disp (matriz_A);

  # Resuelvo por Doolitle
  [L, U] = crout (matriz_A, precision);
  solucion = resolver_LU (L, U, matriz_B, precision);
  
endfunction


# Resuelvo por Doolitle con 3 decimales
solucion1 = resolver_con_doolitle (datos, 3);

# Resuelvo por Crout con 3 decimales
solucion2 = resolver_con_crout (datos, 3);

# Resuelvo por Doolitle con 6 decimales
solucion3 = resolver_con_doolitle (datos, 6);

# Resuelvo por Crout con 6 decimales
solucion4 = resolver_con_crout (datos, 6);

# Resuelvo por Doolitle con 20 decimales
solucion5 = resolver_con_doolitle (datos, 20);

# Resuelvo por Crout con 20 decimales
solucion6 = resolver_con_crout (datos, 20);

% Rango de cÃ¡lculo de la variable x
x = 0:0.001:100;
y1 = (solucion1(4)*(power(x,3))) + (solucion1(3)*(power(x,2)) + (solucion1(2)*x) + (solucion1(1)));
y2 = (solucion2(4)*(power(x,3))) + (solucion2(3)*(power(x,2)) + (solucion2(2)*x) + (solucion2(1)));
y3 = (solucion3(4)*(power(x,3))) + (solucion3(3)*(power(x,2)) + (solucion3(2)*x) + (solucion3(1)));
y4 = (solucion4(4)*(power(x,3))) + (solucion4(3)*(power(x,2)) + (solucion4(2)*x) + (solucion4(1)));
y5 = (solucion5(4)*(power(x,3))) + (solucion5(3)*(power(x,2)) + (solucion5(2)*x) + (solucion5(1)));
y6 = (solucion6(4)*(power(x,3))) + (solucion6(3)*(power(x,2)) + (solucion6(2)*x) + (solucion6(1)));
plot (x,y1);
hold on;
plot (x,y2);
plot (x,y3);
plot (x,y4);
plot (x,y5);
plot (x,y6);
plot (datos(:,1), datos(:, 2),'b.');
hold off;
legend('Doolitle 3 decimales','Crout 3 decimales','Doolitle 6 decimales','Crout 6 decimales','Doolitle 20 decimales','Crout 20 decimales')
axis ([0 60 0 100])
print -djpg archivo_01.jpg 