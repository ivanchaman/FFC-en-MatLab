%Programa para calcular los Factores Franck-Condion de Morse
%Usando la regla compuesta de Simpson
%Realizo: Lourdes Sandoval
%Fecha: 14 de Abril del 2002
%
clear; clc;
format long e
%Primer estado
we(1)= 1675.355
wexe(1)=13.433
re(1)=1.21252
%Segundo estado
we(2)=1411.210
wexe(2)=12.921
re(2)=1.2864
mu=7.500053859
v1=0
v2=0
n=210
%Inicializacion de variables
ka (1) = we(1)/wexe(1)
ka (2) = we(2)/wexe(2)
cte = 1.21777513710683E-01;
b(1) = cte * sqrt(4 * wexe(1)* mu)
b(2) = cte * sqrt(4 * wexe(2)* mu)
g(1) = gamma(ka(1) - 1)
g(2) = gamma(ka(2) - 1)
lim_a=0.4;
lim_b=3.5;
h=(lim_b-lim_a)/n;
%Generacion de Funcion de Onda para cada estado
%Funonda(1) del primer estado
%constante de normalizacion           
cte_normalizacion=sqrt(b(1)/g(1));
for i=0:v1-1,
    num=(ka(1)-2*i-3.)*(ka(1)-i-1);
    den=(i+1)*(ka(1)-2*i-1);
    cte_normalizacion=cte_normalizacion*sqrt(num/den);
end
Nv1=cte_normalizacion;           
%funonda(2) del segundo estado
%constante de normalizacion           
cte_normalizacion=sqrt(b(2)/g(2));
for i=0:v2-1,
    num=(ka(2)-2*i-3.)*(ka(2)-i-1);
    den=(i+1)*(ka(2)-2*i-1);
    cte_normalizacion=cte_normalizacion*sqrt(num/den);
end
Nv2=cte_normalizacion;           
%Generacion de la funcion a integrar que es funonda(v1)*funonda(2)
r=lim_a;
%ri(0)=lim_a;
ri=zeros();
y=zeros();
for i=1:210, %Numero de Nodos
    p = r - re(1);
    x = ka(1) * exp( (-1) * b(1) * p);
    c2 = exp( ((-1) * x) / 2 );
    c3 = (ka(1) - 2*v1 - 1) / 2;
    c3= x^c3;  
    %c4 = laguerre(x,ka(1),v1);
    if v1 ~= 0,
        sum = 0;
        for j = 1:v1,
        %calculo del coeficiente binomial  
            if ( v1 == j),
                coef_binomial = 1;
            else               %(v1 >= j) & (j >= 0),
                coef_superior = gamma(v1 + 1); 
                coef_inferior = gamma(j+1) * gamma ((v1-j)+1);
                coef_binomial = coef_superior / coef_inferior;
            end    
            aux1 = (-1)^j*coef_binomial * x^(v1 - j);         
            %aux2 = producto(ka(1) - v1,j);
            aux = 1;
            for p = 1:j,
                aux = aux * (ka(1) - v1 - p);
            end;
            aux2=aux;      
            aux1 = aux1 * aux2;
            sum = sum + aux1;
        end
        laguerre = x^v1 + sum;
    else 
        laguerre = 1;
    end
    c4=laguerre;
    %fun_onda1 = ((-1.0)^v1)*Nv1 * c2 * c3 * c4;
    fun_onda1 = Nv1 * c2 * c3 * c4;
    %generacion de la funcion de onda del segundo estado
    p = r - re (2);
    x = ka(2) * exp( (-1) * b(2) * p);
    c2 = exp( ((-1) * x) / 2 );
    c3 = x^((ka(2) - 2*v2 - 1) / 2.);  
    % c4 = laguerre(x,ka(1),v1);
    if v2 ~= 0,
        sum = 0;    
        for j = 1:v2,
           %calculo del coeficiente binomial  
            if ( v2 == j),
                coef_binomial = 1;
            else               %(v2 >= j) & (j >= 0),
                coef_superior = gamma(v2 + 1); 
                coef_inferior = gamma(j+1) * gamma ((v2-j)+1);
                coef_binomial = coef_superior / coef_inferior;
            end    
            aux1 = (-1)^j*coef_binomial * x^(v2 - j);      
            %aux2 = producto(ka(2) - v2,j);
            aux = 1;
            for p = 1:j,
                aux = aux * (ka(2) - v2 - p);
            end;
            aux2=aux;      
            aux1 = aux1 * aux2;
            sum = sum + aux1;
        end
        laguerre = x^v2 + sum;
    else 
        laguerre = 1;
    end
    c4=laguerre;
    %fun_onda2 = ((-1.0)^v2)*Nv2 * c2 * c3 * c4;
    fun_onda2 = Nv2 * c2 * c3 * c4;
    y(i)=fun_onda1*fun_onda2;
    r=r+h;
    ri(i)=r;
end 
plot(ri,y,'g');
y'
%Aplicando la Regla compuesta de Simpson 
s = (y(1) - y(210)) / 2;
i = 1;
while i <= 209 %Numero de Nodos
    s = s + ( 2 * y(i) + y (i + 1) );
    i=i+2;
end;
Simpson = 2 * (lim_b - lim_a) * ( s / ( 3 * 210) )
FFC=Simpson^2

           