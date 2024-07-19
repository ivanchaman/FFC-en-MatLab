
%+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
%+                                                               +
%+    Titulo de Tesis :                                          +
%+    FACTORES  FRANCK  CONDON  PARA EL  POTENCIAL  DE  MORSE    +
%+                                                               +
%+    .M¢dulo       : FFC.PAS                                    +
%+    .Archivo      : FFC_SIMP.INC                               +
%+    .Autor        : Marisol Abad Ibarra.                       +
%+    .Fecha        : 26 - Mayo - 90                             +
%+    .Objetivo     : Obtiene los factores Franck Condon (FFC)   +
%+                    para el potencial de Morse usando la regla +
%+                    compuesta de Simpson.                      +
%+    .M¢dulos                                                   +
%+      Necesarios  : UTILERIA.PAS                               +
%+                                                               +
%+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+


  lim_a    = 0.4;
  lim_b    = 2.5;
  h        = 0.01;


%{===============================================================}
%{Esta funcion calcula los polinomios generalizados (no asociados)
% de Laguerre por medio de la  siguiente  f¢rmula,  obtenida  del
% art¡culo de Tim Chau (1985) :
%
%  L[k,0] (x) = 1
%                                                  j
%                      v                           __
%  L[k,v] (x) = x^v +  ä (-1) * C(v,j) * x^(v-j) * ³³ (k - v - p)
%                     j=1                         p=1
% donde :
%   a^b      = a elevada a la b-‚sima potencia.
%   C(v,j)   = coeficiente binomial de v en j.         }

Function laguerre=La(x,k,v)

  if v <> 0 then
  
    sum = 0;
    for j = 1:v,
    
      aux1 = (-1)^j*coef_binomial(v,j) * x^(v - j);
      aux2 = producto(k - v,j);
      aux1 = aux1 * aux2;
      sum = sum + aux1;
    end;
    laguerre = x^v + sum;
  end else laguerre = 1;


%{===================================}

%{Esta funci¢n obtiene el valor aproximado de la integral de una
% funci¢n continua. Debe darse un vector   y = f(x) y la ôa y dx.
%                                                        õb
% Se usa la regla compuesta de Simpson para la aproximaci¢n de la
% integral:

%       ôa           (b - a)
%       õb f(x) dx = -------(y0 + 4y1 + 2y2 + ......+ 4yn-1 + yn)
%                       3n
%
% Es importante que el valor de n de entrada a la funci¢n sea un
% n£mero par. Cons£ltese la secci¢n 3.1.1.                       }

Function Simpson=Simp(n, a,b, y)
  s = (y[0] - y[n]) / 2;
  i = 1;
  while i <= (n - 1) do
  begin
    s = s + ( 2 * y[i] + y [i + 1] );
    i=i+2;
  end;
  Simpson = 2 * (b - a) * ( s / ( 3 * n) );



%{===============================================================}
%{Funci¢n que obtiene el valor de la funcion de onda (Ec. 3.7),
% para alg£n estado vibracional.                                 }%
%
Function fun_onda=fonda(r,Nv,v1,i) 
  p = r - re [i];
  x = ka [i] * exp( (-1) * B [i] * p);
  c2 = exp( ((-1) * x) / 2 );
  c3 = x^(ka [i] - 2*v1 - 1) / 2.;
  c4 = laguerre(x,ka [i],v1);
  fun_onda = Nv * c2 * c3 * c4;



%{===============================================================}
%{Procedimiento que genera la funci¢n a integrar: í' * í'' en el
% intervalo especificado, dejando los resultados en un vector y. }
%
Procedure genera_fn_integracion(v1,v2 : byte; var i : integer);
var
  r,aux1,aux2 : extended;
begin
  r = lim_a;  i = -1;
  while r <= lim_b do
  begin
    inc(i);
    aux1 = fun_onda(r,Nv1,v1,1);
    aux2 = fun_onda(r,Nv2,v2,2);
    y [i] = aux1 * aux2;
    r = r + h;
  end;
end;


%{===============================================================}
%{Programa Principal para el c lculo de los FFC por medio de la
% ecuaci¢n 3.7                                                   }

%Procedure FFC_SIMP;


  for v2 = iv2 to limv2 do
  begin
    for es = 1 to esp do
      begin write(' '); if datos then write(f,' '); end;
    write(' ',v2:2,' ³ ');
    if datos then  write(f,' ',v2:2,' ³ ');
    for v1 = iv1 to limv1 do
    begin
      Nv1 = cte_normalizacion (v1,1);
      Nv2 = cte_normalizacion (v2,2);
      genera_fn_integracion(v1,v2,n);
      ffc = simpson(n,lim_a,lim_b,y);
      write(sqr(ffc):6:4,' ');
      if datos then write(f,sqr(ffc):6:4,' ');
    end;
    writeln;
    if datos then writeln(f);
  end;
  tiempo_calculo(datos,f);
  espera_tecla;
end;

%{===============================================================}
%{Calcula a la K de la siguiente manera :    K = we / we_xe      }

Function K (i : byte) : extended;
begin
  K = we[i] / we_xe[i];
end;

%{===============================================================}
%{Calcula a la beta de la siguiente manera :
                                     ______________
%                       beta = cte * û 4 * we_xe * æ             }

Function beta (i : byte) : extended;
const
  cte = 1.21777513710683E-01;
begin
  beta = cte * sqrt(4 * we_xe[i] * mu[i]);
end;

%{===============================================================}
%{Funcion que nos proporciona la constante de normalizacion, pero
% dividida entre n factorial.
% Esto es      N(v) = cte_normalizacion * v!     =     A(v) * v!.
% La  funci¢n  es  calculada  recursivamente  por  medio  de  las
% siguientes f¢rmulas:
%                                  á
%                A(0) = sqrt ( ----------)
%                               â(k - 1)
%
%                             (k - 2v - 3)(k - v - 1)
%        A(v+1) = A(v) sqrt (-------------------------)
%                               (v + 1)(k - 2v - 1)
%
% las cuales son equivalentes a la f¢rmula:
%
%                        á (k - 2v - 1)
%         A(v) = sqrt ( ----------------- )
%                          v! â(k - v)
%
% donde â(k - v) = (k-v-1) â(k-v-1).
%
% F¢rmulas obtenidas del art¡culo de Tim Chau (1985)   }

Function cte_normalizacion (n, i)
  if n = 0 then  cte_normalizacion = sqrt( b[i] / g[i] )
  else
    begin
      v = n - 1;
      num = (ka[i] - (2 * v) - 3) * (ka[i] - v - 1);
      den = (v + 1) * (ka[i] - (2 * v) - 1);
      cte_normalizacion= cte_normalizacion(v,i) * sqrt(num/den);
    end;
end;
%{===============================================================}
%{Funcion que calcula la funcion gamma de x.
% El algoritmo usado para  calcular la  funci¢n  gamma  tiene  una
% exactitud de hasta 10 lugares decimales.  Se  usa la  f¢rmula de
% Stirlings para calcular el logartimo natural de la funci¢n gamma
% a la cual se le saca el exponencial para obtener el  valor  real
% de la funci¢n gamma.
% Nota :                  n! = â(n + 1)                          }
%
Function gamma=gam(x)
  if x < 7 then

    fs = 1.0;
    z = x;
    while z < 7.0 do

      x = z;
      fs = fs * z;
      z = z + 1.0;
    end;
    x = x + 1.0;
    fs = -ln(fs);
  end
  else fs = 0;
  z = sqr(1.0/x);
%  { Uso de la f¢rmula de Stirlings}
  loggamma = fs + (x - 0.5) * ln(x) - x + 0.918938533204673 +
              (((-0.000595238095238 * z + 0.000793650793651
                ) * z - 0.002777777777778
               ) * z + 0.083333333333333
              ) / x;
  gamma = exp(loggamma);
end;

%{===============================================================}
%{Esta funcion calcula el coeficiente binomial de dos numeros.
%                                      x!
%          coef_binomial (x,y)  = -------------
%                                  (x - y)! y!
%}

Function coef_binomial=C(x,y)
  if y < 0 then coef_binomial = 0
  else if (y = 0) or (x = y) then coef_binomial = 1
    else if (x >=y) and (y >= 0) then
      coef_superior = gamma(x + 1); %{x!}
      coef_inferior = gamma(y+1) * gamma ((x-y)+1);
      coef_binomial = coef_superior / coef_inferior;
    end;

%{===============================================================}
%{Esta funcion calcula el siguiente producto:
% Si j = 0    ==>  producto = 1.
% Si j >= 1   ==>
%                    j
%                  ____
%     producto =    ³³  (k - p) = (k - 1)(k - 2) ......(k - j)
%                 p = 1


Function producto (k , j)
  aux = 1;
  for p = 1 :j,
      aux = aux * (k - p);
  producto = aux;
end;

%{===============================================================}
%{Este procedimiento calcula las variables, K, beta y gamma; y
%llama a la rutina que dibuja el marco donde se presentan los FFC}
%
%Procedure calcula_variables_iniciales;

  ka [1] = K (1);
  ka [2] = K (2);
  b [1] = beta (1);
  b [2] = beta (2);
  g [1] = gamma(ka [1] - 1);
  g [2] = gamma(ka [2] - 1);
  gotoxy(1,4);
end;
end.

