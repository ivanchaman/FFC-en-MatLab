%Function laguerre=La(x,k,v)

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

  if v ~= 0,
  
    sum = 0;
    for j = 1:v,
    
      aux1 = (-1)^j*coef_binomial(v,j) * x^(v - j);
      aux2 = producto(k - v,j);
      aux1 = aux1 * aux2;
      sum = sum + aux1;
    end
    laguerre = x^v + sum
 else 
    laguerre = 1
end

