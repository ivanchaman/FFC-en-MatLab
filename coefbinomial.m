Function coefbinomial=coefbinomial(x,y)

%{Esta funcion calcula el coeficiente binomial de dos numeros.
%                                      x!
%          coef_binomial (x,y)  = -------------
%                                  (x - y)! y!
%}

if y < 0,
   coefbinomial = 0
elseif (y == 0) | ( x == y),
    coefbinomial = 1
    elseif (x >= y) & (y >= 0),
      coef_superior = gamma(x + 1); 
      coef_inferior = gamma(y+1) * gamma ((x-y)+1);
      coefbinomial = coef_superior / coef_inferior
end


