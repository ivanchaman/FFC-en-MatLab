 %  c4 = laguerre(x,ka(1),v1);
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
      
      %      aux2 = producto(ka(1) - v1,j);
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
