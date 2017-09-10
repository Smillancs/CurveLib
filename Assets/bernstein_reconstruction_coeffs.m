function F = bernstein_reconstruction_coeffs(deg, cont)

F = zeros(2*cont, deg+1);

n = deg;

for i=1:cont
  syms k t;
  if i==1
    b = nchoosek(n,k)*t^k*(1-t)^(n-k);
  else
    b = diff(b, t);
  end
  for k=0:deg
    F(i,k+1) = subs(limit(subs(b), t, 0));
  end
end

for i=1:cont
  syms k t;
  if i==1
    b = nchoosek(n,k)*t^k*(1-t)^(n-k);
  else
    b = diff(b, t);
  end
  for k=0:deg
    F(cont+i,k+1) = subs(limit(subs(b), t, 1));
  end
end

F(abs(F)<10^-10)=0;

pinv(F)

F*pinv(F)

end
