n=5
%cr�ation de la matrice A et du vecteur b
A = zeros(n,n)
for i=1:n-1
  A(i+1,i)=-1
  A(i,i)=2
  A(i,i+1)=-1
end
A(n,n)=2
b = ones(n,1)
%test r�solution du syst�me Ax=b
[L,U]=lu(A)
y = linsolve (L,b)
x = linsolve (U,y)

  
  