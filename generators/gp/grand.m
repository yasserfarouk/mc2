function x=grand(Mu,Sigma)
    [V,D]=eig(Sigma);
    A=V*(D.^(1/2));
    x = A * randn(size(A,2),size(Mu,2)) + Mu;
    x=real(x);
end