function x=charroot(z, G)
[p,lp]=size(G);

Mat=eye(p);
for j=1:lp/p
    Mat=z^j*G(:,1+(j-1)*p:j*p);
end
x=det(Mat)^2;
