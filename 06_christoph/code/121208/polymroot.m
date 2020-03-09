function a=polymroot(c)
[rows,cols]=size(c);
n0 = rows/cols - 1;
a = -c(n0*cols+1:rows,:)/c(1:cols,:);
a1 = [];
for i=1:n0-1
    a0 = -c(i*cols+1:(i+1)*cols,:)/c(1:cols,:);
    a1 =[a1 a0];
end;
a = [a1  a];
a=[a; [eye(cols*(n0-1))  zeros(cols*(n0-1),cols)]];
a = eig(a);

