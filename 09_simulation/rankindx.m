function rankindx=rankindx(a,b) % introduced for Gauss compatibility 
[T,N]=size(a);
aux=[a, (1:T)'];
[aux]=sortrows(aux,b);
aux=[aux, (1:T)'];
[aux]=sortrows(aux,N+1);
rankindx=aux(:,end);