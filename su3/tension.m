load out.txt
k=5; %cutoff
n=8; %lattice size
N=20;%num ofconfigurations

sigma=ones(1,N);
for numrun = 1:N
   x=1:k;
   v=ones(1,k);
   for r=1:k
      y=out( (numrun-1)*n+r , 1:k );
      fm=fit(x',y','exp1');
      v(r)=abs(fm.b);
   end
   fm=fit(x',v','exp1');
   sigma(numrun) = abs(fm.b);
end
