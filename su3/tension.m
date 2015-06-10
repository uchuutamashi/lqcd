load out.txt
k=5; %cutoff
sigma=ones(1,20);
for numrun = 1:20
   x=1:k;
   v=ones(1,k);
   for r=1:k
      y=out( (numrun-1)*8+r , 1:k );
      fm=fit(x',y','exp1');
      v(r)=abs(fm.b);
   end
   fm=fit(x',v','exp1');
   sigma(numrun) = abs(fm.b);
end
