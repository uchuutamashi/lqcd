load out.txt
sigma=ones(1,20);
for numrun = 1:20
   x=1:5;
   v=ones(1,5);
   for r=1:5
      y=out( (numrun-1)*8+r , 1:5 );
      fm=fit(x',y','exp1');
      v(r)=abs(fm.b);
   end
   fm=fit(x',v','exp1');
   sigma(numrun) = abs(fm.b);
end
