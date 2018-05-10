function sum=spherical_harmonic(N1,BNC,l,m,x)
sum=0;
legendre=zeros(N1,N1);
spherical=zeros(N1,N1,N1);
legendre(1,1)=1;
legendre(2,1)=0;
legendre(2,2)=1;
for n=2:N1-1
    for j=1:n+1
        if j>1
            legendre(n+1,j)=(2*(n-1)+1)*1.0/(n)*legendre(n,j-1)-(n-1)*1.0/n*legendre(n-1,j);
        else
            legendre(n+1,j)=-(n-1)*1.0/n*legendre(n-1,j);
        end
    end
end

% for a=1:BNC
 %   for b=1:N1
        spherical(1:BNC,1,1:N1)=legendre(1:BNC,1:N1);
  %  end
%end
for g=1:BNC
    for k=2:g
        for i=1:N1-1
            spherical(g,k,i)=i*spherical(g,k-1,i+1);
        end
    end
end


       
if m<=l
    if m==l
       sum=spherical(l,m,1);
       sum=sum*(1-x*x)^((m-1)/2.0)*(-1)^(m-1)*((2*l-1)/(4*pi)*factorial(l-m)*1.0/factorial(l+m-2))^0.5;
     else
      sum=spherical(l,m,l-m+1);
      i=l-m;
      while i>0
           sum=sum*x;
	       sum=sum+spherical(l,m,i);
           i=i-1;
      end
          sum=sum*(1-x*x)^((m-1)/2.0)*(-1)^(m-1)*((2*l-1)/(4*pi)*factorial(l-m)*1.0/factorial(l+m-2))^0.5;
          %sum 连带勒让德函数
    end
end
end
 