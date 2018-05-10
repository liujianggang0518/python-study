%求区间[-1,1]上的高斯点和高斯权重
function  gauss=get_gaussian(n)

gauss=zeros(n,2);

if rem(n,2)==1
    m=(n+1)/2;
else
    m=n/2;
end
    [x,w]=gauss_legendre(n);
	
    for i=1:m
	    gauss(n+1-i,1)=x(m-i+1);%[0,1]上的高斯点
	    gauss(i,1)=-x(m-i+1);%[-1,0]上的高斯点
            gauss(n+1-i,2)=w(m-i+1);%[0,1]上的高斯权重
            gauss(i,2)=w(m-i+1);  %[-1,0]上的高斯权重    
	end
	
end
