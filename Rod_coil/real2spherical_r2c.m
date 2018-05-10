%通过FFT求g(lm)
function qy=real2spherical_r2c(N1,N3C,BNC,N2C,qB0,YPlm)
mycomplex=zeros(1,N3C);
real_t=zeros(1,N3C);
temp2=zeros(N2C,N3C);
qy=zeros(N2C+1,N3C);

for i=1:N2C 
    for j=1:N3C
        real_t(j)=qB0(i,j);
    end
    a=1/(N3C)*fftn(real_t);%固定sita，对g(sita,phi)的phi方向进行傅里叶变换
    mycomplex(1,:)=real(a);%g(sita,phi)傅里叶变换之后的实部
    mycomplex(2,:)=imag(a);%g(sita,phi)傅里叶变换之后的虚部
    for m=1:N3C
        temp2(i,m,1)=mycomplex(1,m);%g_k(sita)实部
	temp2(i,m,2)=mycomplex(2,m);%g_k(sita)虚部
    end
end	

gauss=get_gaussian(N2C);%调用高斯点和权重

%[YPlm,~]=prepare(N1,BNC,N3C,N2C);%调用连带勒让德函数

for l=1:BNC+1
    for m=1:l

        pp2(1)=0;
        pp2(2)=0;
        for  i=1:floor((N2C+1)/2)
             aa=gauss(i,2)* YPlm(((l-1)*l/2+m),i);%对sita进行高斯积分
             pp2(1)=pp2(1)+aa*temp2(i,m,1);
             pp2(2)=pp2(2)+aa*temp2(i,m,2);
        end
        for  i=floor((N2C+1)/2)+1:N2C
             i11=N2C+1-i;
             aa=gauss(i11,2)* YPlm(((l-1)*l/2+m),i);
             pp2(1)=pp2(1)+aa*temp2(i,m,1);
             pp2(2)=pp2(2)+aa*temp2(i,m,2);
        end
 
        pp2(1)=pp2(1)*2*pi;%对phi进行积分
        pp2(2)=pp2(2)*2*pi;
        qy(l,m)=funxushu(pp2(1),pp2(2));%实部与虚部求和，即g(lm)
        if m>1
            qy(l,N3C+2-m)=(-1)^m*conj(qy(l,m));
        end
    end
end

    
        
            
			
			
	
				
			
               
				
           
     
	
