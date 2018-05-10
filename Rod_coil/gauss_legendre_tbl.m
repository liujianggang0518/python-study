%对sita方向求高斯点及其对应的高斯权重（剖分节点数>10）
function [x,w]=gauss_legendre_tbl(n,eps)
m=floor((n+1)/2);
t0=(1.0-(1.0-1.0/n)/(8.0*n*n));
t1=1.0/(4.0*n+2.0);
for i=1:m
    pi=3.1415926535;
    x0=cos(pi*(i*4-1)*t1)*t0;
    j=0;
    dw=1;
    dx=1;
    while (abs(dx)>eps||abs(dw))>eps && j<=100
        P_1=1.0;
        P0=x0;
        if n<1024
            k=2;
            while k<=n
                l=1/k;
                P_2=P_1;
                P_1=P0;
                t2=x0*P_1;
                P0=t2+(1-l)*(t2-P_2);
                k=k+1;
            end
        else
            k=2;
            while k<1024
                l=1/k;
                P_2=P_1;
                P_1=P0;
                t2=x0*P_1;
                P0=t2+(1-l)*(t2-P_2);
                k=k+1;
            end
            k=1024;
            while k<=n
                P_2=P_1;
                P_1=P0;
                t2=x0*P_1;
                t3=((k-1)*1.0)/(k*1.0);
                P0=t2+t3*(t2-P_2);
                k=k+1;
            end
        end     
                dpdx=((x0*P0-P_1)*n)/(x0*x0-1.0);
                x1=x0-P0/dpdx;                
                w1=2.0/((1.0-x1*x1)*dpdx*dpdx);
                if j==0
                    w0=2.0/((1.0-x0*x0)*dpdx*dpdx);
                end
                dx=x0-x1;
                dw=w0-w1;
                x0=x1;
                w0=w1;
                j=j+1;
    end
            x(m+1-i)=x1;
            w(m+1-i)=w1;
end
end
