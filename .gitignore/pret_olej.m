function main
    close all;
    y = [1100
         70];
    h = 0.0001;
    t = 0:h:5;
    
    %czêœæ 1
    y=euler_modified(y,h,t);
    
    y(1,length(y))
    y(2,length(y))
    
    plot(t,y(1,:), t, y(2,:));
    legend('Tb','Tw')
    xlabel('t [s]')
    ylabel('Tb,Tw')
    
    %czêœæ 2
    delta_T=[-1500 -1000 -300 -50 -1 1 20 50 200 400 1000 2000];
    approx_h=[];
    spline_h=[];
    for i = 1:length(delta_T)
        approx_h(i)=myapprox(delta_T(i),0);
        spline_h(i)=myspline(delta_T(i),0)
    end
    figure
    plot(delta_T,approx_h,delta_T,spline_h);
    legend('approx','spline')
    xlabel('delta T [C]')
    ylabel('h [w * m^2]')
    
end

function y = euler(y,h,t)
    for i = 1:length(t)-1
        y(:,i+1) = y(:,i) + h*f(y(:,i));
    end
end

function y = euler_modified(y,h,t)
    for i = 1:length(t)-1
        y_next=y(:,i) + h*f(y(:,i));
        y(:,i+1) = y(:,i) + h*f(y_next); 
    end
end

function dy = f(y)
    h=160;
    A=0.0109;
    mb=0.2;
    mw=2.5;
    cb=3.85;
    cw =4.1813;
    Tw=y(2);
    Tb=y(1);
    
    dy=[(Tw-Tb)*h*A / mb*cb
        (Tb-Tw)*h*A / mw*cw];
end

function h=myapprox(Tb,Tw)
   a=calculate_approx_factors;
   x=Tb-Tw;
   h=a(1) + a(2)*x + a(3)*x^2;
end

function h=myspline(Tb,Tw)
   alfa=-0.119;
   beta=0.09;
   n=40;
   interval_len=(1000+1000)/n;
   nodes=-1150:50:1150;
   nodes_values=[]';

   for i=4:length(nodes)-3
       nodes_values(i-3)=myapprox(nodes(i),0);
   end
   
   c=calculate_spline_factors(alfa,beta,n,interval_len,nodes_values');
   deltaT=Tb-Tw;
   h=0;
   for i = 1:length(c)
        h = h + c(i) * fi(i,deltaT,nodes,interval_len);
   end
end

function h=spline_matlab(Tb,Tw)
   nodes=-1000:50:1000;
   nodes_values=[];

   for i=1:length(nodes)
       nodes_values(i)=myapprox(nodes(i),0);
   end
   
   h=spline(nodes,nodes_values,[Tb-Tw]);
end

function a=calculate_approx_factors
   x=[-1500 -1000 -300 -50 -1 1 20 50 200 400 1000 2000]';
   y=[178 176 168 161 160 160 160.2 161 165 168 174 179]';
   m=get_approx_matrix(x);
   mty=m'*y;
   mtm=m'*m;
   a=mtm\mty;
end

function y=fi(i,x,nodes,h)
   y=0;
   if(i<3)
       i=i+2;
   end
   if(nodes(i-2)<= x & x <=nodes(i-1))
       y=(x-nodes(i-2))^3;
   end
   if(nodes(i-1)<= x & x <=nodes(i))
       y=(x-nodes(i-2))^3 - 4*(x-nodes(i-1))^3;
   end
   if(nodes(i)<= x & x <=nodes(i+1))
       y=(nodes(i+2)-x)^3 - 4*(nodes(i+1)-x)^3;
   end
   if(nodes(i+1)<= x & x <=nodes(i+2))
       y=(nodes(i+2)-x)^3;
   end
   y=(1/h^3)*y;
end

function m=get_approx_matrix(x)
    for i = 1:length(x)
        m(i,:)=[1 x(i) x(i)^2];
    end
end

function c=calculate_spline_factors(alfa,beta,n,h,y) 
    m=get_spline_matrix(n);
    temp_c=y\m;
    c(1)=temp_c(2)-h*alfa/3;
    for i=2:length(temp_c)+1
        c(i)=temp_c(i-1);
    end
    c(length(c)+1)=temp_c(length(temp_c)-1)+h*beta/3;
end

function m=get_spline_matrix(n)
    for i = 1:n+1
       for j=1:n+1
           m(i,j)=0;
       end
       if(i==1)
        m(1,1)=4; m(1,2)=2;
       end
       if(i>1 && i<n+1)
        m(i,i-1)=1; m(i,i)=4; m(i,i+1)=1;  
       end
       if(i==n+1)
        m(i,n)=2; m(i,n+1)=4;
       end
    end
end 

