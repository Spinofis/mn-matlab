function main
    close all;
    x=-5:5;
    y=[];
    
    for i=1:length(x)
        y(i)=myspline(x(i));
    end
    
    plot(x,y);
end

function h=myspline(x)
   alfa=-10;
   beta=10;
   n=10;
   interval_len=(5+5)/n;
   nodes=-8:1:8;
   nodes_values=[];

   for i=4:length(nodes)-3
       nodes_values(i-3)=nodes(i)^2;
   end
   
   nodes_values(1)=nodes_values(1)+(interval_len*alfa)/3;
   nodes_values(length(nodes_values))=nodes_values(length(nodes_values))-(interval_len*beta)/3;
   
   c=calculate_spline_factors(alfa,beta,n,interval_len,nodes_values');
   h=0;
   for i = 1:length(c)
        h = h + c(i) * fi(i,x,nodes,interval_len);
   end
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
   y=y/h^3;
end

function c=calculate_spline_factors(alfa,beta,n,h,y) 
    m=get_spline_matrix(n);
    temp_c=y\m;
    c(1)=temp_c(2)-(h*alfa/3);
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