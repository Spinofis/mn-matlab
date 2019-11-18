function main
    a=1;
    b=3;
    m=4;
    h=(b-a)/m;
    integral=0;
    
    for i=a:2*h:b
        if(i+2*h<=b)
           integral=integral + simpson([i,i+h,i+2*h],h);
        end
    end
    
     integral
end

function y=f(x)
    y=1/(1 + x^2);
end

function y=simpson(x,h)
    y=(h/3) * (f(x(1)) + 4*f(x(2)) + f(x(3)));
end