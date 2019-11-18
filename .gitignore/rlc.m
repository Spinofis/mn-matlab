function rlc
    close all;
    y = [0
         0];
    %t = [0];
    h = 0.00001;
    t = 0:h:3;
    
    for i = 1:length(t)-1
        % t(i+1) = t(i) + h;
        % t(i+1) = h*i;
        % Jawny Euler
        %y(:,i+1) = y(:,i) + h*f(t(i), y(:,i)); 
        
        % Niejawny Euler
        y(:,i+1) = y(:,i) + h*f(t(i+1), y(:,i+1)); 
        
        % Ulepszony Euler
        %y_1_2 = y(:,i) + h/2*f(t(i), y(:,i)); 
        %y(:,i+1) = y(:,i) + h*f(t(i)+h/2, y_1_2); 
        
        % Metoda Heuna (tzw. metoda trapez?w)
        %k1 = f(t(i), y(:,i));
        %k2 = f(t(i+1), y(:,i) + h*k1);
        %y(:,i+1) = y(:,i) + h*(k1+k2)/2; 
        
        % Metoda klasyczna Runge-Kutta 4
        %k1 = f(t(i), y(:,i));
        %k2 = f(t(i)+h/2, y(:,i) + h/2*k1);
        %k3 = f(t(i)+h/2, y(:,i) + h/2*k2);
        %k4 = f(t(i+1), y(:,i) + h*k3);
        %y(:,i+1) = y(:,i) + h*(k1+2*k2+2*k3+k4)/6;
        
    end
    plot(t,y(1,:), t, y(2,:))
    legend('uC','iL')
    xlabel('t [s]')
    ylabel('uC,iL')
    iL = y(2,:);
    s = 0;
    for j=1:size(iL,2)
        s = s + h*iL(j) * e( t(j) );
    end
    s
    E = e(t);
    figure 
    plot(t,e(t))
    s = h*E*iL'
end

function v = e(t)
    %e = 1;
    F = 1;
    w = 2*pi*F;
    v = sin(w*t);
end

function dy = f(t, y)
    uC = y(1);
    iL = y(2);
    C = 0.01;
    L = 0.1;
    R = 1;
    
    dy = [  1/C * iL
            1/L*(e(t) - R*iL - uC)
        ];
end