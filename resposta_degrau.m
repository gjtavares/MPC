function [g] = resposta_degrau(Num,Den,atraso,N)

    g = zeros(N,1);
    auxDen = [Den zeros(1,N)];
    auxNum = [zeros(1,atraso) Num zeros(1,N)];

    for j=1:N+1
        aux1 = 0;
        aux2 = 0;

        for i=2:j
            aux1 = auxDen(i)*g(j-i+1)+aux1;
        end

        for i=1:j
            aux2 = auxNum(i)+aux2;
        end

        g(j) = aux2-aux1;

    end
    
    g=g(2:N+1,1);

end

