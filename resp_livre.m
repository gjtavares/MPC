function [f]=resp_livre(hor_pred,delta_u,y_medido,g,a_rsalto)

p=hor_pred;

u_pas1=delta_u(1:a_rsalto,1);

for t = 1:p
    
    for e=1:a_rsalto
        
        gk(1,e)=(g(t+e)-g(e));
    end
    
    f(t,1) = y_medido+(gk*u_pas1);
    
end 
       