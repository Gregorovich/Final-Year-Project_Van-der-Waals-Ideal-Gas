% LV coexistence from Pg = Pl and Gg = Gl
hold on
k = 1;
fprintf('LV Coexistence from Pg = Pl and Gg = Gl\n')
subplot(1,2,2)
for i=0:t_end-1
    
    if rem(i,update)==0
        fprintf('.')
    end
    
    T = (i+1)/t_end;
    Vr1 = solns_PG(i*2+1);
    Vr3 = solns_PG(i*2+2);
    
    if isreal(V1) && isreal(V3)
        if k == 1
            plot(Vr1,T,'b.')
            plot(Vr3,T,'b.')
            k = 2;
        elseif k == 2
            plot(Vr1,T,'g.')
            plot(Vr3,T,'g.')
            k = 3;
        elseif k == 3
            plot(Vr1,T,'r.')
            plot(Vr3,T,'r.')
            k = 4;
        elseif k == 4
            plot(Vr1,T,'c.')
            plot(Vr3,T,'c.')
            k = 5;
        elseif k == 5
            plot(Vr1,T,'m.')
            plot(Vr3,T,'m.')
            k = 6;
        elseif k == 6
            plot(Vr1,T,'y.')
            plot(Vr3,T,'y.')
            k = 7;
        elseif k == 7
            plot(Vr1,T,'k.')
            plot(Vr3,T,'k.')
            k = 8;
        elseif k == 8
            plot(Vr1,T,'w.')
            plot(Vr3,T,'w.')
            k = 1;
        end
    end
end
title('Liquid Vapour Coexistence')
xlabel('Relative Volume')
ylabel('Relative Temperature')
fprintf('\n')
hold off