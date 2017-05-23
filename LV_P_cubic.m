% LV Coexistence from solving Pr/Tr/Vr cubic
hold on
k = 1;
fprintf('LV Coexistence: Pr/Tr/Vr Cubic\n')
subplot(2,2,3)
for i=0:t_end-1
    
    Tr = (i+1)/t_end;
    
    for j=1:t_end
        if rem(j,update)==0
            fprintf('.')
        end
        
        V1 = solns_Vr(3*t_end*i+3*(j-1)+1);
        V3 = solns_Vr(3*t_end*i+3*(j-1)+3);
        
        fprintf('V1 = %f\tV3 = %f\tTr = %\n',V1,V3,Tr)
        
        if k == 1
            plot(V1,Tr,'b.')
            plot(V3,Tr,'b.')
            k = 2;
        elseif k == 2
            plot(V1,Tr,'g.')
            plot(V3,Tr,'g.')
            k = 3;
        elseif k == 3
            plot(V1,Tr,'r.')
            plot(V3,Tr,'r.')
            k = 4;
        elseif k == 4
            plot(V1,Tr,'c.')
            plot(V3,Tr,'c.')
            k = 5;
        elseif k == 5
            plot(V1,Tr,'m.')
            plot(V3,Tr,'m.')
            k = 6;
        elseif k == 6
            plot(V1,Tr,'y.')
            plot(V3,Tr,'y.')
            k = 7;
        elseif k == 7
            plot(V1,Tr,'k.')
            plot(V3,Tr,'k.')
            k = 8;
        elseif k == 8
            plot(V1,Tr,'w.')
            plot(V3,Tr,'w.')
            k = 1;
        end
    end
end
title('Isotherm Solutions')
xlabel('Relative Volume')
ylabel('Relative Temperature')
fprintf('\n')
hold off