% Isotherms
hold on
k = 1;
subplot(2,2,1)
fprintf('\nPlotting Graphs\n')
fprintf('Isotherms\n')
for i=0:grph_end-1
    
    Tr = i/(grph_end/6*5);% ==i/(grph_end/12*10) : Some Tr's > 1 wanted for curve
    
    if rem(i,update)==0
        fprintf('.')
    end
    
    for j=0:grph_end-1
        
        %Pr_grph(j+1) = -(8*Tr*Vr_grph(j+1)^2 - 9*Vr_grph(j+1) + 3) / (Vr_grph(j+1)^2 - 3*Vr_grph(j+1)^3);
        Pr_grph(j+1) = 8*pr_grph(j+1)*Tr/(3-pr_grph(j+1)) - 3*pr_grph(j+1)^2;
        
        fprintf('\ni=%.0f j=%.0f k=%0.f\tTr = %f\tVr = %f\t Pr = %f\t\n',i,j,k,Tr,Vr_grph(j+1),Pr_grph(j+1))
        
        
        
        
        if k>4
            if k>6
                if     k==7
                    plot(Vr_grph,Pr_grph(i+1),'k.')
                elseif k==8
                    plot(Vr_grph,Pr_grph(i+1),'w.')
                end
            else
                if     k==5
                    plot(Vr_grph,Pr_grph(i+1),'m.')
                elseif k==6
                    plot(Vr_grph,Pr_grph(i+1),'y.')
                end
            end
        else
            if k>2
                if     k==3
                    plot(Vr_grph,Pr_grph(i+1),'r.')
                elseif k==4
                    plot(Vr_grph,Pr_grph(i+1),'c.')
                end
            else
                if     k==1
                    plot(Vr_grph,Pr_grph(i+1),'b.')
                elseif k==2
                    plot(Vr_grph,Pr_grph(i+1),'g.')
                end
            end
        end
        
        
    end
    
        if k == 1
            k = 2;
        elseif k == 2
            k = 3;
        elseif k == 3
            k = 4;
        elseif k == 4
            k = 5;
        elseif k == 5
            k = 6;
        elseif k == 6
            k = 7;
        elseif k == 7
            k = 8;
        elseif k == 8
            k = 1;
        end
    
end
title('Isotherms')
xlabel('Relative Volume')
ylabel('Relative Pressure')
fprintf('\n')
hold off