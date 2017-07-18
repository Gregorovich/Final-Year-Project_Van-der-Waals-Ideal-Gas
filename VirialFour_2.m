%% Input settings
n_intgtn = 1000                 ;% number of integration intervals to use
n_monte_carlo_its = 1000000000   ;% number of hit and miss trails to perform 
T_init  = 1.0                   ;% 1st point in temperature interval 
T_fin  = 12.0                  ;% final point in temperature interval
T_n = 30                  ;% Number of temperature points to use
potential_type = 1                  ;% potential model (1 = LJ, 2 = SplineLJ,4 = MN family etc, 9 = hard disk, 10 = sq well)
r_cut = 10                       ;% LJ Cut
epsilon = 1.0                       ;% epsilon
sigma = 1.0                     ;% sigma
lamda = 1.5                     ;% lamda
    
    switch potential_type
        case 1    % Cut LJ Potential
            r_cut2 = r_cut*2                  ;
        case 2    % Splined-LJ Potential
            r_cut = 1.737051855              ;
            r_cut2 = r_cut*2                  ;
            r_cut_in = 1.24445506         ;
            r_cut_in2 = r_cut_in*2    ;
            spline_r2 = -4.864890083        ;
            spline_r3 = -3.2920028          ;
        case 9    % Hard Sphere Potential
            r_cut = 1                        ;
            r_cut2 = r_cut^2                  ;
        case 10   % Square Well Potential
            epsilon = 1                         ;
            sigma = 1                       ;
            r_cut = lamda*sigma              ;
            r_cut2 = r_cut*2                  ;
    end

%% Preallocation
T=zeros(T_n:1);
B2=zeros(T_n:1);



%% And this is where the code starts.
rand_num = -2 ;% Random Number Generator

T_width = (T_fin-T_init)/(T_n-1) ;% Temp Width
if T_width == 0
    T_width = 1;
end

switch potential_type % Potential type
    case '1'
        fprintf('Cut LJ Potential\r')
    case '2'
        fprintf('Splined-LJ Potential\r')
    case '9'
        fprintf('Hard Sphere Potential\r')
    case '10'
        fprintf('Square Well Potential\r')
end

for i = 1:T_n % loop over all temperatures
    T(i) = T_init + T_width*(i-1);
    T_req=T(i);
    
    ForSquareWell = exp(epsilon/T_req) - 1 ;% no idea. For Square Well...  exp(Epsilon/Temp)
    
    MaxR_Treq = r_cut/n_intgtn ;% MaxRadius/number_intevals
    HR = MaxR_Treq/2 ;% half? diameter thing?
    XSUM = 0;
    
    
    % Gaus Legendre quadrature (T_req,XSUM)
    
    WEIGHT = [0.1527533871307250, 0.1527533871307250, 0.1491729864726030, 0.1491729864726030, 0.1420961093183820, 0.1420961093183820, 0.1316886384491760, 0.1316886384491760, 0.1181945319615180, 0.1181945319615180, 0.1019301198172400, 0.1019301198172400, 0.0832767415767048, 0.0832767415767048, 0.0626720483341091, 0.0626720483341091, 0.0406014298003869, 0.0406014298003869, 0.0176140071391521, 0.0176140071391521];
    
    ABSCISSA = [-0.0765265211334973, 0.0765265211334973, -0.2277858511416450, 0.2277858511416450, -0.3737060887154190, 0.3737060887154190, -0.5108670019508270, 0.5108670019508270, -0.6360536807265150, 0.6360536807265150, -0.7463319064601500, 0.7463319064601500, -0.8391169718222180, 0.8391169718222180, -0.9122344282513250, 0.9122344282513250, -0.9639719272779130, 0.9639719272779130, -0.9931285991850940, 0.9931285991850940];
    
    for j = 1:n_intgtn % loop over all integration intervals
        
        CR = (j-1)*MaxR_Treq + HR  ;% Centre of radius interval
        
        disp(CR)
        
        for k=1:20
            
            X = ABSCISSA(k) ;
            W = WEIGHT(k)   ;
            
            XS = CR + X*HR  ;
            RSQ = XS*XS     ;
            
            fprintf('XS = %f, RSQ = %f, ',XS,RSQ)
            
            % Mayer F Function: RSQ, T_req, FIJ
            switch potential_type
                case 1    % Cut LJ Potential
                    RRR = 1/RSQ     ;
                    RR6 = RRR^3     ;
                    RR12 = RR6^2    ;
                    UIJ = 4*epsilon*(RR12 - RR6)    ;
                    FIJ = exp(-UIJ/T_req) - 1    ;
                    fprintf('UIJ = %f, FIJ = %f, W = %f, HR = %f\n',UIJ,FIJ,W,HR)
                case 2    % Splined-LJ Potential
                    if RSQ <= r_cut_in2
                        RRR = 1/RSQ     ;
                        RR6 = RRR^3     ;
                        RR12 = RR6^2    ;
                        UIJ = 4*epsilon*(RR12-RR6)  ;
                    else
                        R = sqrt(RSQ)   ;
                        RDIF = R - r_cut ;
                        RDIF2 = RDIF^2  ;
                        RDIF3 = RDIF2*RDIF ;
                        UIJ = spline_r2*RDIF2 + spline_r3*RDIF3 ;% WHAT IS < spline_r2 > ?
                    end
                    FIJ = exp(-UIJ/T_req) - 1    ;
                    fprintf('UIJ = %f, FIJ = %f, W = %f, HR = %f\n',UIJ,FIJ,W,HR)
                case 9    % Hard Sphere Potential
                    if RSQ < 1
                        FIJ = -1    ;
                    else
                        FIJ = 0     ;
                    end
                    fprintf('UIJ = %f, FIJ = %f, W = %f, HR = %f\n',UIJ,FIJ,W,HR)
                case 10   % Square Well Potential
                    if RSQ < 1
                        FIJ = -1    ;
                    elseif RSQ >= 1 && RSQ < lamda*2 % In FORTRAN it says " lamda**2 "
                        FIJ = ForSquareWell   ;
                    elseif RSQ >= lamda*2
                        FIJ = 0     ;
                    end
                    fprintf('UIJ = %f, FIJ = %f, W = %f, HR = %f\n',UIJ,FIJ,W,HR)
            end
            XSUM = XSUM + FIJ*RSQ*W*HR  ;
            if XSUM > 100
                disp('XSUM > 10')
            end
            
        end
        
        % XSUM = XSUM + FIJ * (CR + X*HR)^2 * HR * W
        
        fprintf('XSUM = %f, XSUM = %f\n',XSUM,FIJ*RSQ*W*HR)
        
        XSUM = -XSUM*2*pi           ;
        % XSUM = B2
    end
    
    fprintf('i = %f, B2 = %f\n',i,XSUM) ;
    
end
