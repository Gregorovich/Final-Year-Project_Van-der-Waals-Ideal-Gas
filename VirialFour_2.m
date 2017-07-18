% This is where the doc goes.

%% Input settings
NUM_INTS = 1000                 ;% number of integration intervals to use
NUM_HIT_AND_MISS = 1000000000   ;% number of hit and miss trails to perform 
TREQ_0  = 1.0                   ;% 1st point in temperature interval 
TREQ_N  = 12.0                  ;% final point in temperature interval
NUM_TEMPS = 30                  ;% Number of temperature points to use
IPOT_TYPE = 1                  ;% potential model (1 = LJ, 2 = SplineLJ,4 = MN family etc, 9 = hard disk, 10 = sq well)
RCUT = 10                       ;% LJ Cut
EPS = 1.0                       ;% epsilon
SIGMA = 1.0                     ;% sigma
LAMDA = 1.5                     ;% lamda

    %% More Inputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  N            number of atoms
    %  NCYCL:       number of MOnte Carlo Cycles 
    %  N_POWER      Value of n used in m-n potential for plate-plate interactions
    %  M_POWER      Value of m used in m-n potential for plate-plate interactions
    %  DENSITY      Number density of plate 
    %  SIGMA:       sigma for repulsive ball-plate potential 
    %  EPS:         epsilon for repulsive ball-plate potential 
    %  ISTART:  0 = fresh simulation
    %           1 = start from previous config 
    %           2 = start from previous config nad apply a time reversal map
    %  ISUB         block averaging/screen write interval
    %  MOVWRITE     movie snapshot interval
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch IPOT_TYPE
        case 1    % Cut LJ Potential
            CUTSQ = RCUT*2                  ;
        case 2    % Splined-LJ Potential
            RCUT = 1.737051855              ;
            CUTSQ = RCUT*2                  ;
            RCUT_INNER = 1.24445506         ;
            RCUT_INNER_SQ = RCUT_INNER*2    ;
            SPLINE_R2 = -4.864890083        ;
            SPLINE_R3 = -3.2920028          ;
        case 9    % Hard Sphere Potential
            RCUT = 1                        ;
            CUTSQ = RCUT^2                  ;
        case 10   % Square Well Potential
            EPS = 1                         ;
            SIGMA = 1                       ;
            RCUT = LAMDA*SIGMA              ;
            CUTSQ = RCUT*2                  ;
    end

%% Preallocation
TEMP=zeros(NUM_TEMPS:1);
B2=zeros(NUM_TEMPS:1);



%% And this is where the code starts.
IDUM = -2 ;% Random Number Generator - xran.f90

TW = (TREQ_N-TREQ_0)/(NUM_TEMPS-1) ;% Temp Width
if TW == 0
    TW = 1;
end

switch IPOT_TYPE % Potential type
    case '1'
        fprintf('Cut LJ Potential\r')
    case '2'
        fprintf('Splined-LJ Potential\r')
    case '9'
        fprintf('Hard Sphere Potential\r')
    case '10'
        fprintf('Square Well Potential\r')
end

for i = 1:NUM_TEMPS % loop over all temperatures
    TEMP(i) = TREQ_0 + TW*(i-1);
    TREQ=TEMP(i);
    
    FSW = exp(EPS/TREQ) - 1 ;% no idea. For Square Well...  exp(Epsilon/Temp)
    
    WR = RCUT/NUM_INTS ;% MaxRadius/number_intevals
    HR = WR/2 ;% half? diameter thing?
    XSUM = 0;
    
    
    % CALL Gaus_Legendre_16pt_quad(TREQ,XSUM)
    
    %WEIGHT = [0.152753387, 0.152753387, 0.149172986, 0.149172986, 0.142096109, 0.142096109, 0.131688638, 0.131688638, 0.118194532, 0.118194532, 0.10193012, 0.10193012, 0.083276742, 0.083276742, 0.062672048, 0.062672048, 0.04060143, 0.04060143, 0.017614007, 0.017614007];
    WEIGHT = [0.1527533871307250, 0.1527533871307250, 0.1491729864726030, 0.1491729864726030, 0.1420961093183820, 0.1420961093183820, 0.1316886384491760, 0.1316886384491760, 0.1181945319615180, 0.1181945319615180, 0.1019301198172400, 0.1019301198172400, 0.0832767415767048, 0.0832767415767048, 0.0626720483341091, 0.0626720483341091, 0.0406014298003869, 0.0406014298003869, 0.0176140071391521, 0.0176140071391521];
    
    %ABSCISSA = [-0.076526521, 0.076526521, -0.227785851, 0.227785851, -0.373706089, 0.373706089, -0.510867002, 0.510867002, -0.636053681, 0.636053681, -0.746331906, 0.746331906, -0.839116972, 0.839116972, -0.912234428, 0.912234428, -0.963971927, 0.963971927, -0.993128599, 0.993128599];
    ABSCISSA = [-0.0765265211334973, 0.0765265211334973, -0.2277858511416450, 0.2277858511416450, -0.3737060887154190, 0.3737060887154190, -0.5108670019508270, 0.5108670019508270, -0.6360536807265150, 0.6360536807265150, -0.7463319064601500, 0.7463319064601500, -0.8391169718222180, 0.8391169718222180, -0.9122344282513250, 0.9122344282513250, -0.9639719272779130, 0.9639719272779130, -0.9931285991850940, 0.9931285991850940];
    
    for j = 1:NUM_INTS % loop over all integration intervals
        
        CR = (j-1)*WR + HR  ;% Centre of radius interval
        
        disp(CR)
        
        %DO J = 1,16
        %  X = ABSCISSA(J)
        %  W = WEIGHT(J)
        %  XS = CR + X*HR
        %  RSQ = XS*XS
        %  CALL MAYER_F_FUNCTION(RSQ,TREQ,FIJ)
        %  XSUM = XSUM + FIJ*RSQ*W*HR
        %END DO
        
        for k=1:20
            
            X = ABSCISSA(k) ;
            W = WEIGHT(k)   ;
            
            XS = CR + X*HR  ;
            RSQ = XS*XS     ;
            
            fprintf('XS = %f, RSQ = %f, ',XS,RSQ)
            
            % Mayer F Function: RSQ, TREQ, FIJ
            switch IPOT_TYPE
                case 1    % Cut LJ Potential
                    RRR = 1/RSQ     ;
                    RR6 = RRR^3     ;
                    RR12 = RR6^2    ;
                    UIJ = 4*EPS*(RR12 - RR6)    ;
                    FIJ = exp(-UIJ/TREQ) - 1    ;
                    fprintf('UIJ = %f, FIJ = %f, W = %f, HR = %f\n',UIJ,FIJ,W,HR)
                case 2    % Splined-LJ Potential
                    if RSQ <= RCUT_INNER_SQ
                        RRR = 1/RSQ     ;
                        RR6 = RRR^3     ;
                        RR12 = RR6^2    ;
                        UIJ = 4*EPS*(RR12-RR6)  ;
                    else
                        R = sqrt(RSQ)   ;
                        RDIF = R - RCUT ;
                        RDIF2 = RDIF^2  ;
                        RDIF3 = RDIF2*RDIF ;
                        UIJ = SPLINE_R2*RDIF2 + SPLINE_R3*RDIF3 ;% WHAT IS < SPLINE_R2 > ?
                    end
                    FIJ = exp(-UIJ/TREQ) - 1    ;
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
                    elseif RSQ >= 1 && RSQ < LAMDA*2 % In FORTRAN it says " LAMDA**2 "
                        FIJ = FSW   ;
                    elseif RSQ >= LAMDA*2
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