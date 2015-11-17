%function baseline_data = nlib_pseudorange_baseline_solver(ecef_A, measurmentsA, measurmentsB, b_t_corr )
function baseline_data = nlib_pseudorange_baseline_solver(varargin)
ecef_A = varargin{1} ;
measurmentsA = varargin{2} ;
measurmentsB = varargin{3} ;
if nargin>3
    t_delta = varargin{4} ;
else
    t_delta = 0 ;
end

easyLib = getFullPath('..\\easy') ;
addpath(easyLib) ;

N = length(measurmentsA) ;
K = length(measurmentsB) ;
BL = 0 ; % number of baselines computed

baseline_data = zeros(5,max(N,K)) ;

for n=1:N
    recvA_msr = measurmentsA{n} ;
    % get rcvA measurement time
    TowA = recvA_msr{1}.msrTow ;
    % find corresponded measurment for recvB
    k = -1 ;    
    for p=1:K
        recvB_msr = measurmentsB{p} ;
        TowB = recvB_msr{1}.msrTow ;
        if round(TowA)==round(TowB)
            k = p ;
            break ;
        end
    end
    if k>0
        % corresponding measurment found for receiver 2
        recvB_msr = measurmentsB{k} ;
        obsA = [] ;
        obsB = [] ;
        Eph = [] ;
        NUMSAT = 0 ;
        % proceed with sattelite list
        for g=1:length(recvA_msr)
            s1 = -1 ;
            for s=1:length(recvB_msr)
                if recvA_msr{g}.svId==recvB_msr{s}.svId
                    % found the same sattelite for receiver B
                    s1 = s ;
                    break ;
                end
            end
            if s1>0
                NUMSAT = NUMSAT + 1 ;
                obsA(NUMSAT) = recvA_msr{g}.prMes ;
                obsB(NUMSAT) = recvB_msr{s1}.prMes ;
                Eph(:,NUMSAT) = eph2easy(recvA_msr{g}.s_eph, NUMSAT ) ;
            end
        end
        if NUMSAT>3
            BL = BL + 1 ;
            [omc, base] = baseline2(ecef_A(2:4),obsA(:),obsB(:),(1:NUMSAT)',TowA,Eph,t_delta) ;

            baseline_data(1, BL) = TowA ;
            baseline_data(2:4, BL) = base(1:3) ;
        end
        
    end
end

baseline_data = baseline_data(:,1:BL) ;

rmpath(easyLib) ;        

