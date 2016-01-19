% Baseline solver using L1 measurments only
function baseline_data = nlib_L1_baseline_solver(measurmentsA, measurmentsB)

BL = 0 ; % number of baselines computed

N = length(measurmentsA) ; % number of receiver A measurements
K = length(measurmentsB) ; % number of receiver B measurements

baseline_data = zeros(5,max(N,K)) ;

L1_state = [] ;

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
        sat_list = [] ;
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
                sat_list(NUMSAT) = recvA_msr{g}.svId ;
                obsA(NUMSAT) = recvA_msr{g}.prMes ;
                obsB(NUMSAT) = recvB_msr{s1}.prMes ;
                Eph(:,NUMSAT) = eph2easy(recvA_msr{g}.s_eph, NUMSAT ) ;
            end
        end
        if NUMSAT>3
            BL = BL + 1 ;
	    measurementTime = (TowA+TowB)/2 ;

	    L1_state = nlib_L1_baseline(L1_state, measurementTime, svIndices, Eph, obsA, obsB) ;

            baseline_data(1, BL) = measurementTime ;
            baseline_data(2:4, BL) = base(1:3) ;
        end
        
    end
end
