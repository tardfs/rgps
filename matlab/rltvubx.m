function rltvubx()
clc ;
fout = fopen('rltv_out.txt','w+t') ;
recv1File = '..\data\RS_twin04_01.mat' ;
recv2File = '..\data\RS_twin04_02.mat' ;
outFile = '..\data\rltv_out.mat' ;
load(recv1File) ;
recv1_measurments = measurments_queue ;
load(recv2File) ;
recv2_measurments = measurments_queue ;

N = length(recv1_measurments) ;
coarseApos = zeros(N,3) ;
coarseBpos = zeros(N,3) ;
coarseBaseline = zeros(N,4) ;
fineBaseline = zeros(N,4) ;

for idx1=1:length(recv1_measurments)
    recv1_msr_set = recv1_measurments{idx1} ;
    % get rcv1 measurement time
    rcv1Tow = recv1_msr_set{1}.msrTow ;
    % find closest measurment for recv2
    idx2 = -1 ;
    minDelta = 1e100 ; 
    for n=1:length(recv2_measurments)
        recv2_msr_set = recv2_measurments{n} ;
        rcv2Tow = recv2_msr_set{1}.msrTow ;
        if abs(rcv1Tow-rcv2Tow)<minDelta
            minDelta = abs(rcv1Tow-rcv2Tow) ;
            idx2 = n ;
        end
    end
    if idx2>0
        recv2_msr_set = recv2_measurments{idx2} ;
        fprintf(fout, '<DIF:> TOW_A:%f TOW_B:%f\n', recv1_msr_set{1}.msrTow, ...
            recv2_msr_set{1}.msrTow ) ;
        [x_a,y_a,z_a] = pvt_solver(fout,recv1_msr_set) ;
        [x_b,y_b,z_b] = pvt_solver(fout,recv2_msr_set) ;
        
        coarseApos(idx1,:) = [x_a,y_a,z_a] ;
        coarseBpos(idx1,:) = [x_b,y_b,z_b] ;
        
        coarseBaseline(idx1,1) = x_b-x_a ;
        coarseBaseline(idx1,2) = y_b-y_a ;
        coarseBaseline(idx1,3) = z_b-z_a ;
        coarseBaseline(idx1,4) = sqrt((x_b-x_a)^2+ (y_b-y_a)^2+ (z_b-z_a)^2) ;
        
        fprintf(fout, '<COARSE BASELINE:> <%7.2f,%7.2f,%7.2f> <R><%7.2f>\n', ...
            x_b-x_a, y_b-y_a, z_b-z_a, ...
            coarseBaseline(idx1,4) ) ;
        
        % find coarse middle point
        xm = (x_a+x_b)/2 ;
        ym = (y_a+y_b)/2 ;
        zm = (z_a+z_b)/2 ;       
        
        K = zeros(length(recv1_msr_set),4) ;
        R = zeros(length(recv1_msr_set),1) ;
        valid_indices = zeros(length(recv1_msr_set),1) ;
        % proceed with sattelite list
        for g=1:length(recv1_msr_set)
             recv1_msr = recv1_msr_set{g} ;
            % find similar sattelite for second receiver
            idx_msr = -1 ;
            for s=1:length(recv2_msr_set)
                if recv1_msr.svId==recv2_msr_set{s}.svId
                    idx_msr = s ;
                end
            end
            if idx_msr>0
                valid_indices(g) = 1 ;
                recv2_msr = recv2_msr_set{idx_msr} ;
                
                ro_m = (sqrt((recv1_msr.sat_x-x_a)^2+(recv1_msr.sat_y-y_a)^2+(recv1_msr.sat_z-z_a)^2)+ ...
                    sqrt((recv1_msr.sat_x-x_b)^2+(recv1_msr.sat_y-y_b)^2+(recv1_msr.sat_z-z_b)^2))/2 ;
                
                kx = ( (x_b-x_a)/2 - recv1_msr.sat_x)/ro_m ;
                ky = ( (y_b-y_a)/2 - recv1_msr.sat_y)/ro_m ;
                kz = ( (z_b-z_a)/2 - recv1_msr.sat_z)/ro_m ;
                
                K(g,1) = kx ;
                K(g,2) = ky ;
                K(g,3) = kz ;
                K(g,4) = 1 ;
                R(g) = recv2_msr.prMes - recv1_msr.prMes ;
            end
        end
        % get dif solution
        if nnz(valid_indices)>=4
            K = K(valid_indices~=0,:) ;
            R = R(valid_indices~=0) ;
            difSolution = pinv(K)*R ;

            fineBaseline(idx1,1:3) = difSolution(1:3) ;
            fineBaseline(idx1,4) = sqrt(difSolution(1:3).'*difSolution(1:3)) ;


            fprintf(fout, '<FINE BASELINE:> <%7.2f,%7.2f,%7.2f> <R><%7.2f>\n\n', ...
                difSolution(1), difSolution(2), difSolution(3), ...
                 fineBaseline(idx1,4) ) ;
        end
         
    end
end
save(outFile,'coarseBaseline','fineBaseline','coarseApos','coarseBpos') ;
fclose(fout) ;

function [ux,uy,uz] = pvt_solver(fout,measurments_queue)
    numSat = length(measurments_queue) ;
    satpos = zeros(3,numSat) ;
    obs = zeros(numSat,1) ;
    
    %sat_list = [] ;
    %gnssId = 0 ;
    %n = 0 ;
    fprintf(fout, '<PVT SOLVER SAT LIST>\n') ;
    %while length(sat_list)<numSat
    for n=1:length(measurments_queue)
%         m_index = get_actual_measurment(measurments_queue,sat_list,gnssId) ;
%         if m_index<1
%             break ;
%         end
        fprintf(fout, '\t<svId>%2d', measurments_queue{n}.svId ) ;
        fprintf(fout, '\t\t<SAT POS>(%7.2f,%7.2f,%7.2f)\n', measurments_queue{n}.sat_x,measurments_queue{n}.sat_y,measurments_queue{n}.sat_z) ;
        sat_distance = sqrt(measurments_queue{n}.sat_x*measurments_queue{n}.sat_x+measurments_queue{n}.sat_y*measurments_queue{n}.sat_y+measurments_queue{n}.sat_z*measurments_queue{n}.sat_z) ;
        fprintf(fout, '\t\t\t\t<SAT DISTANCE>%7.2f, <PSEUDORANGE>%7.2f\n', sat_distance, measurments_queue{n}.prMes ) ;
        satpos(1,n) = measurments_queue{n}.sat_x ;
        satpos(2,n) = measurments_queue{n}.sat_y ;
        satpos(3,n) = measurments_queue{n}.sat_z ;
        obs(n) = measurments_queue{n}.prMes ;
    end
    
    settings.c    = 299792458 ;    % The speed of light, [m/s]
    settings.useTropCorr = 1 ;     % Use troposphere correction
    [pos, el, az, dop] = leastSquarePos(satpos, obs, settings) ;
    ux = pos(1) ;
    uy = pos(2) ;
    uz = pos(3) ;
    [phi, lambda, h] = cart2geo(pos(1), pos(2), pos(3), 5 ) ;
    fprintf(fout,'<NAVIGATION SOLUTION><ECEF><%f,%f,%f>, <GEO><%f,%f> <HEIGHT><%f>\n\n', ...
        pos(1), pos(2), pos(3), phi, lambda, h ) ;


