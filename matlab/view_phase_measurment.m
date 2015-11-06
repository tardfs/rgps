%function view_phase_measurment

recv1File = '..\data\RS_matv_1400mm_680mm_01.mat' ;
load(recv1File) ;
recv1_measurments = measurments_queue ;
svId = 4 ;
cpMes = zeros(length(recv1_measurments),1) ;
for n=1:length(recv1_measurments)
    ms = recv1_measurments{n} ;
    for k=1:length(ms)
        if ms{k}.svId==svId
            cpMes(n) = ms{k}.cpMes ;
        end
    end
end