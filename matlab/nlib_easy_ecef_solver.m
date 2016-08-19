function easy_data = nlib_easy_ecef_solver(measurments)
easyLib = getFullPath('..\\easy') ;
addpath(easyLib) ;

N = length(measurments) ;
easy_data = zeros(5,N) ;
for idx1=1:N
    recv1_msr_set = measurments{idx1} ;
    [ux,uy,uz,t_cor] = easy_pvt_solver( recv1_msr_set) ;
    easy_data(:,idx1) = [recv1_msr_set{1}.msrTow; ux;uy;uz; t_cor] ;
end

rmpath(easyLib) ;

function [ux,uy,uz, time_corr_meters ] = easy_pvt_solver(measurments_queue)
    
numSat = length(measurments_queue) ;
obs = zeros(numSat,1) ;

Eph = zeros(21,length(measurments_queue)) ;
for n=1:length(measurments_queue)
    obs(n) = measurments_queue{n}.prMes ;        
    Eph(:,n) = eph2easy(measurments_queue{n}.s_eph, n ) ; % convert to easy ephemeris
end
    
[easy_pos, el, gdop] = recpo_ls(obs,(1:length(obs)),measurments_queue{1}.msrTow,Eph) ;
%[phi,lambda,h] = togeod(6378137,298.257223563,easy_pos(1),easy_pos(2),easy_pos(3)) ;

ux = easy_pos(1) ;
uy = easy_pos(2) ;
uz = easy_pos(3) ;
time_corr_meters = easy_pos(4) ;
    



