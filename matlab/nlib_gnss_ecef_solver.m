function ecef_data = nlib_gnss_ecef_solver(measurments)
gnssLib = getFullPath('..\\softgnss') ;
addpath(gnssLib) ;

N = length(measurments) ;
ecef_data = zeros(5,N) ;
for idx1=1:N
    recv1_msr_set = measurments{idx1} ;
    [ux,uy,uz,gdop] = gnss_pvt_solver( recv1_msr_set) ;
    ecef_data(:,idx1) = [recv1_msr_set{1}.msrTow; ux;uy;uz; gdop(1)] ;
end

rmpath(gnssLib) ;


function [ux,uy,uz,gdop, phi,lambda,h, el,az ] = gnss_pvt_solver(measurments_queue)
    
numSat = length(measurments_queue) ;
satpos = zeros(3,numSat) ;
obs = zeros(numSat,1) ;
    
for n=1:length(measurments_queue)
    satpos(1,n) = measurments_queue{n}.sat_x ;
    satpos(2,n) = measurments_queue{n}.sat_y ;
    satpos(3,n) = measurments_queue{n}.sat_z ;
    obs(n) = measurments_queue{n}.prMes + measurments_queue{n}.sat_clk_corr*299792458 ;
end

settings.c    = 299792458 ;    % The speed of light, [m/s]
settings.useTropCorr = 1 ;     % Use troposphere correction
[pos, el, az, dop] = leastSquarePos(satpos, obs, settings) ;

[phi, lambda, h] = cart2geo(pos(1), pos(2), pos(3), 5 ) ;

ux = pos(1) ;
uy = pos(2) ;
uz = pos(3) ;
gdop = dop ;
