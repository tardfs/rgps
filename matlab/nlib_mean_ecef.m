function mean_ecef = nlib_mean_ecef(ecef_data)
mean_ecef = zeros(5,1) ;
mean_ecef(2:4) = mean(ecef_data(2:4,:),2) ;