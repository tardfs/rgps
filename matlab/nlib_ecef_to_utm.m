function utm_data = nlib_ecef_to_utm(ecef_base, ecef_delta)
N = size(ecef_delta,2) ;
utm_data = zeros(4,N) ;
for n=1:N
    utm_data(1,n) = ecef_delta(1,n) ;
    [utm_E,utm_N,utm_U] = topo_enu(ecef_base,ecef_delta(2:4,n)) ;
    utm_data(2:4,n) = [utm_E;utm_N;utm_U] ;
end
