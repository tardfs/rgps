function norm_data = nlib_ecef_norm2(ecef_data)
N = size(ecef_data,2) ;
norm_data = zeros(1,N) ;
for n=1:N
    norm_data(n) = norm(ecef_data(2:4,n)) ;
end
