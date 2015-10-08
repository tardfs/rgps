function baseline_data = nlib_coarse_baseline_solver(ecef_A, ecef_B)

N = size(ecef_A,2) ;
K = size(ecef_B,2) ;
BS = 0 ;
baseline_data = zeros(5,N) ;
lastk = 1 ;
for n=1:N
    for k=lastk:K
        if round(ecef_A(1,n))==round(ecef_B(1,k))
            lastk=k ;
            BS = BS + 1 ;
            baseline_data(1,BS) = ecef_A(1,n) ;
            baseline_data(2:4,BS) = ecef_B(2:4,k)-ecef_A(2:4,n) ;
        end
    end
end

baseline_data = baseline_data(:,1:BS) ;
