function nlib_num_equal_tow(data1,data2)

N = size(ecef_A,2) ;
K = size(ecef_B,2) ;
BS = 0 ;
lastk = 1 ;
for n=1:N
    for k=lastk:K
        if round(ecef_A(1,n))==round(ecef_B(1,k))
            lastk=k ;
            BS = BS + 1 ;
        end
    end
end
