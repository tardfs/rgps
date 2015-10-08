function geodetic_data = nlib_ecef2geodetic(ecef_data)
easyLib = getFullPath('..\\easy') ;
addpath(easyLib) ;

N = size(ecef_data,2) ;
geodetic_data = zeros(5,N) ;
for n=1:N
    [latt, longt] = togeod(6378137,298.257223563, ecef_data(2,n),ecef_data(3,n),ecef_data(4,n) ) ;    
    geodetic_data(1,n) = ecef_data(1,n) ;
    geodetic_data(2,n) = latt ;
    geodetic_data(3,n) = longt ;
end

rmpath(easyLib) ;
