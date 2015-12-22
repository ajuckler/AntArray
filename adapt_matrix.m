function antarray = adapt_matrix(antarray, M, x, y, z, yc, zc)
% Adapt amplitude and phase from the input matrix to focus its beam at the 
% given position wrt to the given centre

if nargin < 6
    yc = 0;
    zc = 0;
end;

if size(M,1) ~= size(antarray.M,1) || size(M,2) ~= size(antarray.M, 2)
    error('Matrix sizes do not match');
end;

c0 = 299792458;     % Speed of light in free space
lambda = c0/antarray.freq;

s = antarray.spacing;

x = x/1000;
y = y/1000;
z = z/1000;
yc = yc/1000;
zc = zc/1000;

rho_g = sqrt(x^2+y^2+z^2);
rho = rho_g;

L = length(antarray.M)*antarray.spacing;
gamma = rho_g/2/L^2*lambda;
disp(['Gamma: ' mat2str(gamma)]);

if mod(size(M,1), 2) == 0
    turnover = size(M,1)/2;

    for i=1:size(M,1)
        z_el = (turnover-i)*s + s/2 + zc;
        for j=1:size(M,2)
            y_el = (j-turnover)*s - s/2 + yc;
            
            if M(i,j) == 1
                phi = 2*pi/lambda * (sqrt(x^2+(y_el-y)^2+(z_el-z)^2) - rho);
                antarray.M(i,j) = exp(-1j*phi);
            end;
        end;
    end;
else
    turnover = (size(M,1)-1)/2;
    
    for i=1:size(M,1)
        z_el = (turnover-i+1)*s + zc;
        for j=1:size(M,2)
            y_el = (j-turnover-1)*s + yc;
            
            if M(i,j) == 1
                phi = 2*pi/lambda * (sqrt(x^2+(y_el-y)^2+(z_el-z)^2) - rho);
                antarray.M(i,j) = exp(-1j*phi);
            end;
        end;
    end;
end;

end