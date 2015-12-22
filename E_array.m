function [E_x, E_y, E_z] = E_array(antarray, x, y, z)
% @brief
% Compute the electric field components of a dipole antenna array
%
% Compute the electric field components at frequency f of a dipole antenna
% array with inter-element spacing s at a given position
%
% @param    M       Array containing the amplitude and phase of the current
%                   on the constituting antenna elements
% @param    f       Frequency of operation in Hz
% @param    s       Inter-element spacing in m
% @param    l       Antenna element length in m
% @param    x,y,z   Position wrt the array centre in m
%
% ANTENNAS ARE PLACED IN YZ PLANE

M_in = antarray.M;
s_in = antarray.spacing;
f_in = antarray.freq;
l_in = antarray.el_len;


if size(antarray.M,1) ~= size(antarray.M,2)
    error('M should be a square matrix');
end

E_x = 0;
E_y = 0;
E_z = 0;

if mod(size(M_in,1), 2) == 0
    turnover = size(M_in,1)/2;

    for i=1:size(M_in,1)
        z_el = (turnover-i)*s_in;
        z_el = z_el + s_in/2;
        zz = z - z_el;
        for j=1:size(M_in,2)
            y_el = (j-turnover)*s_in - s_in/2;

            yy = y - y_el;
            [E_x_tmp, E_y_tmp, E_z_tmp] = E_dipole(l_in, M_in(i,j), f_in, x, yy, zz);

            E_x = E_x + E_x_tmp;
            E_y = E_y + E_y_tmp;
            E_z = E_z + E_z_tmp;
        end;
    end;
else
    turnover = (size(M_in,1)-1)/2;
    
    for i=1:size(M_in,1)
        z_el = (turnover-i+1)*s_in;
        zz = z - z_el;
        for j=1:size(M_in,2)
            y_el = (j-turnover-1)*s_in;
            yy = y - y_el;
            
            [E_x_tmp, E_y_tmp, E_z_tmp] = E_dipole(l_in, M_in(i,j), f_in, x, yy, zz);

            E_x = E_x + E_x_tmp;
            E_y = E_y + E_y_tmp;
            E_z = E_z + E_z_tmp;
        end;
    end;
end;

end