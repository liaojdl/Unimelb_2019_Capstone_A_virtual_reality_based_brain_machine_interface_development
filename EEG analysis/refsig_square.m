function y=refsig_square(f, S, T, H)

% f-- the fundermental frequency
% S-- the sampling rate
% T-- the number of sampling points
% H-- the number of harmonics

for i=1:H
   for j=1:T
    t= j/S;
    sin_square = sin(2*pi*(i*f)*t);
    cos_square = cos(2*pi*(i*f)*t);
    
    if (sin_square > 0)
        y(2*i-1,j)=1;
    elseif (sin_square < 0)
        y(2*i-1,j)=-1;
    else
        y(2*i-1,j) = 0;
    end
    
    if (cos_square > 0)
        y(2*i,j)=1;
    elseif (cos_square < 0)
        y(2*i,j)=-1;
    else
        y(2*i,j)=0;
    end
   end
end
