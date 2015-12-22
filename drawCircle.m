function arr = drawCircle(dims, els)

arr = zeros(dims);
r = els^2;

if mod(dims, 2) == 0
    turnover = dims/2;
else
    turnover = (dims+1)/2;
end;

for i=els:-1:0
    y = round(sqrt(r-i^2));

    yy = turnover-y;
    xx = turnover+i;

    if xx <= dims && yy > 0
        arr(yy, xx)=1;
    end;
end
for i=els:-1:0
    x = round(sqrt(r-i^2));
    
    xx = turnover+x;
    yy = turnover-i;
    
    if xx <= dims && yy > 0
        arr(yy, xx)=1;
    end;
end;

end