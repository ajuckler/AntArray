%
% Copyright 2016, Antoine JUCKLER
%

function arr = drawCircle(dims, r, t)
% Generate array of which elements are forming a disk
%
% INPUT:
%   dims:   Dimensions of the array
%   r:      External radius of the disk
%   t:      [optional] Internal radius of the disk


if nargin < 3 || t == 0
    t = [];
elseif t >= r
    error('Parameter t should be smaller than r');
end;

arr = zeros(dims);

turnover = dims/2 + 0.5;

for i=1:min(r, floor(turnover))
    maxval = sqrt(r^2-i^2);
    if maxval < 1e-6
        continue;
    end;
    if maxval <= floor(turnover)
        arr(floor(turnover)+i-1, floor(turnover):round(turnover+maxval)) = 1;
    else
        arr(floor(turnover)+i-1, floor(turnover):end) = 1;
    end;
end;

arr(arr'==1)=1;

arr(arr(end:-1:1,:)==1)=1;
arr(arr(:,end:-1:1)==1)=1;

if ~isempty(t)
    arr1 = drawCircle(dims, t);
    arr(arr1==1)=0;
end;

% for i=r:-1:r-t
%     y = round(sqrt(r^2-i^2));
% 
%     yy = turnover-y;
%     xx = turnover+i;
% 
%     if xx <= dims && yy > 0
%         arr(yy, xx)=1;
%     end;
% end
% for i=r:-1:r-t
%     x = round(sqrt(r^2-i^2));
%     
%     xx = turnover+x;
%     yy = turnover-i;
%     
%     if xx <= dims && yy > 0
%         arr(yy, xx)=1;
%     end;
% end;

% arr(arr(end:-1:1,:)==1)=1;
% arr(arr(:,end:-1:1)==1)=1;

end