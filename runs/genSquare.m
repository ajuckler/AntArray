function tmp = genSquare(ratio, els, its)
tmp = zeros(els);

half = els/2;
tmp1 = tmp;
tmp2 = tmp;
tmp3 = tmp;

start = round(els/4*(1-ratio));
off = 0;
if mod(round(half*ratio),2) == 1
    disp('PLOP');
    off = 1;
end;
start2 = start+half+1+off+ceil(half*ratio);

k_vec = zeros(1,its);
for i=1:length(k_vec)
   k_vec(i) = round((i-1)*(1+ratio))+1; 
end

for j=1:length(k_vec)
    k = k_vec(j);
    i0 = ceil((j+off)*(1+ratio)/2);
    ie = half+1+off-j;
    
    for i=i0:ie
        ln = start+i+round(k/2)-1;
        rw = start2-ie-ceil(i*ratio);
        tmp(ln,rw) = 1;
    end;
    tmp1 = tmp(end:-1:1,:)';
    tmp2 = tmp(end:-1:1, end:-1:1);
    tmp3 = tmp(:,end:-1:1)';
end;
tmp(tmp1==1)=1;
tmp(tmp2==1)=1;
tmp(tmp3==1)=1;

end