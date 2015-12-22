function antarray = genArrayStruct(M, f, l, s)

f = f*10^6;
l = l/1000;
s = s/1000;

antarray.M = M;
antarray.freq = f;
antarray.el_len = l;
antarray.spacing = s;

end