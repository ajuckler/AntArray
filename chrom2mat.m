function C = chrom2mat(chrom)
% Transforms a chromosome to its matrix representation

chrom_sz = sqrt(length(chrom));
if abs(round(chrom_sz) - chrom_sz) > 10^-6
    error 'Invalid chromosome length';
else
    chrom_sz = round(chrom_sz);
end;

C = zeros(chrom_sz*4);
chrom = repmat(chrom, 2, 1);
chrom = reshape(chrom, chrom_sz*2, chrom_sz);

C(1:end/2, 1:2:end/2) = chrom;
C(1:end/2, 2:2:end/2) = chrom;
C(C(end:-1:1,:)==1) = 1;
C(C(:,end:-1:1)==1) = 1;

end