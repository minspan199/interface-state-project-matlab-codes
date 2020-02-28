% I = imread('tt.tif');
clear J
I = TT;
n = size(I);
% I = I/max(max(I));
J = I;
for k = 1:1:n(2)
    tmp(k) = sum(I(:,k));
    if max(I(:,k))~=0
        t = max(I(:,k))
        J(:,k) = I(:,k)/t;
    end
end

SUM = max(tmp);

for k = 1:1:n(2)
    tmp(k) = sum(I(:,k));
    if sum(I(:,k))~=0
        co = SUM/tmp(k);
        I(:,k) = I(:,k)*co;
    end
end

%     co = SUM/tmp(k)  
%     I(:,k) = sum(I(:,k))*max(I(:,k));


% Max0 = 
imwrite(I,'normalized.tif');
figure;imagesc(I,[0 1.2]);
figure;imagesc(J,[0 1.2]);
colormap jet
colorbar