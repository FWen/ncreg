function [ X ] = v_2_color_imag(v, h, w)

X = zeros(h,w,3);
for k=1:3 %3 channels
    X(:,:,k) = (reshape(v(:,k),h,w));
end

X = uint8(X);

end

