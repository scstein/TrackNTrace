function [rp,rs,tp,ts] = Fresnel(w1,n1,n2)

w1 = w1(:).'; n1 = n1(:).'; n2 = n2(:).';

if length(n1)==1 && length(n2)==1
    w2 = sqrt(n2^2-n1^2+w1.^2);
    w2(imag(w2)<0) = conj(w2(imag(w2)<0));
    rp = (w1*n2^2-w2*n1^2)./(w1*n2^2+w2*n1^2);
    rs = (w1-w2)./(w1+w2);
    tp = 2*n1*n2*w1./(w1*n2^2+w2*n1^2);
    ts = 2*w1./(w1+w2);
elseif length(n1)==length(n2)+2
    if length(n1)==2
        [rp,rs,tp,ts] = Fresnel(w1,n1(1),n1(2));
    else
        n = n1; d = [0 n2 0];
        w(1,:) = w1;
        for j = 2:length(n)    
            w(j,:) = sqrt(n(j)^2-n(1)^2+w1.^2);
            w(j,imag(w(j,:))<0) = conj(w(j,imag(w(j,:))<0));        
        end
        
        j = length(n);
        M11 = (w(j,:)./w(j-1,:)*n(j-1)/n(j)+n(j)/n(j-1))/2;
        M12 = (-w(j,:)./w(j-1,:)*n(j-1)/n(j)+n(j)/n(j-1))/2;
        M21 = M12;
        M22 = M11;
        for j=length(n)-1:-1:2
            M11 = exp(-1i*w(j,:)*d(j)).*M11;
            M12 = exp(-1i*w(j,:)*d(j)).*M12;
            M21 = exp(1i*w(j,:)*d(j)).*M21;        
            M22 = exp(1i*w(j,:)*d(j)).*M22;                
            N11 = (w(j,:)./w(j-1,:)*n(j-1)/n(j)+n(j)/n(j-1))/2;
            N12 = (-w(j,:)./w(j-1,:)*n(j-1)/n(j)+n(j)/n(j-1))/2;
            N21 = N12;
            N22 = N11;
            tmp11 = N11.*M11+N12.*M21;
            tmp12 = N11.*M12+N12.*M22;
            tmp21 = N21.*M11+N22.*M21;
            tmp22 = N21.*M12+N22.*M22;
            M11 = tmp11;
            M12 = tmp12;
            M21 = tmp21;
            M22 = tmp22;
        end
        
        rp = M21./M11;
        tp = 1./M11;
        
        j = length(n);
        M11 = (w(j,:)./w(j-1,:)+1)/2;
        M12 = (-w(j,:)./w(j-1,:)+1)/2;
        M21 = M12;
        M22 = M11;
        for j=length(n)-1:-1:2
            M11 = exp(-1i*w(j,:)*d(j)).*M11;
            M12 = exp(-1i*w(j,:)*d(j)).*M12;        
            M21 = exp(1i*w(j,:)*d(j)).*M21;        
            M22 = exp(1i*w(j,:)*d(j)).*M22;                
            N11 = (w(j,:)./w(j-1,:)+1)/2;
            N12 = (-w(j,:)./w(j-1,:)+1)/2;
            N21 = N12;
            N22 = N11;
            tmp11 = N11.*M11+N12.*M21;
            tmp12 = N11.*M12+N12.*M22;
            tmp21 = N21.*M11+N22.*M21;
            tmp22 = N21.*M12+N22.*M22;
            M11 = tmp11;
            M12 = tmp12;
            M21 = tmp21;
            M22 = tmp22;
        end
        rs = M21./M11;
        ts = 1./M11;
    end
else
    error('Wrong input');
end

