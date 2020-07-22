function [ corr ] = CrossCorrelation( image1, image2 )
    
    corr = fftshift( real( ifft2(fft2(image1).*conj(fft2(image2))) ) );
    
end

