function [ guess ] = FitInitialGuess( corrfunc, peakfindsize )

    n = size(corrfunc,1)/2;
       
    [y, x] = ZeroPaddingFFT(corrfunc,peakfindsize,10);
    guess(4)=int16(ceil(y));
    guess(3)=int16(ceil(x));
    
    guess(1)=corrfunc(guess(4),guess(3));
    guess(5) = mean(mean(corrfunc)); % background
    cropsize = 64;
    corr = corrfunc(guess(4)-cropsize/2+1:guess(4)+cropsize/2,guess(3)-cropsize/2+1:guess(3)+cropsize/2);
    [Wy, Wx] = find( ismember( abs((corr./double(guess(1)) - exp(-1))), min(min(abs(corr./double(guess(1)) - exp(-1)))) ) );
    Wx = mod(Wx,size(corr,2));
    guess(2) = mean(( (Wx - cropsize/2 ).^2  + (Wy - cropsize/2 ).^2   ).^(1/2));
    if(guess(2)<1 || guess(2)>10)
        guess(2) = 5;
    end
    
    guess=double(guess);
end

