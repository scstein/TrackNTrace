function [ y, x ] = GaussianFit( corrfunc )
    % Edit CT, 2020: Switched solver to levenberg-marquardt which is here
    % faster with equal results and disabled output.

    
    [ guess ] = FitInitialGuess( corrfunc, 48 );
    
    radiu = round( 2.5 * guess(2) );
    originy = guess(4);
    originx = guess(3);
    fitdata = corrfunc(originy-radiu:originy+radiu,originx-radiu:originx+radiu);
    [X,Y] = meshgrid(1:1:size(fitdata,2),1:1:size(fitdata,1));
    grid = [X Y];
    
    xoffset = guess(3) - radiu - 1;
    yoffset = guess(4) - radiu - 1;
    guess(3) = radiu + 1;
    guess(4) = radiu + 1;
%     [opt,norm] = lsqcurvefit(@gaussmodel2dwxy, guess, grid, fitdata);
    opt = lsqcurvefit(@gaussmodel2dwxy, guess, grid, fitdata,[],[],optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','Display','off'));

    x = opt(3) + xoffset;
    y = opt(4) + yoffset;
end

