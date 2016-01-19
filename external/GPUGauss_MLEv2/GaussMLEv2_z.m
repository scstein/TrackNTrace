% \file GaussMLEv2_z.m
% \author Peter Relich
% \date June 02, 2011

function[P, CRLB LL] = GaussMLEv2_z(datastack,PSFSigma,iterations,Ax,Ay,Bx,By,gamma,d)

    PSFSigma_y= PSFSigma;
 
    NV = 5;
    Nfits = length(datastack(1,1,:));
    sz = length(datastack(:,1,1));
    PSFSigma_x = PSFSigma; % first PSFSigma input defaults to x direction
    
    P = zeros(Nfits,NV);
    CRLB = zeros(Nfits,NV);

    % Theta is: [x, y, N, bg, z]
    theta = zeros(NV,Nfits);

    maxjump = [1e0, 1e0, 1e2, 2e0, 1e-1];
    g = [1.0, 1.0, 0.5, 1.0, 1.0];
        
    % Initial Values
    [theta(1,:), theta(2,:)] = CenterofMass2D( sz, datastack );
    [Nmax, theta(4,:)] = GaussFMaxMin2D( sz, PSFSigma_x, datastack );
    theta(3,:) = max( 0.0, (Nmax-theta(4,:))*2*pi*PSFSigma_x*PSFSigma_y*sqrt(2.0) );
    theta(5,:) = zeros(1,Nfits);
    
    M = zeros(NV,NV,Nfits);
    Minv = zeros(NV,NV,Nfits);
    
    for kk = 1:iterations % Main iterative loop
        
        % Initialize
        NR_Numerator = zeros(NV,Nfits);
        NR_Denominator = zeros(NV,Nfits);
        %dudt = zeros(NV,Nfits);
        
        for ii = 1:sz
            for jj = 1:sz
                [PSFx, PSFy, dudt, d2udt2] = DerivativeIntGauss2Dz( ii-1, jj-1, theta, PSFSigma_x, PSFSigma_y, Ax, Ay, Bx, By, gamma, d );
                
                model = theta(4,:) + theta(3,:).*PSFx.*PSFy;
                data = permute( datastack(ii,jj,:),[1,3,2] );
                
                % Calculating Remaining Derivatives
                dudt(3,:) = PSFx.*PSFy;
                d2udt2(3,:) = zeros(1,Nfits);
                dudt(4,:) = ones(1,Nfits);
                d2udt2(4,:) = zeros(1,Nfits);
                
                cf = zeros(1,Nfits);
                df = zeros(1,Nfits);
                if model > 10e-3
                    cf = data./model-1;
                    df = data./model.^2;
                end
                cf = min(cf, 10e4);
                df = min(df, 10e4);
                
                for ll = 1:NV
                    NR_Numerator(ll,:) = NR_Numerator(ll,:) + dudt(ll,:).*cf;
                    NR_Denominator(ll,:) = NR_Denominator(ll,:) + d2udt2(ll,:).*cf - (dudt(ll,:).^2).*df;
                end
            end
        end
        
        % The update
        if kk < 3
            for ll = 1:NV
                theta(ll,:) = theta(ll,:) - g(ll)*min( max( NR_Numerator(ll,:)./NR_Denominator(ll,:), -maxjump(ll) ),maxjump(ll) );
            end
        else
            for ll = 1:NV
                theta(ll,:) = theta(ll,:) - min( max( NR_Numerator(ll,:)./NR_Denominator(ll,:), -maxjump(ll) ), maxjump(ll) );
            end
        end
        
        % Other Constraints
        theta(3,:) = max( theta(3,:), 1.0 );
        theta(4,:) = max( theta(4,:), 0.01 );
    end
    
    % Calculating the CRLB and LogLikelihood
    Div = 0.0;
    
    for ii = 1:sz
        for jj = 1:sz
            [PSFx, PSFy, dudt] = DerivativeIntGauss2Dz( ii-1, jj-1, theta, PSFSigma_x, PSFSigma_y, Ax, Ay, Bx, By, gamma, d );
            
            model = theta(4,:) + theta(3,:).*PSFx.*PSFy;
            data = permute(datastack(ii,jj,:),[1,3,2]);
            
            % Calculating Remaining Derivatives
            dudt(3,:) = PSFx.*PSFy;
            dudt(4,:) = ones(1,Nfits);
            
            % Building the Fisher Information Matrix
            for nn = 1:Nfits
                M(:,:,nn) = M(:,:,nn) + dudt(:,nn) * dudt(:,nn)' / model(nn);
                
                % LogLikelyhood
                if model(nn) > 0
                    if data(nn) > 0
                        Div = Div + data(nn)*log(model(nn)) - model(nn) - data(nn)*log(data(nn)) + data(nn);
                    else
                        Div = Div - model(nn);
                    end
                end
            end
        end
    end
    
    % Matrix inverse (CRLB = F^-1) and output assignments
    for nn = 1:Nfits
        Minv(:,:,nn) = inv(M(:,:,nn));
        %Diag = diag(Minv);
    end
    
    % Write to global arrays
    for nn = 1:Nfits
        for kk = 1:NV
            P(nn,kk) = theta(kk,nn);
            CRLB(nn,kk) = Minv(kk,kk,nn);
        end
    end
    LL = Div;
    
end
