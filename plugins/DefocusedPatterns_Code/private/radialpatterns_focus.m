function [err, bim, cim, sim, xc, yc, bc, cc, sc, len, imm, imm_c, mask, molim] = radialpatterns_focus(tag, pattern, pixel, theta, phi, NA, n0, n, n1, d0, d, d1, lamex, focpos, pic)
% radialpatterns_focus computes patterns for the givel pixel size and the
% excitation wavelength and then identifies molecules by calling
% FindPatternSep.m
    if nargin<14 || isempty(pic)
        pic=0;
    end

    if nargin<13 || isempty(focpos)
        focpos=0;
    end

    close all
    nn = round(0.6/pixel);
    resolution = [lamex/0.02 lamex/0.001];
    rhofield = [-lamex/resolution(1)/2 nn*pixel*1.1];
    zfield = [0 0.01];
    fd = 3e3;

    over = inf;
    maxm = 3;
    atf  = [];
    ring = 'cos(psi).*rad'; 
    % subplot(2,1,1)
    % FocusImage2D(rho,z,cat(3,fxc,fxs))
    % axis normal
    % colorbar
    % subplot(2,1,2)
    % FocusImage2D(rho,z,cat(3,fzc,fzs))
    % colorbar
    % axis normal
    % 

    [xx,yy] = meshgrid(-nn:nn,-nn:nn);
    rr = pixel*sqrt(xx.^2+yy.^2);
    psi = angle(xx + 1i*yy);
    mask = zeros(2*nn+1,2*nn+1,numel(theta)*numel(focpos));
    c=0;
    
    for i=1:numel(focpos)
        [fxc{i}, fxs{i}, fyc{i}, fys{i}, fzc{i}, fzs{i}, rho{i}, ~] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos(i), atf, resolution, ring, maxm);            %#ok<*AGROW>
        [fxc2{i}, fxs2{i}, fyc2{i}, fys2{i}, fzc2{i}, fzs2{i}] = RotateEMField(fxc{i}, fxs{i}, fyc{i}, fys{i}, fzc{i}, fzs{i}, pi/2);
        fxc{i} = fxc{i} + fxc2{i};
        fxs{i} = fxs{i} + fxs2{i};
        fyc{i} = fyc{i} + fyc2{i};
        fys{i} = fys{i} + fys2{i};
        fzc{i} = fzc{i} + fzc2{i};
        fzs{i} = fzs{i} + fzs2{i};
        

        for k=1:numel(theta)
            c=c+1;
                 mask(:,:,c) = interp1(rho{i}(:,1),(fxc{i}(:,1,1)*sin(theta(k))*cos(phi(k))+fyc{i}(:,1,1)*sin(theta(k))*sin(phi(k))+fzc{i}(:,1,1)*cos(theta(k))) ,rr,'cubic',0);
            for j=1:maxm
                 mask(:,:,c) = mask(:,:,c) + interp1(rho{i}(:,1),(fxc{i}(:,1,j+1)*sin(theta(k))*cos(phi(k))+fyc{i}(:,1,j+1)*sin(theta(k))*sin(phi(k))+fzc{i}(:,1,j+1)*cos(theta(k))) ,rr,'cubic',0).*cos(j*psi) + ...
                 interp1(rho{i}(:,1),(fxs{i}(:,1,j)*sin(theta(k))*cos(phi(k))+fys{i}(:,1,j)*sin(theta(k))*sin(phi(k))+fzs{i}(:,1,j)*cos(theta(k))) ,rr,'cubic',0).*sin(j*psi);
            end
        end

    end


    mask = abs(mask).^2;
%     bck = Disk(nn);

%     if pic==1;
%         col = ceil(sqrt(size(mask,3)));
%         wdth = size(mask,2);
%         hght = size(mask,1);
%         immask = zeros(ceil(size(mask,3)/col)*hght,col*wdth);
%         for j=1:size(mask,3)
%             immask(fix((j-1)/col)*hght+1:(fix((j-1)/col)+1)*hght,mod(j-1,col)*wdth+1:(mod(j-1,col)+1)*wdth) = mask(:,:,j)./max(max(mask(:,:,j)));
%         end
%     mim(immask)
%     end

    if strcmpi(pattern,'radial')
        tresh = 1;
    else
        tresh = 0.8;
    end

   [err, bim, cim, sim, xc, yc, bc, cc, sc, len, imm, imm_c] = FindPatternSep(tag,mask,mask>0,[],[],tresh,[],0);

    if pic==1
        figure
        CombineImages(cat(3,tag,imm_c),1,2)
    end

%     if pic==1
%         figure
%         set(gcf,'name',['pixel size = ',num2str(pixel*1e3),' nm'],'NumberTitle','off')
    if ~isempty(len)
        nx = 2*round(sqrt(numel(len)/2)); 
        ny = ceil(numel(len)/nx);
        molim = zeros(2*(2*nn+1)*ny,(2*nn+1)*nx); 
        tag_ext=zeros(size(tag,1)+2.*nn+1,size(tag,2)+2.*nn+1);
        tag_ext(nn+1:end-nn-1,nn+1:end-nn-1)=tag;
        imm_c_ext=zeros(size(imm_c,1)+2.*nn+1,size(imm_c,2)+2.*nn+1);
        imm_c_ext(nn+1:end-nn-1,nn+1:end-nn-1)=imm_c;
        for jx=1:nx
            for jy=1:ny
                if (jy-1)*nx+jx<=numel(len)
                    j = (jy-1)*nx+jx;
                    molim((2*jy-2)*(2*nn+1)+(1:(2*nn+1)),(2*nn+1)*(jx-1)+1:(2*nn+1)*jx) = tag_ext(yc(j)+(0:2*nn),xc(j)+(0:2*nn));
                    molim((2*jy-1)*(2*nn+1)+(1:(2*nn+1)),(2*nn+1)*(jx-1)+1:(2*nn+1)*jx) = imm_c_ext(yc(j)+(0:2*nn),xc(j)+(0:2*nn));
                end
            end
        end
    else
        molim=zeros(size(tag));
    end
%         mim(molim)
%         for jy=1:ny-1
%             line([0.5 (2*nn+1)*nx+0.5],(jy*(2*nn+1)*2-0.5)*[1 1],'color','yellow')
%         end
%         for jy=1:ny
%             line([0.5 (2*nn+1)*nx+0.5],((2*jy-1)*(2*nn+1)-0.5)*[1 1],'color','yellow','linewidth',0.3)
%         end
%         for jx=1:nx-1
%             line(((2*nn+1)*jx+0.5)*[1 1],[0.5 ny*2*(2*nn+1)+0.5],'color','yellow','linewidth',0.3)
%         end
end