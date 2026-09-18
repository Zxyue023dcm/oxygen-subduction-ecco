function [OS,OSw,OSeddy,OSlavoX,OSlavoY] = oxygen_flux(oceanlon,oceanlat,Hml,Hst,evel1,evel2,nvel1,nvel2,wvelb2,estar1,estar2,nstar1,nstar2,prho,prho10,oxy,oxy10)
%% Calculate monthly oxygen-subduction flux
% Inputs are the physical and dissolved-oxygen fields used in the three
% original demo scripts. Outputs are in mol O2 m^-2 month^-1.

Slavo_x=nan(size(Hml));
Slavo_y=nan(size(Hml));
Seddy=nan(size(Hml));

[lonDistance,latDistance]=distcoordinate(oceanlon(:,1),oceanlat(1,:)');
deltaLon=diff(lonDistance,1,1);

%% Lateral induction: u1*grad(Hml) + u2*grad(Hst)
lonBackward=[lonDistance(1,:)-deltaLon(1,:);lonDistance];
latBackward=[latDistance(1,:);latDistance];
lonBackward=[nan(size(lonBackward,1),1) lonBackward];
latBackward=[nan(size(latBackward,1),1) latBackward];

gHmlx1=nan(size(Hml));gHmly1=nan(size(Hml));
gHstx1=nan(size(Hml));gHsty1=nan(size(Hml));

for i=1:size(Hml,3)
    disp(i)
    hml=[Hml(end,:,i);Hml(:,:,i)];
    hml=[nan(size(hml,1),1) hml];
    [~,gHmlx,gHmly]=divdffs(hml,lonBackward,hml,latBackward,'backward');

    hst=[Hst(end,:,i);Hst(:,:,i)];
    hst=[nan(size(hst,1),1) hst];
    [~,gHstx,gHsty]=divdffs(hst,lonBackward,hst,latBackward,'backward');

    % Average paired longitude bands, then fill their shared boundaries.
    for a=1:size(Hml,1)/2
        gHmlx1(2*a-1:2*a,:,i)=repmat(mean(gHmlx(2*a-1:2*a,:), ...
            1,'omitnan'),[2 1]);
        gHstx1(2*a-1:2*a,:,i)=repmat(mean(gHstx(2*a-1:2*a,:), ...
            1,'omitnan'),[2 1]);
    end
    for a=1:size(Hml,1)/2-1
        gHmlx1(2*a:2*a+1,:,i)=repmat(mean(gHmlx1(2*a:2*a+1,:,i), ...
            1,'omitnan'),[2 1]);
        gHstx1(2*a:2*a+1,:,i)=repmat(mean(gHstx1(2*a:2*a+1,:,i), ...
            1,'omitnan'),[2 1]);
    end
    for a=1:size(Hml,1)/2
        gHmlx1(2*a-1:2*a,:,i)=repmat(mean(gHmlx1(2*a-1:2*a,:,i), ...
            1,'omitnan'),[2 1]);
        gHstx1(2*a-1:2*a,:,i)=repmat(mean(gHstx1(2*a-1:2*a,:,i), ...
            1,'omitnan'),[2 1]);
    end

    % Apply the same averaging to paired latitude bands.
    for b=1:size(Hml,2)/2
        gHmly1(:,2*b-1:2*b,i)=repmat(mean(gHmly(:,2*b-1:2*b), ...
            2,'omitnan'),[1 2]);
        gHsty1(:,2*b-1:2*b,i)=repmat(mean(gHsty(:,2*b-1:2*b), ...
            2,'omitnan'),[1 2]);
    end
    for b=1:size(Hml,2)/2-1
        gHmly1(:,2*b:2*b+1,i)=repmat(mean(gHmly1(:,2*b:2*b+1,i), ...
            2,'omitnan'),[1 2]);
        gHsty1(:,2*b:2*b+1,i)=repmat(mean(gHsty1(:,2*b:2*b+1,i), ...
            2,'omitnan'),[1 2]);
    end
    for b=1:size(Hml,2)/2
        gHmly1(:,2*b-1:2*b,i)=repmat(mean(gHmly1(:,2*b-1:2*b,i), ...
            2,'omitnan'),[1 2]);
        gHsty1(:,2*b-1:2*b,i)=repmat(mean(gHsty1(:,2*b-1:2*b,i), ...
            2,'omitnan'),[1 2]);
    end

    Slavo_x(:,:,i)=gHmlx1(:,:,i).*evel1(:,:,i) ...
        +gHstx1(:,:,i).*evel2(:,:,i);
    Slavo_y(:,:,i)=gHmly1(:,:,i).*nvel1(:,:,i) ...
        +gHsty1(:,:,i).*nvel2(:,:,i);
end

%% Convert lateral volume flux to oxygen flux
conversionFactor=10^(-6)*3600*24*365/12;

OSlavoX=Slavo_x;
positive=Slavo_x>0;
negative=Slavo_x<0;
OSlavoX(positive)=Slavo_x(positive).*oxy10(positive).*prho10(positive)*conversionFactor;
OSlavoX(negative)=Slavo_x(negative).*oxy(negative).*prho(negative)*conversionFactor;

OSlavoY=Slavo_y;
positive=Slavo_y>0;
negative=Slavo_y<0;
OSlavoY(positive)=Slavo_y(positive).*oxy10(positive).*prho10(positive)*conversionFactor;
OSlavoY(negative)=Slavo_y(negative).*oxy(negative).*prho(negative)*conversionFactor;

OSlavoX=-OSlavoX;
OSlavoY=-OSlavoY;

%% Vertical velocity contribution
OSw=wvelb2;
positive=wvelb2>0;
negative=wvelb2<0;
OSw(positive)=wvelb2(positive).*oxy10(positive).*prho10(positive)*conversionFactor;
OSw(negative)=wvelb2(negative).*oxy(negative).*prho(negative)*conversionFactor;
OSw=-OSw;

%% Eddy-induced contribution: div(ustar1*Hml + ustar2*Hst)
lonGradient=[lonDistance(1,:)-deltaLon(1,:);lonDistance; ...
    lonDistance(end,:)+deltaLon(1,:)];
latGradient=[latDistance(1,:);latDistance;latDistance(end,:)];

for i=1:size(Hml,3)
    disp(i)
    estarh=estar1(:,:,i).*Hml(:,:,i)+estar2(:,:,i).*Hst(:,:,i);
    nstarh=nstar1(:,:,i).*Hml(:,:,i)+nstar2(:,:,i).*Hst(:,:,i);
    es=[estarh(end,:);estarh;estarh(1,:)];
    ns=[nstarh(end,:);nstarh;nstarh(1,:)];
    Seddy(:,:,i)=divdffs(es,lonGradient,ns,latGradient,'gradient');
end

OSeddy=Seddy;
positive=Seddy>0;
negative=Seddy<0;
OSeddy(positive)=Seddy(positive).*oxy10(positive).*prho10(positive)*conversionFactor;
OSeddy(negative)=Seddy(negative).*oxy(negative).*prho(negative)*conversionFactor;
OSeddy=-OSeddy;

%% Total oxygen-subduction flux
OS=OSlavoX+OSlavoY+OSw+OSeddy;
end
