function [r, err] = cond_xy(x,y,condn)
if ~isempty(condn)
    condn=[condn ones(size(x))];
    
    by = regress(y,condn);
    bx = regress(x,condn);
    yfit=0; xfit=0;
    for icnt =1 : length(bx)
        yfit = by(icnt)*condn(:,icnt)+yfit;
        xfit = bx(icnt)*condn(:,icnt)+xfit;
    end
    
    ry= y-yfit;
    rx= x-xfit;
    
else
    ry= y;
    rx= x;
end

[r, err]=calc_err(rx, ry);