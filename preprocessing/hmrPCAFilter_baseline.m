%% baseline PCA filter based on hmrPCAFilter from Homer2
% do PCA on baseline data, project real data onto PCA coordinate system, 
% remove nSV components. 
%
%removed wavelength data for simplicity
%removed var expl routine, take only numerical SV/components
%integrated baseline code snipplets from "bottom"
%
%this is only applicable to this specific case, additional
%functions removed for simplicity - handle with care!
%
% KK August 2024

function [yc, svs, nSV] = hmrPCAFilter_baseline(yData, yBaseline, SD, nSV )

if ~exist('nSV')
    disp('USAGE: [yc,svs] = hmrPCAFilter( y, SD, nSV )');
    return
end

ml = SD.MeasList;
nMeas = size(ml,1);


if ~isfield(SD,'MeasListAct')
    SD.MeasListAct = ones(nMeas,1);
end

% do the PCA
%ndim = ndims(y);

    % PCA on Concentration
    lstAct = find(SD.MeasListAct==1);
    lstAct = lstAct( find(ml(lstAct,4)==1) );
    ybas = yBaseline(:,:,lstAct);
    yd = yData(:,:,lstAct);
    yc = yData;
    for iConc = 1:2
        yb = squeeze(ybas(:,iConc,:));
        yb = detrend(yb);
        y = squeeze(yd(:,iConc,:));
        y = detrend(y);
        
        %baseline PCA
        c = yb.' * yb;
        [v,s,foo] = svd(c);
    %    u = y*v*inv(s);
    
    
        svs(:,iConc) = diag(s)/sum(diag(s));
        
        if nSV(iConc)<1 % find number of SV to get variance up to nSV
            svsc = svs(:,iConc);
            for idx = 2:size(svs,1)
                svsc(idx) = svsc(idx-1) + svs(idx,iConc);
            end
            ev = diag(svsc<nSV(iConc));
            nSV(iConc) = find(diag(ev)==0,1)-1;
        else
            ev = zeros(size(svs,1),1);
            ev(1:nSV(iConc)) = 1;
            ev = diag(ev);
        end
        lst = 1:nSV(iConc);
        u = y*v*inv(s);
        yc(:,iConc,lstAct) = y - u(:,lst)*s(lst,lst)*v(:,lst)';
    end
    
    yc(:,3,:) = yc(:,1,:) + yc(:,2,:);
    


if isempty(nSV)
    nSV = 0;
end

