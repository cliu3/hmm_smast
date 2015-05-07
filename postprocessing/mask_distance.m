function DistT=mask_distance(mask)
% compute distance map from matrix center with obstacles
% Chang Liu
import java.util.LinkedList
% mask=zeros(17);
% mask(6,1:5)=1;
% mask(2:5,5)=1;

[rows,cols]=size(mask);
ic=floor((rows+1)/2);
jc=floor((cols+1)/2);

DistT=nan(rows,cols);
DistT(ic,jc)=0;

q=LinkedList();
q.add([ic,jc]);

while any(isnan(DistT(mask(:)==0)))
    try
        cell=q.remove();
    catch
        break
    end
    i=cell(1);
    j=cell(2);
    for ii=-1:1
        for jj=-1:1
            if (ii==0 && jj==0)
                continue
            end
            
            if (i+ii > 0 && j+jj > 0 && i+ii <= rows && j+jj <= cols && mask(i+ii,j+jj)==0)
                if ( abs(ii)+abs(jj)>1 && mask(i,j+jj)==1 && mask(i+ii,j)==1 )
                    continue
                end
                if(ii==0 || jj==0)
                    dist_nbn = 1.0;
                else
                    dist_nbn = 1.414;
                end
                if ~(DistT(i,j)+dist_nbn>=DistT(i+ii,j+jj))
                    DistT(i+ii,j+jj)=DistT(i,j)+dist_nbn;
                    q.add([i+ii,j+jj]);
                end
            end
        end
    end
    
    
    
    
end



end

