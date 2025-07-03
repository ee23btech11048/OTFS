
function [nu_ml_tilda]=Gen_delay_time_channel_vectors_OTFSvariants(N,M,l_max,gs,L_cp,variant)
nu_ml_tilda=zeros(N,M,l_max+1);
for n=1:N
    for m=1:M
        for ell=1:(l_max+1)
            if(strcmp(variant,'RZP'))
                if(not(m<ell && n==1))
                    nu_ml_tilda(n,m,ell)=gs(ell,m+(n-1)*M);
                end
            elseif(strcmp(variant,'RCP'))
                nu_ml_tilda(n,m,ell)=gs(ell,m+(n-1)*M);
            elseif(strcmp(variant,'CP'))
                nu_ml_tilda(n,m,ell)=gs(ell,m+(n-1)*(M+L_cp));
            elseif(strcmp(variant,'ZP'))
                if(m>=ell)
                    nu_ml_tilda(n,m,ell)=gs(ell,m+(n-1)*M);
                end
            end
        end
    end
end
end