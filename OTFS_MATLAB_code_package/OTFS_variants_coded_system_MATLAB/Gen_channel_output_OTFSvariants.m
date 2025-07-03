function r=Gen_channel_output_OTFSvariants(N,M,L_set,gs,s,sigma_2,length_CP,variant)
r=zeros(N*M,1);
noise= sqrt(sigma_2/2)*(randn(size(s)) + 1i*randn(size(s)));
if(strcmp(variant,'ZP'))
    for q=0:N*M-1
        for l=L_set
            if(q>=l)
                r(q+1)=r(q+1)+gs(l+1,q+1)*s(q-l+1);    % see Section 4.5.4, Chapter 4, [R1]
            end
        end
    end
elseif(strcmp(variant,'CP'))
    for q=0:N*M-1
        n=floor(q/M);
        m=mod(q,M);
        for l=L_set
            r(q+1)=r(q+1)+gs(l+1,m+n*(M+length_CP)+1)*s(n*M+mod(m-l,M)+1);
        end
    end
elseif(strcmp(variant,'RZP'))
    for q=0:N*M-1
        for l=L_set
            if(q>=l)
                r(q+1)=r(q+1)+gs(l+1,q+1)*s(mod(q-l,N*M)+1);
            end
        end
    end
elseif(strcmp(variant,'RCP'))
    for q=0:N*M-1
        for l=L_set
            r(q+1)=r(q+1)+gs(l+1,q+1)*s(mod(q-l,N*M)+1);
        end        
    end
else
    msg = 'Choose a valid OTFS variant: (RZP / RCP / CP / ZP)';
    error(msg)
end
r=r+noise;
end