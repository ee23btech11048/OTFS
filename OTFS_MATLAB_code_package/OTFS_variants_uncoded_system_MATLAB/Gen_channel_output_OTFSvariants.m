%References
%  [R1]. Y. Hong, T. Thaj, E. Viterbo, ``Delay-Doppler Communications: Principles and Applications'', Academic Press, 2022, ISBN:9780323850285

function r=Gen_channel_output_OTFSvariants(N,M,L_set,gs,s,sigma_2,variant,length_CP)
r=zeros(N*M,1);
noise= sqrt(sigma_2/2)*(randn(size(s)) + 1i*randn(size(s)));
if(strcmp(variant,'ZP')) %Section 4.5.4, Chapter 4, [R1]
    for q=0:N*M-1
        for l=L_set
            if(q>=l)
                r(q+1)=r(q+1)+gs(l+1,q+1)*s(q-l+1);    %see Matlab code 18, Appendix C.4 in [R1]
            end                                        %Remark: slight difference from code 18: ZP is conisdered to be part of the frame so the frame size is NxM (type(c) in Fig 4.13 in [R1])
        end
    end
elseif(strcmp(variant,'CP')) %Section 4.5.3, Chapter 4, [R1]
    for q=0:N*M-1
        n=floor(q/M);
        m=mod(q,M);
        for l=L_set
            r(q+1)=r(q+1)+gs(l+1,m+n*(M+length_CP)+1)*s(n*M+mod(m-l,M)+1);  %see Matlab code 17, Appendix C.4 in [R1] 
        end                                                                 %Remark: CP is NOT conisdered to be part of the frame and is discared, so the transmit frame size is Nx(M+L_cp) (type (b) in Fig 4.13 in [R1])
    end
elseif(strcmp(variant,'RZP')) %Section 4.5.1, Chapter 4, [R1]
    for q=0:N*M-1
        for l=L_set
            if(q>=l)
                r(q+1)=r(q+1)+gs(l+1,q+1)*s(mod(q-l,N*M)+1); %see Matlab code 15, Appendix C.4 in [R1]
            end
        end
    end
elseif(strcmp(variant,'RCP')) %Section 4.5.2, Chapter 4, [R1]
    for q=0:N*M-1
        for l=L_set
            r(q+1)=r(q+1)+gs(l+1,q+1)*s(mod(q-l,N*M)+1); %see Matlab code 16, Appendix C.4 in [R1]
        end        
    end
else
    msg = 'Choose a valid OTFS variant: (RZP / RCP / CP / ZP)';
    error(msg)
end
r=r+noise;
end