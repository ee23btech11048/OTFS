function [hEnc,hDec,hDec_coded_soft,hDec_coded_hard]=LDPC_system_objects(LDPC_rate,LDPC_codeword_length)

H_parity=select_LDPC_parity_matrix(LDPC_rate,LDPC_codeword_length);

hEnc = comm.LDPCEncoder('ParityCheckMatrix',sparse(logical(H_parity)));

hDec = comm.LDPCDecoder('ParityCheckMatrix',sparse(logical(H_parity)), 'IterationTerminationCondition','Parity check satisfied');

hDec_coded_soft = comm.LDPCDecoder('ParityCheckMatrix',sparse(logical(H_parity)), 'IterationTerminationCondition','Parity check satisfied', 'DecisionMethod', 'Soft decision', 'OutputValue', 'Whole codeword');
hDec_coded_soft.FinalParityChecksOutputPort = 1;

hDec_coded_hard = comm.LDPCDecoder('ParityCheckMatrix',sparse(logical(H_parity)), 'IterationTerminationCondition','Parity check satisfied', 'DecisionMethod', 'Hard decision', 'OutputValue', 'Whole codeword');
hDec_coded_hard.FinalParityChecksOutputPort = 1;

end