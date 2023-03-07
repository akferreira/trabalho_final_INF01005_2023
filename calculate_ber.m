function [ber] = calculate_ber(NA,modData, userData,Eb_N0_lin,pskdemod)
ber = zeros(length(Eb_N0_lin),1);

for i = 1:length(Eb_N0_lin)
    n = NA(i)*complex(randn(1, size(modData,1)), randn(1, size(modData,1)))*sqrt(0.5);
    recebido = modData + n'; 
    demod = step(pskdemod,recebido);
    errors = sum(demod ~= userData);
    ber(i) = errors / size(userData,1); % contagem de erros e cálculo do BER
end