clear;
close;
floor_ber = 10^-12;
Eb_N0_dB = 0:1:12; %faixa de Eb/N0
Eb_N0_lin = 10 .^ (Eb_N0_dB/10); %faixa de Eb/N0 linearizada
ber_base = zeros(size(Eb_N0_lin)); %prÃ©-alocaÃ§Ã£o do vetor de BER
ber_reed = zeros(size(Eb_N0_lin)); %prÃ©-alocaÃ§Ã£o do vetor de BER
Eb = 1; % energia por bit para a modulaÃ§Ã£o BPSK utilizada
rng(1993);     % Seed random number generator for repeatable results
M = 8;         % Modulation order
bps = log2(M); % Bits per symbol
N = 255;         % RS codeword length
K = 223;         % RS message length
chunk_size = 250*log2(N+1);
Es = 1;
rate = K/N;

num_b = bps*K*log2(N+1)*250; %nÃºmero de bits a serem simulados
bits = complex(2*randi(2, 1, num_b)-3, 0); %bits aleatÃ³rios modulados em BPSK (parte real em 1 e -1)


Eb_8psk_base = Es / bps;
NP_8psk_base = Eb_8psk_base ./ Eb_N0_lin;
NA_8psk_base = sqrt(NP_8psk_base);
Eb_8psk_RS = Es /(bps*rate);
NP_8psk_RS = Eb_8psk_RS ./ Eb_N0_lin;
NA_8psk_RS = sqrt(NP_8psk_RS);

rsEncoder =  comm.RSEncoder('BitInput',true, 'CodewordLength',N, 'MessageLength',K);
rsDecoder = comm.RSDecoder('BitInput',true, 'CodewordLength',N, 'MessageLength',K);
pskmod = comm.PSKModulator('ModulationOrder',M, 'SymbolMapping','Gray', 'PhaseOffset',0, 'BitInput',true);
pskdemod = comm.PSKDemodulator('ModulationOrder',M,'SymbolMapping','Gray','PhaseOffset',0, 'BitOutput',true, 'DecisionMethod','Hard decision');
data = randi([0 1], num_b, 1);

modData = step(pskmod,data);
encodedData = step(rsEncoder,data);
encodedModData = step(pskmod,encodedData);

%ber
for i = 1:length(Eb_N0_lin)
    n = NA_8psk_base(i)*complex(randn(1, num_b/bps), randn(1, num_b/bps))*sqrt(0.5);
    recebido = modData + n'; 
    demod = step(pskdemod,recebido);
    errors = sum(demod ~= data);
    ber_base(i) = errors / num_b; % contagem de erros e cÃ¡lculo do BER
end


%ber reed
for i = 1:length(Eb_N0_lin)
    len_noise = size(encodedModData,1);
    n = NA_8psk_RS(i)*complex(randn(1, len_noise), randn(1, len_noise))*sqrt(0.5); %vetor de ruÃ­do complexo com desvio padrÃ£o igual a uma posiÃ§Ã£o do vetor NA
    recebido = encodedModData+n';
    demod = step(pskdemod,recebido); 
    decoded = step(rsDecoder,demod);
    errors = sum(decoded ~= data);
    ber_reed(i) = errors / num_b; % contagem de erros e cÃ¡lculo do BER
end

ber_reed(ber_reed == 0) = floor_ber;
ber_theoretical = berawgn(Eb_N0_dB,'psk',8,'nondiff');
semilogy(Eb_N0_dB, ber_base, 'x',Eb_N0_dB, ber_reed, 'o', Eb_N0_dB, ber_theoretical, 'r', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('Eb/N0 (dB)');
ylabel('BER');
ylim([floor_ber 10^0]);
legend('Simulado','Reed','TeÃ³rico');



