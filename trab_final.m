clear;
close all;
floor_ber = 10^-12; %Menor valor do eixo y
Eb_N0_dB = 0:1:12; %faixa de Eb/N0
Eb_N0_lin = 10 .^ (Eb_N0_dB/10); %faixa de Eb/N0 linearizada
ber_base8 = zeros(size(Eb_N0_lin)); %vetor de BER para 8psk
ber8_c1 = zeros(size(Eb_N0_lin)); %codigo convolucional1
ber8_c2 = zeros(size(Eb_N0_lin)); %codigo convolucional2
ber8_c3 = zeros(size(Eb_N0_lin)); %codigo convolucional3

ber_base16 = zeros(size(Eb_N0_lin)); %vetor de BER para 16PSK
ber16_c1 = zeros(size(Eb_N0_lin)); %codigo convolucional1
ber16_c2 = zeros(size(Eb_N0_lin)); %codigo convolucional2
ber16_c3 = zeros(size(Eb_N0_lin)); %codigo convolucional2
Eb = 1; % energia por bit para a modulaÃ§Ã£o BPSK utilizada
rng(1993);     % seed randomica
M = 8;         % Modulation order 1
M2 = 16; %Modulation order 2
bps = log2(M); % 8PSK bps
bps2 = log2(M2); %16PSK bps
Es = 1;
num_b = 1020000/8; %número de bits a serem simulados

%Trellis para os códigos de convolução
trellis1 = poly2trellis(7,[171 133]);
r1 = 1/2;
trellis2 = poly2trellis(9,[753 561]);
r2 = 1/2;
trellis3 = poly2trellis(15,[46321 51271 63667 70535]);
r3 = 1/4;



%Vetores de ruídos para cada caso de simulação
NA_8psk_base = noiseA(Es,bps, Eb_N0_lin,1);
NA_8psk_c1 =  noiseA(Es,bps, Eb_N0_lin,r1);
NA_8psk_c2 =  noiseA(Es,bps, Eb_N0_lin,r2);
NA_8psk_c3 =  noiseA(Es,bps, Eb_N0_lin,r3);

NA_16psk_base =  noiseA(Es,bps2, Eb_N0_lin,1);
NA_16psk_c1 =  noiseA(Es,bps2, Eb_N0_lin,r1);
NA_16psk_c2 = noiseA(Es,bps2, Eb_N0_lin,r2);
NA_16psk_c3 = noiseA(Es,bps2, Eb_N0_lin,r3);

%Moduladores/Demoduladores 8PSK e 16PSK
pskmod = comm.PSKModulator('ModulationOrder',M, 'SymbolMapping','Gray', 'PhaseOffset',0, 'BitInput',true);
pskdemod = comm.PSKDemodulator('ModulationOrder',M,'SymbolMapping','Gray','PhaseOffset',0, 'BitOutput',true, 'DecisionMethod','Hard decision');
pskmod2 = comm.PSKModulator('ModulationOrder',M2, 'SymbolMapping','Gray', 'PhaseOffset',0, 'BitInput',true);
pskdemod2 = comm.PSKDemodulator('ModulationOrder',M2,'SymbolMapping','Gray','PhaseOffset',0, 'BitOutput',true, 'DecisionMethod','Hard decision');

%Dados a serem enviados. Um não codificado e três codificados por convolução
data = randi([0 1], num_b, 1);
convdata1 = convenc(data,trellis1);
convdata2 = convenc(data,trellis2);
convdata3 = convenc(data,trellis3);

%Dados modulados. (3 convoluidos +1 puro)* 2 tipos de modulação
modData = step(pskmod,data);
modData2 = step(pskmod2,data);
modData_conv_1 = step(pskmod,convdata1);
modData_conv_2 = step(pskmod,convdata2);
modData_conv_3 = step(pskmod,convdata3);
modData2_conv_1 = step(pskmod2,convdata1);
modData2_conv_2 = step(pskmod2,convdata2);
modData2_conv_3 = step(pskmod2,convdata3);


%ber 8PSK puro
ber_base8 = calculate_ber(NA_8psk_base,modData, data,Eb_N0_lin,pskdemod);

%ber 8PSK C1
ber8_c1 = calculate_ber_conv(NA_8psk_c1,modData_conv_1, trellis1,data,Eb_N0_lin,pskdemod);

%ber 8PSK C2
ber8_c2 = calculate_ber_conv(NA_8psk_c2,modData_conv_2, trellis2,data,Eb_N0_lin,pskdemod);
disp('8psk C3')

%ber 8PSK C3
ber8_c3 = calculate_ber_conv(NA_8psk_c3,modData_conv_3, trellis3,data,Eb_N0_lin,pskdemod);

disp('16 psk')

%ber 16PSK
for i = 1:length(Eb_N0_lin)
    n = NA_16psk_base(i)*complex(randn(1, num_b/bps2), randn(1, num_b/bps2))*sqrt(0.5);
    recebido = modData2 + n'; 
    demod = step(pskdemod2,recebido);
    errors = sum(demod ~= data);
    ber_base16(i) = errors / num_b; % contagem de erros e cÃ¡lculo do BER
end

%ber 16PSK C1
for i = 1:length(Eb_N0_lin)
    n = NA_16psk_c1(i)*complex(randn(1, (num_b/bps2)/r1), randn(1, (num_b/bps2)/r1))*sqrt(0.5);
    recebido = modData2_conv_2 + n'; 
    demod = step(pskdemod2,recebido);
    decoded = vitdec(demod,trellis1,4,'trunc','hard');
    errors = sum(decoded ~= data);
    ber16_c1(i) = errors / num_b; % contagem de erros e cÃ¡lculo do BER
end

%ber 16PSK C2
for i = 1:length(Eb_N0_lin)
    n = NA_16psk_c2(i)*complex(randn(1, (num_b/bps2)/r2), randn(1, (num_b/bps2)/r2))*sqrt(0.5);
    recebido = modData2_conv_2 + n'; 
    demod = step(pskdemod2,recebido);
    decoded = vitdec(demod,trellis2,4,'trunc','hard');
    errors = sum(decoded ~= data);
    ber16_c2(i) = errors / num_b; % contagem de erros e cÃ¡lculo do BER
end

%{
disp('16psk C3')
%ber 16PSK C3
for i = 1:length(Eb_N0_lin)
    disp(i);
    n = NA_16psk_c3(i)*complex(randn(1, (num_b/bps2)/r3), randn(1, (num_b/bps2)/r3))*sqrt(0.5);
    recebido = modData2_conv_3 + n'; 
    demod = step(pskdemod2,recebido);
    decoded = vitdec(demod,trellis3,4,'trunc','hard');
    errors = sum(decoded ~= data);
    ber16_c3(i) = errors / num_b; % contagem de erros e cÃ¡lculo do BER
end
%}

%Vetor ber teorico para cada modulação (não codificado)
ber_theoretical_8 = berawgn(Eb_N0_dB,'psk',8,'nondiff');
ber_theoretical_16 = berawgn(Eb_N0_dB,'psk',16,'nondiff');
figure
semilogy(Eb_N0_dB, ber_base8,'x', Eb_N0_dB, ber8_c1,'o',Eb_N0_dB, ber8_c2,'*',Eb_N0_dB, ber8_c3,'diamond',Eb_N0_dB, ber_theoretical_8, 'r', 'LineWidth', 2, 'MarkerSize', 10);
legend('Simulado 8psk','Conv K=7','Conv K=9','Conv K=15','Teorico 8psk');
ylim([floor_ber 10^0]);
figure
semilogy(Eb_N0_dB, ber_base16, 'x', Eb_N0_dB, ber16_c1,'o',Eb_N0_dB, ber16_c2,'*',Eb_N0_dB, ber16_c3,'diamond', Eb_N0_dB, ber_theoretical_16, 'g', 'LineWidth', 2, 'MarkerSize', 10);
legend('Simulado 16psk','Conv K=7','Conv K=9','Conv K=15','Teorico 16psk');
xlabel('Eb/N0 (dB)');
ylabel('BER');
ylim([floor_ber 10^0]);



