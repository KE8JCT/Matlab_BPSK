%The following is an attempt to modulate and demodulate binary phase-shift
%keying
string1 = 'usma';
b = flip(de2bi(double(string1)),2);

b = reshape(b.',1, []);

B = (b + .5)*10;
c = ones(1, (length(B)*10 + 1));
C = c

for d = 0:(length(B)-1)
    for e = 1:10
        C(e+(d*10)) = B((d+1));
    end
end



tt = 0:1/length(b):10
%F = C(:)
%xx = sin(2*pi*C.*tt);
%yy = sin(2*pi*F*tt);
%{
%To plot two B-FSK Signals
tiledlayout(2, 1)
nexttile
plot(tt, xx)
title('Binary FSK signal')
xlabel('TIME (sec)')
nexttile
plot(tt, yy)
title('Binary FSK signal - Blocky')
xlabel('TIME (sec)')
%}

%Binary Phased-Shift Keying

for d = 0:(length(b)-1)
    for e = 1:10
        c(e+(d*10)) = b((d+1));
    end
end

for i = 1:length(c)
    if c(i) == 0
        c(i) = -1;
    end
    
end

%add awgn



zz = c.*sin(2*pi*1.4*tt);
zn = awgn(zz, 20) %SNR of 20
xx = sin(2*pi*.789*tt+pi)


%{
plot(tt, zz)
hold on
plot(tt, c)
%yline(1:(10/28):10)
%plot(tt, abs(zz))
%stem(tt, ipt)
title('Binary Phased-Shift Keying')
xlabel('TIME (sec)')
xlim([0 10])
ylim([-1.2 1.2])
hold off
%}

%Demodulation

windowSize = 5; 
coef_a= (1/windowSize)*ones(1,windowSize);
coef_b = 1;
zf = filter(coef_a,coef_b,zn)
z2 = abs(zf);
[pks, locs] = findpeaks(z2);

T = length(pks);

ts = 0:(1/T):10
test_signal = sin(pi*T/10*ts) %clock/carrier frequency

product_detector = test_signal.*zn %product detector

%Input low-pass filter here

%threshold it:
bitstream = zeros(1,length(product_detector));
for i=1:length(product_detector)
   
    if product_detector(i) > 0
        bitstream(i) = 1 
    end
    
    if product_detector(i) <= 0
         bitstream(i) = 0
    end
    %now to reverse the "upsampling" and get an actual binary stream
end

binary_data = bitstream(5:10:end)

%Now get the values into separate vectors and convert from binary to ascii
%to text

%ascii is always 7 bits long

final_binary = [];
temp_vector = [];
r = rem(length(binary_data), 7)

for i=1:length(binary_data)
  temp_vector = [temp_vector, binary_data(i)]
  if (mod(i,7) == 0)
      final_binary = [final_binary; temp_vector]
      temp_vector = [];
  end
end

%decimal_string = num2str(final_binary);
%decimal_string = bin2dec(decimal_string);
decimal_string = num2str(final_binary);
decimal_string = bin2dec(decimal_string);
output_vector = [];

for i = 1:length(decimal_string)
    output_vector = [output_vector; char(decimal_string(i))]
end

%{
tiledlayout(2,1)
nexttile
plot(test_signal)
hold on
plot(zz)
plot(-test_signal)
legend('clock signal','received signal', 'negative clock')
hold off
nexttile
plot(zz)
%}

%{
xnor_flag = [];
loc_flag = [];
for i=1:2:length(zz) %for the whole length of the signal
   if zz(i) == test_signal(i) %compare received signal to the clock signal
      xnor_flag = [xnor_flag 0];
      loc_flag = [loc_flag 0];
   end
   if zz(i) ~= test_signal(i)
       xnor_flag = [xnor_flag zz(i)];
       loc_flag = [loc_flag i];
   end
    
end
%}


tiledlayout(2,1)
nexttile
%{
plot(zz)
hold on
plot(test_signal)
hold off
nexttile
%}

plot(bitstream)
hold on
plot(c)
plot(zz)
xlim([0 280])
ylim([-1.2 1.2])
legend('bitstream','original data', 'transmitted BPSK Signal')
title('Transmission of BPSK')
hold off
nexttile

plot(zn, 'r')
hold on
plot(c, 'g')
plot(zf, 'b')
legend('Received Signal','Binary Bit Stream','Filtered RX Signal')
xlim([0 280])
ylim([-1.2 1.2])
title('Reception of BPSK - SNR 20')
hold off
