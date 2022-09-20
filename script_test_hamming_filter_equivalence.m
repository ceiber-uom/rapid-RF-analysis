
raised_cos = @(w) 0.54 - 0.46*cos(2*pi*w);
hamming_window = @(L) raised_cos(linspace(0,1,L));

dual_win = @(L,a) sqrt( raised_cos(linspace(0,1,ceil(L))) .* ...
                        raised_cos(linspace(0-a,1+a,ceil(L))) ); 

u_ = @(f) f/sum(f(:)); 

n = 2e3; 
wave = linspace(1,1,n) .* exp(2i*pi*rand(1,n)); 
wave(1) = 1; 
wave = ifft(wave,'symmetric');

with_ = @(k) conv(wave, u_(k), 'same');

w_h3 = with_(hamming_window(3));
w_dw = with_(dual_win(10,1.2));

butter_gain = @(fc) 1./sqrt(1 + (linspace(0,1,n)/fc).^2); 
w_bf = ifft(fft(wave) .* butter_gain(0.6),'symmetric');


fn = linspace(-0.5,0.5,n);

cla, hold on
plot(fn, abs(fftshift(fft(wave))))
plot(fn, abs(fftshift(fft(w_bf))))
plot(fn, abs(fftshift(fft(w_h3))))
plot(fn, abs(fftshift(fft(w_dw))))

%%


w_tgt = ifft(fft(wave) .* butter_gain(0.6),'symmetric');

target = fft(w_tgt);

with_ = @(k) conv(wave, u_(k), 'same');

gof = @(p) mean(abs(log(fft(with_(dual_win(p(1),p(2))))./target))); 
% p_bwf = fmincon(gof, [5 10], [],[],[],[],[4 -11],[11 11]);

c = turbo(201); 
w = linspace(-1,1,201); 

style = {'linewidth',1.2,'color'};

cla, hold on
plot(fn, abs(fftshift(fft(wave))))
plot(fn, abs(fftshift(target)), style{:}, 'k');

for ii = 1:numel(w)
    y = with_(dual_win(3,w(ii)));
    plot(fn, abs(fftshift(fft(y))), style{:}, [c(ii,:) 0.5]);
end

