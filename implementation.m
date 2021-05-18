%{
    Term-Paper ELL319 
    "Pathological speech signal analysis and classification using empirical mode decomposition"
    Group - 10
    Members - Prem Prakash, Siddhant Haritwal, Abhinav Kumar, Ankit kumar Meena, Ashish Bhaskar
    TA - Ayush Tripathi
%}
clc;
close all;
warning off;

universal = [];
% energy_array_imfs = [];
% spread_inst_temporal_energy_density_array = [];
% deviation_inst_temporal_energy_density_array = [];
% spread_of_value_inst_freq_array = [];
% deviation_of_instanteneous_spectral_energy_array = [];

file = uigetfile('*.wav','Grab','Multiselect','on');

for g=1:length(file)
    %samples = [1,16000];
    [yx,fs] = audioread(cell2mat(file(1,g)));
    y1 = resample(yx,16000,fs);
    y = y1(16001:40000);
    
%{    
    %play sound
    % sound(y,fs);

    %plot of input signal in time domain
    t=(0:length(y)-1)/fs;
    figure;
    plot(t,y); 
    title('Input audio signal in time domain'); 
    xlabel('time (seconds)'); 
    ylabel('Amplitude'); 
    
    %plot of input signal in frequency domain
    nfft = 32000; %length for fourier transform
    f = linspace(0,32000,nfft);
    yy = fft(y,nfft);
    figure;
    plot(f,yy);
    title('Input audio signal in freq. domain');
    xlabel('frequncy domain');
    ylabel('Time Samples');

    %plot of EMD
    emd(y);
%}
   
    %Emperical Mode Decomposition (EMD) 
    [imf,residual,info] = emd(y);
    
    % 2.7 Feature extraction
    r_c = size(imf);
    
    for i = 1:r_c(2)  %loop for the number of imf
        e = 0;
        s = 0;
        d = 0;
        sq = 0;
        D_ISED = 0;
        x_ni = zeros(length(imf),1);
        imf_fn = imf(1:length(imf),i);
        %d_temp = diff(imf);
        ht = hilbert(imf_fn);
        theta = [];%initialize
        for k = 1:length(imf_fn)
            theta(k,1) = atan(ht(k)/imf_fn(k));
        end
        ai = []; %initialize
        for b = 1:length(imf_fn)
            ai(b,1) = sqrt(imf_fn(b)*imf_fn(b) + ht(b)*ht(b));
        end
        omega = diff(theta);
        s_omega=trapz(omega); %trapz gives the area
        for n = 1 : (r_c(1)-1) %loop for the nubmer of points in each imf
            e = e + imf_fn(n)*imf_fn(n);
            new_entry = ai(n)*ai(n) - conj(ai(n))*conj(ai(n));
            s = s + new_entry;
            new_entry1 = n*ai(n);
            d = d + new_entry1;
            new_entry2 = omega(n)*omega(n) - conj(omega(n))*conj(omega(n));
            sq = sq + new_entry2;
            new_entry3 = ai(n)*omega(n);
            D_ISED = D_ISED + new_entry3;
            x_ni(n) = (ai(n)*exp(1j*s_omega));
        end

        s = s/r_c(1);
%         energy_array_imfs(g,i) = e;
%         spread_inst_temporal_energy_density_array(g,i) = s;
        d = d/r_c(1);
%         deviation_inst_temporal_energy_density_array(g,i) = d;
        sq = sq/r_c(1);
%         spread_of_value_inst_freq_array(g,i) = sq;
        D_ISED = D_ISED/r_c(1);
%         deviation_of_instanteneous_spectral_energy_array(g,i) = D_ISED;
        
        universal(g,i) = e;
        universal(g,i+10) = s;
        universal(g,i+20) = d;
        universal(g,i+30 ) = sq;
        universal(g,i+40) = D_ISED;
        universal(g,i+50) = real(sum(x_ni));
        
    end

    % sum() do the column wise sum of the matrix
    % sum(A,dim) returns the sum along dimension dim

    % %Use the 'Display' name-value pair to output a table showing the number of sifting iterations, the relative tolerance, and the sifting stop criterion for each IMF
    %imf = emd(y,'Display',1);
    % %Hilbert spectrum of a signal sampled at sapmling frequency 16000
    %hht(imf,16000);
    
end

writematrix(universal, 'Healthy_B1.csv');
