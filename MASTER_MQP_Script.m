%% Initialize variables
tic
%User Inputs
a = input('How many basilar membrane locations? (minimum 5) ');
b = input('How many inner hair cells per basilar membrane location? ');
c = input('How many neurons in each population? ');
%Duration?
%Fs?
%SignalIn?
NumFilters = a;
NumCells = b; 
NumNeurons = c; 

% NumFilters = 5;
% NumCells =3; 
% NumNeurons = 5; 

%% Input Signal

%%%%READING PRERECORDED AUDIO SIGNAL%%%%%%%
% Input Signal / Define Variables
[SignalIn Fs] = audioread("AppleOrange.wav"); %reads audio signal stores
% datapoints and sample frequenct
SignalIn = SignalIn(:,1)';% takes first column to isolate datapoints 
SignalIn = 10.*SignalIn(floor(0.5*end):end);
dur = length(SignalIn)/Fs; %Duration of time signal in seconds
t = [0:(1/Fs):(dur-(1/Fs))];
n = length(t); % Length of Simulation (samples)

%%%%%%READING ABRITRARY WAVEFORM%%%%%%%%%%

% dur = 0.91; % sets duration of arbitrary signal
% Fs = 44100; %sets sampling frequency of arbitrary signal
% t = [0:(1/Fs):(dur-(1/Fs))]; %Time vector based on dur and Fs
% SignalIn = sin(2*pi*2.*t); % signal vector
% n = length(t); % Length of Simulation (samples)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Placeholder vectors Spike Train Based on Time
D = cell(NumCells,NumFilters, NumNeurons);  % Spike train recording 
TT = cell(NumCells,NumFilters, NumNeurons); % Time Recording
XX = cell(NumCells,NumFilters, NumNeurons);% Stimulus Record
for i = 1:NumFilters
    for k = 1: NumCells
        for j = 1:NumNeurons
            D{k,i,j} = zeros(1,length(t));
            TT{k,i,j} = zeros(1,length(t)); % Time Recording
            XX{k,i,j} = zeros(1,length(t)); % Stimulus Record
        end
    end
end
%% Basilar Membrane
% Create Gamma Filter Bank and Filter Parameters:
m = hz2mel([200,12000]);% Convert [200, 20000] to mel scale.

melVect = linspace(m(1),m(2),NumFilters);% Generate a row vector of 32 values linearly spaced on the mel scale.

hzVect = mel2hz(melVect);% Get equivalent frequencies in hertz.

%Create Filter 
F = cell(1,NumFilters); %placeholder cell
fc = 1; %placeholder value
bw = 1;%placeholder value
for i = 1: NumFilters
    fc = hzVect(i); % Center Frequency of GT Filter (Hz)
    bw = 0.25*hzVect(i); % Bandwidth of GT Filter (Hz)
    F{i} = makefiltgt(fc,bw,Fs); % Create GT Filter
end 

%Apply Filter
FilterValues = cell(1,NumFilters); %FilterValues is a cell that holds the output signals for each filter created
for i = 1: NumFilters
    y_realtime = zeros(n,1); % Placeholder for Output Signal
    for j = 1:length(SignalIn)
        [out F{i}] = applyfilt(SignalIn(j),F{i});
        y_realtime(j) = out;
    end
   FilterValues{i} = y_realtime;
end


%Cascade 1

% Create Duplicate Filters
F2 = cell(1,NumFilters); %placeholder cell
fc = 1; %placeholder value
bw = 1;%placeholder value
for i = 1: NumFilters
    fc = hzVect(i); % Center Frequency of GT Filter (Hz)
    bw = 0.25*hzVect(i); % Bandwidth of GT Filter (Hz)
    F2{i} = makefiltgt(fc,bw,Fs); % Create GT Filter
end 

FilterValues2 = cell(1,NumFilters); %FilterValues is a cell that holds the output signals for each filter created
for i = 1: NumFilters
    y_realtime = zeros(n,1); % Placeholder for Output Signal
    for j = 1:length(SignalIn)
        [out F2{i}] = applyfilt(FilterValues{i}(j),F2{i});
        y_realtime(j) = out;
    end
   FilterValues2{i} = y_realtime;
end

%Compare the filtered signals to the original
filterSum = 0;
for i = 1:NumFilters
    filterSum = filterSum + FilterValues2{1,i};%summing all filters
end


%% Generate ihc
% Initialize an array to store IHC models

%Constant values for IHC characterization
ALPHA = zeros(1,NumFilters);
PHI0 = zeros(1, NumFilters);
J_BIAS = zeros(1,NumFilters);

IHC_populations = cell(NumCells,NumFilters);%Cell to store IHC populations 
%Each column holds all the IHCs in one population 

% Create IHC models for each population
for i = 1:NumFilters
        for j = 1: NumCells
        ALPHA = rand(1,NumFilters); % between 0 and 3
        PHI0 = ones(1, NumFilters);
        J_BIAS = 0.0001*rand(1,NumFilters);
        ALPHA_i = ALPHA(i);
        PHI0_i = PHI0(i);
        J_BIAS_i = J_BIAS(i);
        IHC_populations{j,i} = makeIHC(ALPHA_i, PHI0_i, J_BIAS_i, Fs);
        end 
end

%For loop to create 2D cell of empty vectors length of time (to hold
%Voltage values over the time)
V = cell(1,NumFilters);% cell array for voltage
for i = 1:NumFilters
    for k = 1: NumCells
        V{k,i} = zeros(1,length(t));
    end
end

%Update the voltage in the IHCs overtime and store in V cell
 for i = 1:NumFilters
     for k = 1:NumCells
        for j = 1:length(t)
            IHC = IHC_populations{k,i};
            IHC = updateIHC(FilterValues2{i}(j), IHC);
            IHC_populations{k,i} = IHC;
            % Record membrane voltage (V)
            V{k,i}(j) = IHC.V;
        end
     end
 end

%% Auditory Nerve Fibers
%%%%% Part 1/4: Generating a Population of Neurons %%%%%
% Initializing parameters for Neurons
 AlphaN = zeros(1,NumNeurons);
 PhiN = zeros(1, NumNeurons);
 J_Bias_N = zeros(1, NumNeurons);

% Generate Neurons - Neurogenesis?!? (makelifn)
N = cell(NumCells, NumFilters, NumNeurons); %3D cell to hold neurons (each channel in z direction is population)
for k = 1:NumFilters
    for i = 1:NumCells
        [AlphaN, PHI0_BEST, J_Bias_N] = Script_determineparams(V{i,k}',Fs);
        for j = 1:NumNeurons
            N{i, k, j} = makelifn(AlphaN(j),PHI0_BEST(j),J_Bias_N(j));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2/4: Encoding with a Population of Neurons using updatelifn

% Collect Samples
for i = 1:NumFilters
     for k = 1: NumCells
            for j = 1:NumNeurons
                for l = 1:length(t)
        % Record Time Stamp
        TT{k,i,j}(l) = t(l); 
        % Record Physical Value
        x_t = V{k,i}(l);
        XX{k,i,j}(l) = x_t;
        % Simulate/Update LIFN (updatelifn) and record voltage spikes
        N{k,i,j} = updatelifn(x_t,N{k,i,j},Fs);
        D{k,i,j}(l) = N{k,i,j}.V==1;
                end 
            end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part 3/4: Probe & Determine Decoders For a Population of Neurons
% Parameters for Neural "Probe"
res = 200; % resolution of values implemented in "probing"
minval = -1; % minimum considered value
maxval = 1; % maximum considered value

% Characterize Activation Functions of Population (characterizelifn)
%Store Activation Function in cell A (one for each population of neurons)
%Number of columns of A = NumNeurons
A = cell(NumCells,NumFilters);
for i = 1:NumFilters
     for k = 1: NumCells
         for j = 1:NumNeurons
            [A{k,i}(:,j) X] = characterizelifn(N{k,i,j},minval,maxval,res,Fs);
         end
     end
end

% for i = 1:NumFilters
%     for k = 1:NumCells
%         A{k,i} = 1+A{k,i};
%     end
% end

%Determine Decoders using Activation Functions (one for each population of
%neurons)
PHI = cell(NumCells,NumFilters);
for i = 1:NumFilters
    for k = 1:NumCells
        PHI{k,i} = determinedecoders(A{k,i},X);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 4/4: Decode Value Encoded By the Neural Population

% Obtain Temporal Activation Functions (decodespiketrain)
% Compile temporal activation functions into a matrix At
%   appropriate as input to 'decoderspikerate'
x_hat = cell(NumCells,NumFilters);
for i = 1:NumFilters
    for k = 1: NumCells
            x_hat{k,i} = zeros(length(D{k,i}),NumNeurons); % hard coded, need to change
    end
end

for i = 1:NumFilters
     for k = 1: NumCells
         for j = 1:NumNeurons
             temp = decodespiketrain(N{k,i,j},D{k,i,j},Fs);
            % x_hat{k,i}(:,j) = temp(1:end/2);
            x_hat{k,i}(:,j) = temp(1:length(D{k,i}));
         end
     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decode Recorded Spike Rates (decodespikerate)
x_t_hat = cell(NumCells,NumFilters);

for i = 1: NumFilters
    for k = 1:NumCells
        for j = 1:NumNeurons
        x_t_hat{k,i} = decodespikerate(x_hat{k,i},PHI{k,i}); 
        end
    end
end 

%% 

DecodedSignal = zeros(size(SignalIn,2),1);
for i = 1:NumFilters
    for j = 1:NumCells
        DecodedSignal(:,1) = DecodedSignal(:,1) + x_t_hat{j,i}(:,1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
cutoff_frequency = 50;  
sampling_frequency = Fs;
filtered_signal = DecodedSignal;

for i = 1:10
    filtered_signal = lowpass(filtered_signal,cutoff_frequency,sampling_frequency);

end

    IHCSum = 0;
    for i = 1:NumFilters
        for j = NumCells
        IHCSum = IHCSum + V{j,i};%summing all filters
        end 
    end
    
    %%% Verification 
    CorrelationBasilar = corr(SignalIn',filterSum);
    CorrelationIHC = corr(SignalIn',IHCSum');
    CorrelationDecoded = corr(SignalIn',DecodedSignal);
    CorrelationLowPass = corr(SignalIn',filtered_signal);
    
    Correlation1 = corr(IHCSum',DecodedSignal);
    Correlation2 = corr(IHCSum',filtered_signal);
    
    fprintf('Correlation Between Input Signal and Sum of Basilar Membrane Outputs: %.5f\n', CorrelationBasilar);
    fprintf('Correlation Between Input Signal and Sum of IHC Outputs: %.5f\n', CorrelationIHC);
    fprintf('Correlation Between Input Signal and Decoded Signal: %.5f\n', CorrelationDecoded);
    fprintf('Correlation Between Input Signal and LP Filtered Decoded Signal: %.5f\n', CorrelationLowPass);
    
    %% GRAPHS%%%%%%%%
    %FIGURE 1: Input soundwave
    figure %plot input frequency
    plot(t,SignalIn)
    xlabel('Time (seconds)','FontSize',18)
    ylabel('Amplitude','FontSize',18)

    %FIGURE 2: Plot original signal and sum of all filtered signals 
    figure
    hold on
    %title('Filtered Signal Sum vs Original') % added this for the future for the report
    plot(t,SignalIn,'b')
    plot(t,filterSum,'r')
    xlabel('Time (seconds)','FontSize',18)
    ylabel('Voltage','FontSize',18)

    %FIGURE 3: Plot each individual filter output 
    Function = 'Function ';
    legendLabels = cell(NumFilters, 1);
    Offset = 1;
    figure
    %title('Signal Output from all Filters')
    for i = 1:NumFilters
        %legendLabels{i} = ['Function ' num2str(i)];
        Num = num2str(i);
        hold on
        plot(t,FilterValues2{1,i} + Offset*i)
    end
    xlabel('Time (seconds)','FontSize',18)
    ylabel('Voltage','FontSize',18)
    %legend(legendLabels, 'Location', 'best');

    %FIGURE 4:Plot of all IHCs Updated Voltages 
    figure
    hold on
    for i = 1:NumFilters
         for k = 1: NumCells
             plot(t,V{k,i})
         end
    end
    %title('Individual IHC Voltages')
    xlabel('Time (seconds)','FontSize',18)
    ylabel('Voltage','FontSize',18)

    %FIGURE 5
    NumRows = 5;
    Function = 'Function ';
    legendLabels2 = cell(NumCells, 1);
    Offset = 1;
    figure
    %title('Signal Output from First 5 Locations')
    for i = 1:2
        subplot(NumRows,1,i)
        title(['Basilar Mebrane Location' num2str(i)]);
        for k = 1: NumCells
        %legendLabels2{k} = ['Function ' num2str(k)];
        hold on
        plot(t,V{k,i} + Offset*k)
        xlabel('Time (seconds)','FontWeight','bold')
        ylabel('Voltage','FontWeight','bold')
        hold on
        end
    end

    %legend(legendLabels2, 'Location', 'best');

    %FIGURE 6: IHC ADDition

    figure
    hold on
    plot(t,IHCSum)
    plot(t,SignalIn,'r')
    %title('Sum of IHC Voltages')
    xlabel('Time (seconds)','FontSize',18)
    ylabel('Voltage','FontSize',18)

    % FIGURE 7: Visualize Spikes (of first neuron)
    figure
    hold on
    plot(TT{1}(1:0.1*end),D{1,1,1}(1:0.1*end),'linewidth',2)
    plot(t(1:0.1*end),SignalIn(1:0.1*end))
    %title('Visualize Spikes','FontSize',24)
    xlabel('Time (seconds)','fontsize',18)
    ylabel('Voltage','fontsize',18)

    %Figure 8
    counter = 0;
    NumRows = 5;
    NumCols = 2;
    Function = 'Function ';
    legendLabels3 = cell(NumNeurons, 1);
    figure
    title('Spike Train of first 10 populations')
    for i = 1:5
        for k = 1: 2
        %fprintf('Index %d: i=%d, j=%d\n', counter, i, j);
        counter = counter + 1;

        subplot(NumRows,NumCols,counter)
        title(['Population of Filter ' num2str(i) ' IHC' num2str(k)]);
            for j = 1: NumNeurons
        %legendLabels3{j} = ['Function ' num2str(k)];
        hold on
        plot(TT{i},D{k,i,j})
        xlabel('Time','FontWeight','bold')
        ylabel('Voltage','FontWeight','bold') 
            end 
        end
    end

    %FIGURE 9: Plot Activation functions for given population
    Interval = linspace(minval,maxval,res);
    figure
    plot(Interval,A{1,1})
    %title('Activation Functions')
    xlabel('Stimulation','FontSize',18)
    ylabel('Spike Rate','FontSize',18)

    %%%%Figure 10
    counter1 = 0;
    NumRows = 5;
    NumCols = 2;
    Function = 'Function ';
    %legendLabels3 = cell(NumNeurons, 1);
    figure
    title('Activation functions of first 10 populations')
    for i = 1:5
        for k = 1: 2
        counter1 = counter1 + 1;

        subplot(NumRows,NumCols,counter1)
        title(['Population of Filter ' num2str(i) ' IHC' num2str(k)]);
       % legendLabels3{j} = ['Function ' num2str(k)];
        hold on
        plot(Interval,A{k,i})
        xlabel('Stimulation','FontWeight','bold')
        ylabel('Spike Rate','FontWeight','bold')
        end
    end

    %%Figure 11
    figure
    hold on
    title('DecodedSignal')
    plot(t,DecodedSignal)
    xlabel('Time (seconds)','FontSize',18)
    ylabel('Amplitude','FontSize',18)

    %%%Figure 12
    figure;
    subplot(3,1,1);
    plot(t, SignalIn);
    title('Original Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');

    subplot(3,1,2);
    plot(t, DecodedSignal);
    title('Raw Decoded Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');

subplot(3,1,3);
plot(t, filtered_signal);
title('Low-Pass Filtered Signal');
xlabel('Time (s)');
ylabel('Amplitude');

%Listen to decoded output
soundsc(DecodedSignal,Fs)
soundsc(SignalIn,Fs)
soundsc(filtered_signal,Fs)

audioVector = DecodedSignal;
%audioVector = filtered_signal;
audioVector = audioVector / max(abs(audioVector));
audiowrite('OrangeModelOutput.wav',audioVector,Fs)

    %%Figure 13
figure
hold on
    sgtitle('Simulation 1A','FontWeight','bold')
   
    subplot(3,1,1);
    plot(t, filterSum);
    title('Basilar Mebrane Output');
    xlabel('Time (s)');
    ylabel('Amplitude');
    
    subplot(3,1,2);
    plot(t,IHCSum);
    title('Inner Hair Cell Output');
    xlabel('Time (s)');
    ylabel('Amplitude');

    subplot(3,1,3);
    plot(t, filtered_signal);
    title('Low-Pass Filtered Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');


toc