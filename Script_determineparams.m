function [ALPHA_BEST, PHI0_BEST, J_BIAS_BEST] = Script_determineparams(inputVoltage,Fs)

    dur = length(inputVoltage); % Simulation duration (sec)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%FIND BEST PARAMS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 1/4: Generating a Population of Neurons
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters for Neurons
    
    niter = 100;
    e_best = 0;
    x_t = inputVoltage;
    for ntor = 1:niter
        if mod(ntor,25)==0
            fprintf('Percent Complete: %3.2f\n',100*(ntor/niter))
        end
    
        ALPHA = 10*rand(1,5);
        PHI0 = ones(1,5);
        J_BIAS = 10*rand(1,5);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate Neurons - Neurogenesis?!? (makelifn)
    
        N1 = makelifn(ALPHA(1),PHI0(1),J_BIAS(1));
        N2 = makelifn(ALPHA(2),PHI0(2),J_BIAS(2));
        N3 = makelifn(ALPHA(3),PHI0(3),J_BIAS(3));
        N4 = makelifn(ALPHA(4),PHI0(4),J_BIAS(4));
        N5 = makelifn(ALPHA(5),PHI0(5),J_BIAS(5));
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Part 2/4: Encoding with a Population of Neurons
    
        TT = zeros(dur,1); % Time Recording
        % inputVoltage = zeros(dur,1); % Stimulus Record
        D1 = zeros(dur,1); % Spike Train Vector #1
        D2 = zeros(dur,1); % Spike Train Vector #2
        D3 = zeros(dur,1); % Spike Train Vector #3
        D4 = zeros(dur,1); % Spike Train Vector #4
        D5 = zeros(dur,1); % Spike Train Vector #5
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Encode Physical Value Into Spike Rates (updatelifn)
    
        % Collect Samples
        for itor = 1:dur
    
            % Record Time Stamp
            TT(itor) = itor/Fs;
    
            % Simulate/Update LIFN (updatelifn)
            N1 = updatelifn(x_t(itor),N1,Fs);
            N2 = updatelifn(x_t(itor),N2,Fs);
            N3 = updatelifn(x_t(itor),N3,Fs);
            N4 = updatelifn(x_t(itor),N4,Fs);
            N5 = updatelifn(x_t(itor),N5,Fs);
    
            % Record Spikes
            D1(itor) = N1.V==1;
            D2(itor) = N2.V==1;
            D3(itor) = N3.V==1;
            D4(itor) = N4.V==1;
            D5(itor) = N5.V==1;
    
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Part 3/4: Probe & Determine Decoders For a Population of Neurons
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parameters for Neural "Probe"
    
        res = 200; % resolution of values implemented in "probing"
        minval = -1.0; % minimum considered value
        maxval = 1.0; % maximum considered value
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Characterize Activation Functions of Population (characterizelifn)
    
        %[A1 X] = characterizelifn(N1,minval,maxval,res,Fs);
        %[A2 X] = characterizelifn(N2,minval,maxval,res,Fs);
        %[A3 X] = characterizelifn(N3,minval,maxval,res,Fs);
        %[A4 X] = characterizelifn(N4,minval,maxval,res,Fs);
        %[A5 X] = characterizelifn(N5,minval,maxval,res,Fs);
    
        [A1 X] = responselifn(N1,minval,maxval,res);
        [A2 X] = responselifn(N2,minval,maxval,res);
        [A3 X] = responselifn(N3,minval,maxval,res);
        [A4 X] = responselifn(N4,minval,maxval,res);
        [A5 X] = responselifn(N5,minval,maxval,res);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compile activation functions into a single matrix
        %   appropriate as input to 'determinedecoders'
    
        A = [A1 A2 A3 A4 A5];
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determine Decoders for Population (determinedecoders)
    
        PHI = determinedecoders(A,X);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Part 4/4: Decode Value Encoded By the Neural Population
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Obtain Temporal Activation Functions (decodespiketrain)
    
        x_hat1 = decodespiketrain(N1,D1,Fs);
        x_hat2 = decodespiketrain(N2,D2,Fs);
        x_hat3 = decodespiketrain(N3,D3,Fs);
        x_hat4 = decodespiketrain(N4,D4,Fs);
        x_hat5 = decodespiketrain(N5,D5,Fs);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compile temporal activation functions into a matrix At
        %   appropriate as input to 'decoderspikerate'
    
        At = [x_hat1 x_hat2 x_hat3 x_hat4 x_hat5];
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Decode Recorded Spike Rates (decodespikerate)
    
        x_t_hat = decodespikerate(At,PHI);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Identify New Params, If Better
        e = corr(inputVoltage,x_t_hat(1:length(TT)));
        if e>e_best
            fprintf('Best Correlation So Far: %3.2f\n',e)
            fprintf('     Found At Iteration: %d\n',ntor)
            e_best = e;
            ALPHA_BEST = ALPHA;
            PHI0_BEST = PHI0;
            J_BIAS_BEST = J_BIAS;
        end
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%VISUALIZE BEST PARAMS%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ALPHA = ALPHA_BEST;
    PHI0 = PHI0_BEST;
    J_BIAS = J_BIAS_BEST;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate Neurons - Neurogenesis?!? (makelifn)
    
    N1 = makelifn(ALPHA(1),PHI0(1),J_BIAS(1));
    N2 = makelifn(ALPHA(2),PHI0(2),J_BIAS(2));
    N3 = makelifn(ALPHA(3),PHI0(3),J_BIAS(3));
    N4 = makelifn(ALPHA(4),PHI0(4),J_BIAS(4));
    N5 = makelifn(ALPHA(5),PHI0(5),J_BIAS(5));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2/4: Encoding with a Population of Neurons
    
    TT = zeros(dur,1); % Time Recording
    inputVoltage = zeros(dur,1); % Stimulus Record
    V1 = zeros(dur,1); % Membrane Voltage Record - Neuron #1
    V2 = zeros(dur,1); % Membrane Voltage Record - Neuron #2
    V3 = zeros(dur,1); % Membrane Voltage Record - Neuron #3
    V4 = zeros(dur,1); % Membrane Voltage Record - Neuron #4
    V5 = zeros(dur,1); % Membrane Voltage Record - Neuron #5
    D1 = zeros(dur,1); % Spike Train Vector #1
    D2 = zeros(dur,1); % Spike Train Vector #2
    D3 = zeros(dur,1); % Spike Train Vector #3
    D4 = zeros(dur,1); % Spike Train Vector #4
    D5 = zeros(dur,1); % Spike Train Vector #5
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Encode Physical Value Into Spike Rates (updatelifn)
    
    % Collect Samples
    for itor = 1:dur
        
        % Record Time Stamp
        TT(itor) = itor/Fs;
        
        % Simulate/Update LIFN (updatelifn)
        N1 = updatelifn(x_t(itor),N1,Fs);
        N2 = updatelifn(x_t(itor),N2,Fs);
        N3 = updatelifn(x_t(itor),N3,Fs);
        N4 = updatelifn(x_t(itor),N4,Fs);
        N5 = updatelifn(x_t(itor),N5,Fs);
        
        % Record Membrane Voltage
        V1(itor) = N1.V;
        V2(itor) = N2.V;
        V3(itor) = N3.V;
        V4(itor) = N4.V;
        V5(itor) = N5.V;
    
        % Record Spikes
        D1(itor) = N1.V==1;
        D2(itor) = N2.V==1;
        D3(itor) = N3.V==1;
        D4(itor) = N4.V==1;
        D5(itor) = N5.V==1;
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Visualize Stimulus
    
    % figure
    % plot(TT,inputVoltage,'linewidth',2);
    % xlabel('Time (sec)','fontsize',18)
    % ylabel('Magnitude','fontsize',18)
    % %axis([0 1 -0.2 0.2])
    % axis([0 1 0 1])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Visualize Voltage
    
    % figure
    % plot(TT,V1,'linewidth',2);
    % hold all
    % plot(TT,V2,'linewidth',2);
    % plot(TT,V3,'linewidth',2);
    % plot(TT,V4,'linewidth',2);
    % plot(TT,V5,'linewidth',2);
    % xlabel('Time (sec)','fontsize',18)
    % ylabel('Voltage','fontsize',18)
    % axis([0 0.02 0 1])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Visualize Spikes
    
    % figure
    % plot(TT,D1,'linewidth',2);
    % hold all
    % plot(TT,D2,'linewidth',2);
    % plot(TT,D3,'linewidth',2);
    % plot(TT,D4,'linewidth',2);
    % plot(TT,D5,'linewidth',2);
    % xlabel('Time (sec)','fontsize',18)
    % ylabel('Voltage','fontsize',18)
    % axis([0 0.02 0 1])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 3/4: Probe & Determine Decoders For a Population of Neurons
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters for Neural "Probe"
    
    res = 200; % resolution of values implemented in "probing"
    minval = 0.0; % minimum considered value
    maxval = 1.0; % maximum considered value
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Characterize Activation Functions of Population (characterizelifn)
    
    [A1 X] = characterizelifn(N1,minval,maxval,res,Fs);
    [A2 X] = characterizelifn(N2,minval,maxval,res,Fs);
    [A3 X] = characterizelifn(N3,minval,maxval,res,Fs);
    [A4 X] = characterizelifn(N4,minval,maxval,res,Fs);
    [A5 X] = characterizelifn(N5,minval,maxval,res,Fs);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compile activation functions into a single matrix
    %   appropriate as input to 'determinedecoders'
    
    A = [A1 A2 A3 A4 A5];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine Decoders for Population (determinedecoders)
    
    PHI = determinedecoders(A,X);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Visualize Neuron Response Functions
    
    % figure
    % plot(X,A,'k','linewidth',2)
    % axis([minval maxval 0 600])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 4/4: Decode Value Encoded By the Neural Population
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Obtain Temporal Activation Functions (decodespiketrain)
    
    x_hat1 = decodespiketrain(N1,D1,Fs);
    x_hat2 = decodespiketrain(N2,D2,Fs);
    x_hat3 = decodespiketrain(N3,D3,Fs);
    x_hat4 = decodespiketrain(N4,D4,Fs);
    x_hat5 = decodespiketrain(N5,D5,Fs);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compile temporal activation functions into a matrix At
    %   appropriate as input to 'decoderspikerate'
    
    At = [x_hat1 x_hat2 x_hat3 x_hat4 x_hat5];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Decode Recorded Spike Rates (decodespikerate)
    
    x_t_hat = decodespikerate(At,PHI);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Visualize Temporal Activation Functions
    
    % figure
    % plot(TT,At(1:length(TT),:))
    % xlabel('Time (sec)','fontsize',18);
    % ylabel('x(t)','fontsize',18);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Visualize Decoding
    
    % figure
    % plot(TT,inputVoltage,'linewidth',2)
    % hold all
    % plot(TT,x_t_hat(1:length(TT)).*100,':')
    % xlabel('Time (sec)','fontsize',18);
    % ylabel('x(t)','fontsize',18);
    % %axis([0 1 -0.5 0.5])
    % axis([0 1 -1 2])
    
    fprintf('Best Parameter Values Are:\n')
    ALPHA
    PHI0
    J_BIAS

    ALPHA_BEST = ALPHA;
    PHI0_BEST = PHI0;
    J_BIAS_BEST = J_BIAS;
    
    %eof
end