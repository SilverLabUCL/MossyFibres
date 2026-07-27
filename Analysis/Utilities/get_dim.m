function [varexp,dimmax,varmax] = get_dim(dFF,N_sub,D_max,acquisition_rate,frac_T_train,frac_N_train,N_reps)
    
    if nargin < 4 
        acquisition_rate = [];
    end
    if nargin < 5 || isempty(frac_T_train)
        frac_T_train = 0.8;
    end
    if nargin < 6 || isempty(frac_N_train)
        frac_N_train = 0.8;
    end
    % Number random repeats, including both subsampling population &
    % subsampling time. 
    if nargin < 7 || isempty(N_reps)
        N_reps = 200;
    end
    
    %%
    [N, T] = size(dFF);

    % If N_sub is left empty, use maximum population size
    if isempty(N_sub)
        N_sub = N;
    end
    
    % Cases to throw errors
    if D_max > N_sub
        error('Maximum tested dimensionality cannot be larger than subpopulation size.')
    end

    if N_sub > N
        error('Subpopulation size is too large.')
    end
    
    if ~isempty(acquisition_rate)
        disp('Block shuffle')
    else
        disp('Random')        
    end
    
    %%
    varexp = nan(N_reps,D_max);

    for k = 1:N_reps

        if mod(k/N_reps*10,1) == 0
            disp(['Dimensionality calculation is ',num2str(floor(k/N_reps*10)*10),'% complete.'])
        end
        
        % Randomly subsample the population
        dFF_sub = dFF(randsample(N,N_sub),:);
        
        % For training indices use block shuffled time points (1s default)
        % if no acquisition rate, randomly shuffle
        if ~isempty(acquisition_rate)
            train_ixs = block_shuffle_time(T,acquisition_rate);
            train_ixs = train_ixs(1:round(frac_T_train*T)); 
        else
            train_ixs = randsample(T,round(frac_T_train*T));           
        end
        
        % Randomly sample training population
        ix_x = sort(randsample(N_sub,round(frac_N_train*N_sub)));
        ix_y = setdiff(1:N_sub,ix_x)';

        Neuron_x = dFF_sub(ix_x,:);
        Neuron_y = dFF_sub(ix_y,:);

        [num_dim, ~,~,res,~,~] = peer_predict_dim(Neuron_x,Neuron_y,train_ixs,1:D_max);
        varexp(k,:) = nanmean(res,1);
    end
    
    [varmax,dimmax]=max(nanmean(varexp,1));