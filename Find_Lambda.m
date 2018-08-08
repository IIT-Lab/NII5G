function [lambda_min, lambda_max, lambda, User] = Find_Lambda(L,Chn,P_max_total_dBm,User,Cells,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico,Power_Mode,...
    per_BS_power_constraint,Num_Rx_Ant,noise_dB,intf_2nd,per_BS_power_constraint_matrix)

%%%---------------------------------------------------------------------
%%% This function finds the upper and lower bound of lambda for which the
%%% power constraint meets and breaks
%%%---------------------------------------------------------------------

Total_BS = L_Macro + L_Pico;
lambda_min = 0;
lambda_max = 10;
lambda = zeros(L*Total_BS,1);

%lambda = ones(L*Total_BS,1)*1e3;

switch Power_Mode
    
    %%
    case 1
        %if Power_Mode == 1 %sum power constraint
        
        User = BF_wmmse(L,User,Cells,Chn,lambda_min,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico);
        sumpower = Sum_Power(L,Cells,User,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico);
        
        if dec2dB(sumpower) < P_max_total_dBm
            
            while dec2dB(sumpower) < P_max_total_dBm    %keep decreasing lambda_min until the total power exceeds
                
                lambda_max = lambda_min;
                lambda_min = lambda_min - 10;
                
                User = BF_wmmse(L,User,Cells,Chn,lambda_min,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico);
                sumpower = Sum_Power(L,Cells,User,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico);
            end
            
        else
            User = BF_wmmse(L,User,Cells,Chn,lambda_max,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico);
            sumpower = Sum_Power(L,Cells,User,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico);
            
            while dec2dB(sumpower) > P_max_total_dBm %keep increasing lambda_max until the total power less than the budget
                
                lambda_min = lambda_max;
                
                lambda_max = lambda_max * 2;
                
                User = BF_wmmse(L,User,Cells,Chn,lambda_max,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico);
                sumpower = Sum_Power(L,Cells,User,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico);
                
            end
        end
        fprintf('$$$$$$$$$  lambda_min = %f   makes sum power deficient $$$$$$$$$ \n', lambda_min);
        fprintf('$$$$$$$$$  lambda_max = %f   makes sum power sufficient $$$$$$$$$ \n', lambda_max);
        
        %%
    case 2
        %if Power_Mode == 2 %per BS power constraint
        
        sub_iter = 1;
        previous_wgt_sumrate = 1;
        while 1 > 0
            
            [User, max_diff,min_norm,max_norm] = BF_wmmse(L,User,Cells,Chn,lambda,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico);
            [User, sumrate, curr_wgt_sumrate, SumBackhaul] = rate_DL(L,Num_Rx_Ant,User,Cells,Chn,noise_dB,intf_2nd,L_Macro,L_Pico,Num_TxAnt_Pico);
            [~,per_BS_power] = Sum_Power(L,Cells,User,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico);
            
            power_gap = (per_BS_power - per_BS_power_constraint)./per_BS_power_constraint;
            max_power_gap = max(power_gap);
            
            wgt_sumrate_diff = abs(curr_wgt_sumrate - previous_wgt_sumrate)/curr_wgt_sumrate;
            %        step_size = max_norm/sub_iter;
            %step_size = 10/norm(power_gap);
            if sub_iter < 500
                step_size = 1e3;
            end
            %         else
            %             step_size = min([max_norm*1e2 1e5])/sqrt(sub_iter);
            %         end
            % update lambda using subgradient method
            %      if (max_diff > 1e-2)  || (max_power_gap > 1e-2)
            lambda_new = max(lambda + step_size*(per_BS_power - per_BS_power_constraint),0);
            %      end
            %            fprintf('$$$$$ max BF diff = %f $$$$$\n',max_diff);
            %            fprintf('$$$$$ max power gap = %f $$$$$\n',max((per_BS_power - per_BS_power_constraint)./per_BS_power_constraint));
            
            if norm(lambda_new) == 0 || (wgt_sumrate_diff < 1e-2  &&  abs(max_power_gap) < 2e-2)
                %|| (norm(lambda_new - lambda)/norm(lambda_new) < 1e-2 && wgt_sumrate_diff < 1e-2  &&  max_power_gap < 2e-2 )
                break;
            else
                if mod(sub_iter,1000) == 0
                    fprintf('$$$$$ subgradient method to update lambda iter = %d $$$$$$\n',sub_iter);
                    fprintf('$$$$$ wgt sum rate diff = %f $$$$$\n',wgt_sumrate_diff);
                    fprintf('$$$$$ max power gap = %f $$$$$\n',max_power_gap);
                    fprintf('$$$$$ lambda diff = %f $$$$$\n',norm(lambda_new - lambda)/norm(lambda_new));
                end
                lambda = lambda_new;
                previous_wgt_sumrate = curr_wgt_sumrate;
            end
            
            if sub_iter >= 1e4
                break;
            end
            
            sub_iter = sub_iter + 1;
        end
        
        
        %%
    case 21
        %elseif Power_Mode == 21   %per BS power constraint using gradient descent method to find the optimal dual variable
        lambda = zeros(L,Total_BS); %initial dual variable for PowerMode = 21 using gradient descent algorithm
        
        %compute current dual obj
        [curr_UL_obj,User] = Compute_UL_obj(L,User,Cells,lambda,Chn,L_Macro,L_Pico,Num_TxAnt_Pico,per_BS_power_constraint_matrix);
        
        iter = 1;
        t = 1e-4;
        while 1>0
            %%% compute KZ_UL for each user for current lambda
            %     for l = 1:L
            %         for ik = 1:length(Cells(l).Scheduled_User)
            %             k = Cells(l).Scheduled_User(ik);
            %             User(l,k).Kz_UL = User(l,k).Kz_UL_intf;
            %             for iBS = 1:length(User(l,k).ServingCluster)
            %                 curr_cell = ceil(User(l,k).ServingCluster(iBS)/Total_BS);
            %                 curr_BS = User(l,k).ServingCluster(iBS) - (curr_cell-1)*Total_BS;
            %
            %                 User(l,k).Kz_UL = User(l,k).Kz_UL + lambda(curr_cell,curr_BS)*User(l,k).CoopMatrix(:,:,iBS);
            %             end
            %         end
            %     end
            
            %  delta_lambda = per_BS_power_constraint_matrix; %descent direction of lambda
            delta_lambda = zeros(L,Total_BS);
            %%% compute the descent direction of lambda
            gradient_lambda = per_BS_power_constraint_matrix;
            hessian_lambda = zeros(L,Total_BS);
            for l = 1:L
                for ik = 1:length(Cells(l).Scheduled_User)
                    k = Cells(l).Scheduled_User(ik);
                    H = (User(l,k).beam_rx)' * Get_Chn(Chn, User(l,k).ServingCluster, l, k, L_Macro,L_Pico,Num_TxAnt_Pico)';
                    for iBS = 1:length(User(l,k).ServingCluster)
                        curr_cell = ceil(User(l,k).ServingCluster(iBS)/Total_BS);
                        curr_BS = User(l,k).ServingCluster(iBS) - (curr_cell-1)*Total_BS;
                        
                        gradient = - real( (User(l,k).wgt_mse)^2 * H * pinv(User(l,k).Kz_UL) * User(l,k).ServeBS(iBS).CoopMatrix * pinv(User(l,k).Kz_UL) * H' );% + ...
                        %      2*real( (User(l,k).wgt_mse)^2 * H * pinv(User(l,k).Kz_UL) * pinv(User(l,k).Kz_UL) * User(l,k).ServeBS(iBS).CoopMatrix * (eye(User(l,k).ServingAnt) - User(l,k).Kz_UL * pinv(User(l,k).Kz_UL) ) * H' );
                        
                        gradient_lambda(curr_cell,curr_BS) = gradient_lambda(curr_cell,curr_BS) + gradient;
                        hessian_lambda(curr_cell,curr_BS) = hessian_lambda(curr_cell,curr_BS) + ...
                            real(2*(User(l,k).wgt_mse)^2*H*pinv(User(l,k).Kz_UL)*User(l,k).ServeBS(iBS).CoopMatrix*inv(User(l,k).Kz_UL + eye(User(l,k).ServingAnt))*User(l,k).ServeBS(iBS).CoopMatrix*pinv(User(l,k).Kz_UL)*H');
                        
                        %                   Kz_UL_delta = User(l,k).Kz_UL + t*User(l,k).ServeBS(iBS).CoopMatrix;
                        %                   gradient_delta = - real((User(l,k).wgt_mse)^2*H*pinv(Kz_UL_delta)*User(l,k).ServeBS(iBS).CoopMatrix*pinv(Kz_UL_delta)*H') - gradient;% + ...
                        %               %          2*real( (User(l,k).wgt_mse)^2*H*pinv(Kz_UL_delta)*pinv(Kz_UL_delta)*User(l,k).ServeBS(iBS).CoopMatrix*(eye(User(l,k).ServingAnt) - Kz_UL_delta*pinv(Kz_UL_delta))*H') - gradient;
                        %                     numer_hessian = max(gradient_delta/t,0);
                        %                     hessian_lambda(curr_cell,curr_BS) = hessian_lambda(curr_cell,curr_BS) + numer_hessian;
                    end
                    
                end
            end
            
            
            if norm(lambda) == 0 && max(-gradient_lambda./per_BS_power_constraint_matrix) < 0    %redundant power constraint
                fprintf('$$$ redundant power constraint max power gap is %f $$$ \n',max(-gradient_lambda./per_BS_power_constraint_matrix));
                break;
            end
            
            if norm(lambda) ~= 0 && abs(max(-gradient_lambda./per_BS_power_constraint_matrix)) < 1e-2
                fprintf('$$$ max power gap is %f at iter %d $$$ \n',max(-gradient_lambda./per_BS_power_constraint_matrix),iter);
                break;
            end
            
            
            for l = 1:L
                for iBS = 1:Total_BS
                    if hessian_lambda(l,iBS) <= 1e-7
                        hessian_lambda(l,iBS) = 1e-7;
                    end
                    delta_lambda(l,iBS) = - gradient_lambda(l,iBS)/hessian_lambda(l,iBS);
                    
                    while abs(delta_lambda(l,iBS)) > 1e2
                        delta_lambda(l,iBS) = delta_lambda(l,iBS)/10;
                    end
                    
                end
            end
            
            %%% backtrack line search to find proper step_size
            
            %   step_size = 1;
            %1/sqrt(iter)*10^(-floor(iter/1000));
            
            %    if iter>1e3
            step_size = 1/sqrt(iter)*10^(-floor(iter/1000));
            %    end
            
            alpha = 0.3;
            beta = 0.5;
            
            %       while 1>0   %backtrack line search
            new_lambda = max(lambda + step_size*delta_lambda,0);
            %             if sum(isfinite(new_lambda)) ~= size(new_lambda,2);
            %                 new_lambda
            %             end
            [new_UL_obj,User] = Compute_UL_obj(L,User,Cells,new_lambda,Chn,L_Macro,L_Pico,Num_TxAnt_Pico,per_BS_power_constraint_matrix);
            
            %           if new_UL_obj > curr_UL_obj
            %+ alpha*step_size*trace(gradient_lambda'*delta_lambda);
            %               step_size = step_size*beta;
            %           else
            %               fprintf('$$$ backtrack step size = %f $$$ \n',step_size);
            %               break;
            %           end
            %       end
            
            
            %         if   abs(new_UL_obj - curr_UL_obj)/curr_UL_obj < 1e-2 && max(-gradient_lambda./per_BS_power_constraint_matrix) < 2e-2
            %             %(norm(new_lambda) <= 1e-10 || norm(lambda-new_lambda)/norm(new_lambda) < 1e-2) && max(-gradient_lambda./per_BS_power_constraint_matrix) < 2e-2
            %
            %             %        if sum(gradient_lambda.*gradient_lambda./hessian_lambda) < 1e-2
            %             %          lambda = new_lambda;
            %             fprintf('$$$ max power gap from gradient is %f $$$ \n',max(-gradient_lambda./per_BS_power_constraint_matrix));
            %             [~,User] = Compute_UL_obj(L,User,Cells,lambda,Chn,L_Macro,L_Pico,Num_TxAnt_Pico,per_BS_power_constraint_matrix);
            %             break;
            %         else
            lambda = new_lambda;
            curr_UL_obj = new_UL_obj;
            %        end
            
            iter = iter+1;
            %         if mod(iter,1000) == 0
            %             iter
            %         end
            
            if iter > 1e4
                iter
                break;
            end
            
        end
        
        
    case 3
        % if Power_Mode == 3, two total power constraint, one for Macro
        % BS's the other for Pico BS's
        
        sub_enable = 0;
        
        lambda = zeros(L,Total_BS);
        
        %compute current dual obj
        [curr_UL_obj,User] = Compute_UL_obj(L,User,Cells,lambda,Chn,L_Macro,L_Pico,Num_TxAnt_Pico,per_BS_power_constraint_matrix);
        
        %%% compute the descent direction of lambda, which is also the
        %%% subgradient
        gradient_lambda = per_BS_power_constraint_matrix;
        for l = 1:L
            for ik = 1:length(Cells(l).Scheduled_User)
                k = Cells(l).Scheduled_User(ik);
                H = (User(l,k).beam_rx)' * Get_Chn(Chn, User(l,k).ServingCluster, l, k, L_Macro,L_Pico,Num_TxAnt_Pico)';
                for iBS = 1:length(User(l,k).ServingCluster)
                    curr_cell = ceil(User(l,k).ServingCluster(iBS)/Total_BS);
                    curr_BS = User(l,k).ServingCluster(iBS) - (curr_cell-1)*Total_BS;
                    
                    gradient = - real( (User(l,k).wgt_mse)^2 * H * pinv(User(l,k).Kz_UL) * User(l,k).ServeBS(iBS).CoopMatrix * pinv(User(l,k).Kz_UL) * H' );
                    gradient_lambda(curr_cell,curr_BS) = gradient_lambda(curr_cell,curr_BS) + gradient;
                    
                end
            end
        end
        
        %check sum power constraint for zero lambdas
        if L == 1 % single cell
            
            if gradient_lambda(1)>0 && sum(gradient_lambda(2:end)) > 0 %&& sub_enable
                
                fprintf('$$$ Two Total Power constraints are both redundant $$$ \n');
                
            elseif gradient_lambda(1)>0  %&& sub_enable
                
                fprintf('$$$ Sum Macro Power is met, only do bisection over Pico Lambda $$$ \n');
                
                % look for the bound of pico lambda's
                max_pico_lambda = 10;
                min_pico_lambda = 0;
                while 1>0
                    lambda(2:end) = max_pico_lambda*ones(1,L_Pico);
                    %compute current dual obj
                    [curr_UL_obj,User] = Compute_UL_obj(L,User,Cells,lambda,Chn,L_Macro,L_Pico,Num_TxAnt_Pico,per_BS_power_constraint_matrix);
                    %check if new lambda will make Sum Power for Pico's
                    %satisfied
                    gradient_lambda = per_BS_power_constraint_matrix;
                    for l = 1:L
                        for ik = 1:length(Cells(l).Scheduled_User)
                            k = Cells(l).Scheduled_User(ik);
                            H = (User(l,k).beam_rx)' * Get_Chn(Chn, User(l,k).ServingCluster, l, k, L_Macro,L_Pico,Num_TxAnt_Pico)';
                            for iBS = 1:length(User(l,k).ServingCluster)
                                curr_cell = ceil(User(l,k).ServingCluster(iBS)/Total_BS);
                                curr_BS = User(l,k).ServingCluster(iBS) - (curr_cell-1)*Total_BS;
                                
                                gradient = - real( (User(l,k).wgt_mse)^2 * H * pinv(User(l,k).Kz_UL) * User(l,k).ServeBS(iBS).CoopMatrix * pinv(User(l,k).Kz_UL) * H' );
                                gradient_lambda(curr_cell,curr_BS) = gradient_lambda(curr_cell,curr_BS) + gradient;
                                
                            end
                        end
                    end
                    
                    if gradient_lambda(1) < 0
                        fprintf('$$$ Macro Power is NOT right now !!! $$$\n');
                 %       break;
                    end
                    
                    if sum(gradient_lambda(2:end)) < 0
                        min_pico_lambda = max_pico_lambda;
                        max_pico_lambda = max_pico_lambda*2;
                    else
                        break;
                    end
                end
                
                %find the optimal lambda for Pico BS's
                while 1>0
                    lambda(2:end) = (max_pico_lambda + min_pico_lambda)/2*ones(1,L_Pico);
                    %compute current dual obj
                    [curr_UL_obj,User] = Compute_UL_obj(L,User,Cells,lambda,Chn,L_Macro,L_Pico,Num_TxAnt_Pico,per_BS_power_constraint_matrix);
                    %check if new lambda will make Sum Power for Pico's
                    %satisfied
                    gradient_lambda = per_BS_power_constraint_matrix;
                    for l = 1:L
                        for ik = 1:length(Cells(l).Scheduled_User)
                            k = Cells(l).Scheduled_User(ik);
                            H = (User(l,k).beam_rx)' * Get_Chn(Chn, User(l,k).ServingCluster, l, k, L_Macro,L_Pico,Num_TxAnt_Pico)';
                            for iBS = 1:length(User(l,k).ServingCluster)
                                curr_cell = ceil(User(l,k).ServingCluster(iBS)/Total_BS);
                                curr_BS = User(l,k).ServingCluster(iBS) - (curr_cell-1)*Total_BS;
                                
                                gradient = - real( (User(l,k).wgt_mse)^2 * H * pinv(User(l,k).Kz_UL) * User(l,k).ServeBS(iBS).CoopMatrix * pinv(User(l,k).Kz_UL) * H' );
                                gradient_lambda(curr_cell,curr_BS) = gradient_lambda(curr_cell,curr_BS) + gradient;
                                
                            end
                        end
                    end
                    
                    if gradient_lambda(1) < 0
                        fprintf('$$$ Macro Power is NOT right now !!! $$$\n');
                    %    break;
                    end
                    
   %                 fprintf('$$$ Pico Power Gap = %f  $$$\n',sum(gradient_lambda(2:end))/sum(per_BS_power_constraint_matrix(2:end)));
                    
                    if abs(sum(gradient_lambda(2:end)))/sum(per_BS_power_constraint_matrix(2:end)) < 2e-2
                        break;
                    end
                    
                    if sum(gradient_lambda(2:end)) < 0
                        min_pico_lambda = (max_pico_lambda + min_pico_lambda)/2;
                    else
                        max_pico_lambda = (max_pico_lambda + min_pico_lambda)/2;
                    end
                end
                
            elseif sum(gradient_lambda(2:end)) > 0   %&&  sub_enable
                
                fprintf('$$$ Sum Pico Power is met, only do bisection over Macro Lambda $$$\n');
                % look for the bound of Macro lambda
                max_Macro_lambda = 10;
                min_Macro_lambda = 0;
                while 1>0
                    lambda(1) = max_Macro_lambda;
                    %compute current dual obj
                    [curr_UL_obj,User] = Compute_UL_obj(L,User,Cells,lambda,Chn,L_Macro,L_Pico,Num_TxAnt_Pico,per_BS_power_constraint_matrix);
                    %check if new lambda will make Sum Power for Pico's
                    %satisfied
                    gradient_lambda = per_BS_power_constraint_matrix;
                    for l = 1:L
                        for ik = 1:length(Cells(l).Scheduled_User)
                            k = Cells(l).Scheduled_User(ik);
                            H = (User(l,k).beam_rx)' * Get_Chn(Chn, User(l,k).ServingCluster, l, k, L_Macro,L_Pico,Num_TxAnt_Pico)';
                            for iBS = 1:length(User(l,k).ServingCluster)
                                curr_cell = ceil(User(l,k).ServingCluster(iBS)/Total_BS);
                                curr_BS = User(l,k).ServingCluster(iBS) - (curr_cell-1)*Total_BS;
                                
                                gradient = - real( (User(l,k).wgt_mse)^2 * H * pinv(User(l,k).Kz_UL) * User(l,k).ServeBS(iBS).CoopMatrix * pinv(User(l,k).Kz_UL) * H' );
                                gradient_lambda(curr_cell,curr_BS) = gradient_lambda(curr_cell,curr_BS) + gradient;
                                
                            end
                        end
                    end
                    
                    if sum(gradient_lambda(2:end)) < 0
                        fprintf('$$$ Pico Power is NOT right now !!! $$$\n');
                    %    break;
                    end
                    
                    if gradient_lambda(1) < 0
                        min_Macro_lambda = max_Macro_lambda;
                        max_Macro_lambda = max_Macro_lambda*2;
                    else
                        break;
                    end
                end
                
                %find the optimal lambda for Macro BS
                while 1>0
                    lambda(1) = (max_Macro_lambda + min_Macro_lambda)/2;
                    %compute current dual obj
                    [curr_UL_obj,User] = Compute_UL_obj(L,User,Cells,lambda,Chn,L_Macro,L_Pico,Num_TxAnt_Pico,per_BS_power_constraint_matrix);
                    %check if new lambda will make Sum Power for Pico's
                    %satisfied
                    gradient_lambda = per_BS_power_constraint_matrix;
                    for l = 1:L
                        for ik = 1:length(Cells(l).Scheduled_User)
                            k = Cells(l).Scheduled_User(ik);
                            H = (User(l,k).beam_rx)' * Get_Chn(Chn, User(l,k).ServingCluster, l, k, L_Macro,L_Pico,Num_TxAnt_Pico)';
                            for iBS = 1:length(User(l,k).ServingCluster)
                                curr_cell = ceil(User(l,k).ServingCluster(iBS)/Total_BS);
                                curr_BS = User(l,k).ServingCluster(iBS) - (curr_cell-1)*Total_BS;
                                
                                gradient = - real( (User(l,k).wgt_mse)^2 * H * pinv(User(l,k).Kz_UL) * User(l,k).ServeBS(iBS).CoopMatrix * pinv(User(l,k).Kz_UL) * H' );
                                gradient_lambda(curr_cell,curr_BS) = gradient_lambda(curr_cell,curr_BS) + gradient;
                                
                            end
                        end
                    end
                    
                    if sum(gradient_lambda(2:end)) < 0
                        fprintf('$$$ Pico Power is NOT right now !!! $$$\n');
                     %   break;
                    end
                    
    %                fprintf('$$$ Macro Power Gap = %f  $$$\n',sum(gradient_lambda(1))/sum(per_BS_power_constraint_matrix(1)));
                    
                    if abs(sum(gradient_lambda(1)))/sum(per_BS_power_constraint_matrix(1)) < 2e-2
                        break;
                    end
                    
                    if sum(gradient_lambda(1)) < 0
                        min_Macro_lambda = (max_Macro_lambda + min_Macro_lambda)/2;
                    else
                        max_Macro_lambda = (max_Macro_lambda + min_Macro_lambda)/2;
                    end
                end
                
            else
         %       fprintf('$$$ No idea which lambda should be zero $$$ \n');
         %       fprintf('$$$ Try Pico Lambda be zero and bisection over Macro Lambda \n');
                
                % look for the bound of Macro lambda
                max_Macro_lambda = 10;
                min_Macro_lambda = 0;
                while 1>0 && sub_enable 
                    lambda(1) = max_Macro_lambda;
                    %compute current dual obj
                    [curr_UL_obj,User] = Compute_UL_obj(L,User,Cells,lambda,Chn,L_Macro,L_Pico,Num_TxAnt_Pico,per_BS_power_constraint_matrix);
                    %check if new lambda will make Sum Power for Pico's
                    %satisfied
                    gradient_lambda = per_BS_power_constraint_matrix;
                    for l = 1:L
                        for ik = 1:length(Cells(l).Scheduled_User)
                            k = Cells(l).Scheduled_User(ik);
                            H = (User(l,k).beam_rx)' * Get_Chn(Chn, User(l,k).ServingCluster, l, k, L_Macro,L_Pico,Num_TxAnt_Pico)';
                            for iBS = 1:length(User(l,k).ServingCluster)
                                curr_cell = ceil(User(l,k).ServingCluster(iBS)/Total_BS);
                                curr_BS = User(l,k).ServingCluster(iBS) - (curr_cell-1)*Total_BS;
                                
                                gradient = - real( (User(l,k).wgt_mse)^2 * H * pinv(User(l,k).Kz_UL) * User(l,k).ServeBS(iBS).CoopMatrix * pinv(User(l,k).Kz_UL) * H' );
                                gradient_lambda(curr_cell,curr_BS) = gradient_lambda(curr_cell,curr_BS) + gradient;
                                
                            end
                        end
                    end
                    
                    if gradient_lambda(1) < 0
                        min_Macro_lambda = max_Macro_lambda;
                        max_Macro_lambda = max_Macro_lambda*2;
                    else
                        break;
                    end
                    
                    if max_Macro_lambda > 1e10
                        fprintf('$$$ max Macro Lambda is too large $$$ \n');
                        return;
                    end
                    
                end
                
                %find the optimal lambda for Macro BS for fixed Pico Lambda
                while 1>0  && sub_enable 
                    lambda(1) = (max_Macro_lambda + min_Macro_lambda)/2;
                    %compute current dual obj
                    [curr_UL_obj,User] = Compute_UL_obj(L,User,Cells,lambda,Chn,L_Macro,L_Pico,Num_TxAnt_Pico,per_BS_power_constraint_matrix);
                    %check if new lambda will make Sum Power for Pico's
                    %satisfied
                    gradient_lambda = per_BS_power_constraint_matrix;
                    for l = 1:L
                        for ik = 1:length(Cells(l).Scheduled_User)
                            k = Cells(l).Scheduled_User(ik);
                            H = (User(l,k).beam_rx)' * Get_Chn(Chn, User(l,k).ServingCluster, l, k, L_Macro,L_Pico,Num_TxAnt_Pico)';
                            for iBS = 1:length(User(l,k).ServingCluster)
                                curr_cell = ceil(User(l,k).ServingCluster(iBS)/Total_BS);
                                curr_BS = User(l,k).ServingCluster(iBS) - (curr_cell-1)*Total_BS;
                                
                                gradient = - real( (User(l,k).wgt_mse)^2 * H * pinv(User(l,k).Kz_UL) * User(l,k).ServeBS(iBS).CoopMatrix * pinv(User(l,k).Kz_UL) * H' );
                                gradient_lambda(curr_cell,curr_BS) = gradient_lambda(curr_cell,curr_BS) + gradient;
                                
                            end
                        end
                    end
                    
     %               fprintf('$$$ Macro Power Gap = %f  $$$\n',sum(gradient_lambda(1))/sum(per_BS_power_constraint_matrix(1)));
                    
                    if abs(sum(gradient_lambda(1)))/sum(per_BS_power_constraint_matrix(1)) < 2e-2
                        break;
                    end
                    
                    if sum(gradient_lambda(1)) < 0
                        min_Macro_lambda = (max_Macro_lambda + min_Macro_lambda)/2;
                    else
                        max_Macro_lambda = (max_Macro_lambda + min_Macro_lambda)/2;
                    end
                end
                
                if sum(gradient_lambda(2:end))>0  && sub_enable %optimal Pico Lambda is indeed zero
                    fprintf('$$$ optimal Pico Lambda is indeed zero $$$\n');
                    return;
                else
                  %  fprintf('$$$ optimal Pico Lambda is NOT zero$$$\n');
                  %  fprintf('$$$  Try Macro Lambda be zero and bisection over Pico Lambda  $$$\n');
                    lambda = zeros(1,Total_BS);
                    
                    % look for the bound of pico lambda's
                    max_pico_lambda = 10;
                    min_pico_lambda = 0;
                    while 1>0 && sub_enable
                        lambda(2:end) = max_pico_lambda*ones(1,L_Pico);
                        %compute current dual obj
                        [curr_UL_obj,User] = Compute_UL_obj(L,User,Cells,lambda,Chn,L_Macro,L_Pico,Num_TxAnt_Pico,per_BS_power_constraint_matrix);
                        %check if new lambda will make Sum Power for Pico's
                        %satisfied
                        gradient_lambda = per_BS_power_constraint_matrix;
                        for l = 1:L
                            for ik = 1:length(Cells(l).Scheduled_User)
                                k = Cells(l).Scheduled_User(ik);
                                H = (User(l,k).beam_rx)' * Get_Chn(Chn, User(l,k).ServingCluster, l, k, L_Macro,L_Pico,Num_TxAnt_Pico)';
                                for iBS = 1:length(User(l,k).ServingCluster)
                                    curr_cell = ceil(User(l,k).ServingCluster(iBS)/Total_BS);
                                    curr_BS = User(l,k).ServingCluster(iBS) - (curr_cell-1)*Total_BS;
                                    
                                    gradient = - real( (User(l,k).wgt_mse)^2 * H * pinv(User(l,k).Kz_UL) * User(l,k).ServeBS(iBS).CoopMatrix * pinv(User(l,k).Kz_UL) * H' );
                                    gradient_lambda(curr_cell,curr_BS) = gradient_lambda(curr_cell,curr_BS) + gradient;
                                    
                                end
                            end
                        end
                        
                        if sum(gradient_lambda(2:end)) < 0
                            min_pico_lambda = max_pico_lambda;
                            max_pico_lambda = max_pico_lambda*2;
                        else
                            break;
                        end
                        
                        if max_pico_lambda > 1e10
                            fprintf('$$$ max pico lambda is too large $$$\n');
                            return;
                        end
                        
                    end
                    
                    %find the optimal lambda for Pico BS's
                    while 1>0 &&  sub_enable
                        lambda(2:end) = (max_pico_lambda + min_pico_lambda)/2*ones(1,L_Pico);
                        %compute current dual obj
                        [curr_UL_obj,User] = Compute_UL_obj(L,User,Cells,lambda,Chn,L_Macro,L_Pico,Num_TxAnt_Pico,per_BS_power_constraint_matrix);
                        %check if new lambda will make Sum Power for Pico's
                        %satisfied
                        gradient_lambda = per_BS_power_constraint_matrix;
                        for l = 1:L
                            for ik = 1:length(Cells(l).Scheduled_User)
                                k = Cells(l).Scheduled_User(ik);
                                H = (User(l,k).beam_rx)' * Get_Chn(Chn, User(l,k).ServingCluster, l, k, L_Macro,L_Pico,Num_TxAnt_Pico)';
                                for iBS = 1:length(User(l,k).ServingCluster)
                                    curr_cell = ceil(User(l,k).ServingCluster(iBS)/Total_BS);
                                    curr_BS = User(l,k).ServingCluster(iBS) - (curr_cell-1)*Total_BS;
                                    
                                    gradient = - real( (User(l,k).wgt_mse)^2 * H * pinv(User(l,k).Kz_UL) * User(l,k).ServeBS(iBS).CoopMatrix * pinv(User(l,k).Kz_UL) * H' );
                                    gradient_lambda(curr_cell,curr_BS) = gradient_lambda(curr_cell,curr_BS) + gradient;
                                    
                                end
                            end
                        end
                        
        %                fprintf('$$$ Pico Power Gap = %f  $$$\n',sum(gradient_lambda(2:end))/sum(per_BS_power_constraint_matrix(2:end)));
                        
                        if abs(sum(gradient_lambda(2:end)))/sum(per_BS_power_constraint_matrix(2:end)) < 2e-2
                            break;
                        end
                        
                        if sum(gradient_lambda(2:end)) < 0
                            min_pico_lambda = (max_pico_lambda + min_pico_lambda)/2;
                        else
                            max_pico_lambda = (max_pico_lambda + min_pico_lambda)/2;
                        end
                    end
                    
                    if gradient_lambda(1)/per_BS_power_constraint_matrix(1) > 0 && sub_enable% Macro Lambda is indeed zero
                        fprintf('$$$ Macro Lambda is indeed zero $$$ \n');
                        return;
                    else
                      %  fprintf('$$$ Both Lambds are not zero $$$$ \n');
                        fprintf('$$$ Do Subgradient $$$ \n');
                        
                        sub_iter = 1;
                        lambda = Cells(1).initial_lambda;
                        Total_power_const = [per_BS_power_constraint_matrix(1) sum(per_BS_power_constraint_matrix(2:end))];
                        while 1>0
                            step_size = 100;
                            if sub_iter > 1e4
                                step_size = 10;
                            end
                            %compute current dual obj
                            [curr_UL_obj,User] = Compute_UL_obj(L,User,Cells,lambda,Chn,L_Macro,L_Pico,Num_TxAnt_Pico,per_BS_power_constraint_matrix);
                            %check if new lambda will make Sum Power for Pico's
                            %satisfied
                            gradient_lambda = per_BS_power_constraint_matrix;
                            for l = 1:L
                                for ik = 1:length(Cells(l).Scheduled_User)
                                    k = Cells(l).Scheduled_User(ik);
                                    H = (User(l,k).beam_rx)' * Get_Chn(Chn, User(l,k).ServingCluster, l, k, L_Macro,L_Pico,Num_TxAnt_Pico)';
                                    for iBS = 1:length(User(l,k).ServingCluster)
                                        curr_cell = ceil(User(l,k).ServingCluster(iBS)/Total_BS);
                                        curr_BS = User(l,k).ServingCluster(iBS) - (curr_cell-1)*Total_BS;
                                        
                                        gradient = - real( (User(l,k).wgt_mse)^2 * H * pinv(User(l,k).Kz_UL) * User(l,k).ServeBS(iBS).CoopMatrix * pinv(User(l,k).Kz_UL) * H' );
                                        gradient_lambda(curr_cell,curr_BS) = gradient_lambda(curr_cell,curr_BS) + gradient;
                                        
                                    end
                                end
                            end
                            
                            delta_lambda2 = [gradient_lambda(1) sum(gradient_lambda(2:end))];
                            
                            if norm(lambda) == 0 && delta_lambda2(1) > 0 && delta_lambda2(2) > 0
                               % Cells(1).initial_lambda = lambda;
                                break;
                            elseif abs(max(-delta_lambda2./Total_power_const)) < 2e-2
                               % Cells(1).initial_lambda = lambda;
                                break;
                            else
                                delta_lambda = [delta_lambda2(1) delta_lambda2(2)*ones(1,L_Pico)];
                                lambda = max(lambda-step_size*delta_lambda,0);
                                sub_iter = sub_iter + 1;
                            end
                            
%                             if mod(sub_iter,1000) == 0
%                                 sub_iter
%                             end
                            
                            if mod(sub_iter,1e4) == 0
                                sub_iter
                            end
                            
                            if sub_iter > 1e5
                                return;
                            end

                        end

                    end
                    
                end
                
            end
        end
end


