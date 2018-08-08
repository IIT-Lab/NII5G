function [Scheduled_User UnScheduled_User]= Round_Robin(User_pool,Num_Scheduled_User)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This function performs round robin scheduling to select
%%%% Num_Scheduled_User's from the User_pool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Scheduled_User = zeros(1,Num_Scheduled_User);
UnScheduled_User = User_pool;
if length(User_pool) < Num_Scheduled_User
    fprintf('$$$ There are not enough candidate users to select from $$$ \n');
    return;
else
    for k = 1:Num_Scheduled_User
        user_index = randi(length(UnScheduled_User)); 
        %uniformly choose a integer number from 1:length(UnScheduled_User)
        Scheduled_User(k) = UnScheduled_User(user_index);
        
        UnScheduled_User(user_index) = [];
    end
    Scheduled_User = sort(Scheduled_User);
    
    end

end
