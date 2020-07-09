function [sim_state_cell, sim_emission_cell, sim_time_cell] = stochastic_sim_fun(Q,ss_vec,emission_vec,n_sim,t_sim)

  % state options
  state_option_vec = 1:size(Q,1);
  
  % initialize cell arrays
  sim_state_cell = cell(1,n_sim);
  sim_emission_cell = cell(1,n_sim);
  sim_time_cell = cell(1,n_sim);
  
  % iterate
  parfor n = 1:n_sim
      state_val_vec = [randsample(state_option_vec,1,true,ss_vec)];
      jump_time_vec = [0];
      t_curr = 0;
      while t_curr < t_sim            
          state_curr = state_val_vec(end);   % current state  
          next_jump = exprnd(-1/Q(state_curr,state_curr)); % time til next jump              
          t_next = t_curr + next_jump;

          % randomly choose next state
          weight_vec = Q(state_curr,:);
          weight_vec(state_curr) = 0;
          weight_vec = weight_vec / sum(weight_vec);
          state_next = randsample(state_option_vec,1,true,weight_vec);
          
          % if within alotted time, add state to state vector
          if t_next < t_sim
              state_val_vec = [state_val_vec state_next];
              jump_time_vec = [jump_time_vec t_next];
          end
          t_curr = t_next;
      end
      sim_state_cell{n} = state_val_vec;    
      sim_state_cell{n} = emission_vec(state_val_vec);    
      sim_time_cell{n} = jump_time_vec;
  end  