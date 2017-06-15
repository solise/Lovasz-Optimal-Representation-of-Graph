function beeObject = bee(n,d)

  data_p = zeros(n,d);
  vel_p = zeros(n,d);
  
  beeObject = struct('data',@data_bee,...
                      'set',@set_element,...
                      'setVel',@set_vel,...
                      'dataVel',@vel_bee,...
                      'setData',@set_data);

 function data = data_bee()
    %# Displays the data_p in the list
   data = data_p;
 end

function set_element(vector,index)
    %# Overwrites an element at an index in the list
    data_p(index,:) = vector;
end

function set_data(data)
    %# Overwrites an element at an index in the list
    data_p = data;
end

 function Vel = vel_bee()
    %# Displays the data_p in the list
   Vel = vel_p;
 end

function set_vel(data)
    %# Overwrites an element at an index in the list
    vel_p = data;
end


end