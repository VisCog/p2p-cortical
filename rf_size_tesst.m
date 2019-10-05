f_data = [
0.5002811212722751, 1.1996420006655129
2.0590698901905906, 1.2033711603997688
3.0183245172172466, 1.2056660279285438
4.815779509127838, 1.0903489346077535
6.615529368567204, 1.214271781161436
11.539741368429508, 2.063372766807035
19.10534589390828, 3.277644547968471
28.234328923362895, 4.974125368613102
40.24222326766188, 6.79711076178128
49.36891142958773, 8.254351642551438];

 % Freeman and Simoncelli, 2011
                    cf.slope = 0.1644; % receptive field diameter
                    cf.min = 1.21;
                    cf.intercept = 0;
%                     
%                 elseif strcmp(c.rfmodel, 'smirnakis')
%                     c.slope =  0.08; % in terms of sigma
%                     c.intercept = 0.16;
%                     c.min = 0;
%    end
                
ecc = 0:50;   
plot(f_data(:, 1), f_data(:, 2), 'ro'); hold on
    cf.rfsize = max((cf.slope .* [0:50]), cf.min)+cf.intercept;
    plot(ecc, cf.rfsize, 'r-');
    
    % 'smirnakis')
                  cs.slope =  0.08; % in terms of sigma
                     cs.intercept = 0.16;
                     cs.min = 0; 
      cs.rfsize = max((cs.slope .* [0:50]), cs.min)+cs.intercept;  
       plot(ecc, cs.rfsize*2, 'b-'); 