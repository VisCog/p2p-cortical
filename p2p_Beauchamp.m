% p2p_Beauchamp
%
% to be folded into p2p_c
%
% Updated functions for c2v and v2c with 'squish', plus a fitting function
% 'fitElectrodeGrid' for fitting the Beauchamp electrodes.
%
% Note: I disabled the 'ok' part of c2v_cplx since I'm not using 'isvalidcortex'
% That will have to be reenabled when folding into p2p_c
%
% gmb 3/7/2023

classdef p2p_Beauchamp  % sorry about the lame name
    methods(Static)
        
        function [vx, vy, ok] = c2v_real(c, cx, cy)
            % takes in real x, y numbers on the cortical grid and returns
            % real x, y numbers in visual space
            [z, ok] = p2p_Beauchamp.c2v_cplx(c,cx + sqrt(-1)*cy); % add back p2p_c.
            vx = real(z);
            vy = imag(z);
        end
        
        function w = v2c_cplx(c, z)
            % takes in imaginary numbers in visual space and returns imaginary
            % positions in the cortical grid
            
            lvf = real(z)<0; % find points in the left visual field
            z(lvf) = z(lvf)*exp(-sqrt(-1)*pi);  % rotate lvf points 180 degrees
            w = (c.k*log(z + c.a))-c.shift; % Use the Schwartz!
            w(~lvf) = w(~lvf)*exp(-sqrt(-1)*pi);  % rotate rvf (~lvf) points back 180 degrees
            
            % squish
            w = real(w)+c.squish*sqrt(-1)*imag(w);
            
            
        end
        
        function [cx, cy] = v2c_real(c, vx, vy)
            % takes in real x, y numbers in visual space and returns real
            % x, y positions on the cortical grid
            z = p2p_Beauchamp.v2c_cplx(c,vx + sqrt(-1)*vy);  % add back p2p_c.
            cx = real(z);
            cy = imag(z);
        end
        
        function [w,ok] = c2v_cplx(c, z)
            % takes in imaginary numbers in cortical space and returns
            % imaginary positions in visual space. Not all points in
            % cortical space are valid, so invalid cortical points are
            % returned as NaNs in visual space.  'ok' is a logical where 0
            % is invalid and 1 is valid.
            
            % unsquish
            z = real(z)+ sqrt(-1)*imag(z)/c.squish;
            
            lvf = real(z)>0;  % find points in left visual field (right cortex)
            z(~lvf) = z(~lvf)*exp(-sqrt(-1)*pi); % rotate rvf points 180 degres
            w = exp((z+c.shift)/c.k)-c.a;  % Undo the Schwartz!
            w(lvf) = w(lvf)*exp(-sqrt(-1)*pi);  % rotate lvf points back 180 degrees
            
            % set the invalid points to NaNs
            % ok = p2p_c.isValidCortex(c,real(z),imag(z));
            ok = ones(size(w));  % gmb this needs to be removed
            w(~ok) = NaN;
        end
        
        function [err,predcx,predcy,cx,cy] = fitElectrodeGrid(p,vx,vy)
            % Gives a fit of a projected set of phosphenes onto a 6x4
            % electrode grid by first projecting the phosphenes using the
            % v2c mapping function and then moving and rotating the grid.
            
            
            p.shift = p.k*log(p.a);  % This really should be built in

            % Project phosphenes into cortex (including squishing)
            [cx,cy] = p2p_Beauchamp.v2c_real(p, vx, vy);
            
            % Creat a 6x4 array
            [x0,y0] = meshgrid(-1.5:1.5,-2.5:2.5);
            x0 =x0(:)';
            y0 =y0(:)';
            
            % Expand by dx
            x0 = x0*p.dx;
            y0 = y0*p.dx;
            
            % Rotate by ang and shift by (xc,yc)
            rot = [cos(p.ang) sin(p.ang);-sin(p.ang) cos(p.ang)];
            
            M = rot*[x0;y0] + repmat([p.xc;p.yc],1,24);
            predcx = M(1,:);
            predcy = M(2,:);
            
            % Compare predicted to projected
            err = sum((predcx-cx).^2 + (predcy-cy).^2);
        end
    end
end


