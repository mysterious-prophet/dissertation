%% Calculate Invariant Function
% calculates invariant function 

% inputs: func_name - name of function for invariant calculations
%       : l, m, n - spherical harmonic function parameters
%       : rho - function radius
%       : M, N, P - input image size

function invar_func = calcInvarFunc(func_name, l, m, n, rho, M, N, P)
    if(P == 0)
        [k, mu, ~] = calcFreqCoords(M, N, 0);

        if(strcmp(func_name, 'zerPoly'))
            invar_func = sqrt(2*pi * (2*n + 2)) * (-1i)^n * (-1)^((n-m)/2) * ...
                exp(1i * m * mu) .* besselj(n + 1, k) ./ k;
        elseif(strcmp(func_name, 'sphHarm'))
            invar_func = (-1i)^m * sqrt(2*pi) * exp(1i* m * mu) .* besselj(abs(m), rho * k); 
        end

    elseif(P > 0)
        % calculate omega 3D coordinates 
        [k, mu, psi] = calcFreqCoords(M, N, P);
    
        if(strcmp(func_name, 'zerPoly'))
            invar_func = (-1)^(n) / (1i^(l) * (2*pi)^(l + 1/2)) * ...
                calcSpherHarmFunc(l, m, mu ./(rho * k), psi ./(rho * k)) .* (besselj(2*n + l + 3/2, -rho * k) ./ (-rho * k));
        elseif(strcmp(func_name, 'sphHarm'))
            % spherical besselj(l, rho * k) = sqrt(pi/(2*rho*k)) .* besselj(l + 1/2, rho * k)
            
            % previous version before the new article
            % invar_func = (-1i)^(l) * sqrt(pi) / (4*pi*rho^2 * sqrt(2*rho*k)) .* ...
            %     calcSpherHarmFunc(l, m, mu, psi) .* besselj(l + 1/2, rho * k);

            invar_func = 4*pi * (-1i)^(l) * calcSpherHarmFunc(l, m, mu, psi) .* besselj(l + 1/2, rho * k);

        end
    
        invar_func(isnan(invar_func)) = 0;
        invar_func = conj(invar_func);
    end
end