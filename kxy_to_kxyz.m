function [place_fun, test_fun] = kxy_to_kxyz(place_fun_2d, kz_sign)
    if nargin < 2
        kz_sign = 1; % shining up by default
    end
	function [kx, ky, kz] = ff(N)
		[kx, ky] = place_fun_2d(N);
		kz = kz_sign * sqrt(1 - kx.^2 - ky.^2);
	end
	place_fun = @ff;
    if nargout >= 2
        [~, ~, test_fun] = place_fun_2d(0);
    end
end