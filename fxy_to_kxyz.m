function [place_fun, test_fun] = fxy_to_kxyz(place_fun_2d, f0)
	function [kx, ky, kz] = ff(N)
		[x, y] = place_fun_2d(N);
		r = sqrt(x.^2+y.^2+f0.^2);
		kx = f0 ./ r;
		ky = x ./ r;
		kz = y ./ r;
	end
	place_fun = @ff;
    if nargout >= 2
        [~, ~, test_fun] = place_fun_2d(0);
    end
end