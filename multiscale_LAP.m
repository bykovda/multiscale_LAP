function [X1_final, Y1_final, Z1_final,  X2_final, Y2_final, Z2_final] = ...
	multiscale_LAP(placefun1, placefun2, cost_fun, sign, n_max, n_neighbors, inverse_continuous)

%TODO: Проверять существование паросочетания до(!) вычисления функции стоимости
%TODO: у place_points первым аргументом нужно сделать n

%% Input arguments parsing

cf = cost_fun(1);

if cf == 5 %compatibility 2020-05-06
    warning('Use cost_fun 11 instead of 5')
    cf = 11;
end
if cf == 55 %compatibility 2020-05-06
    warning('Use cost_fun 12 instead of 55')
    cf = 12;
end


if numel(cost_fun) > 1
	param = cost_fun(2); % 1/n or gamma or f0
else
	param = 1;
end
param_squared = param^2;

if nargin < 6 || isempty(n_neighbors)
	n_neighbors = 15;
end

if nargin < 7 || isempty(inverse_continuous)
	inverse_continuous = true;
end
if ~inverse_continuous
	if cf == 11  % changing cost fun if non-symemtriс
		cost_fun(1) = 12;
	end
	[X2_final, Y2_final, Z2_final,  X1_final, Y1_final, Z1_final] = ...
		multiscale_LAP(placefun2, placefun1, cost_fun, sign, n_max, n_neighbors, true);
	return
end

%% Warm-up

t = tic;

switch cf
	case 10		% -(x1 x2 + y1 y2)
		dim1=2;
		dim2=2;
	case 0		% -(x1 x2 + y1 y2 + z1 z2)
		dim1=3;
		dim2=3;
	case 1		% -log(1-param(x1 x2 + y1 y2 + z1 z2))
		dim1=3;
		dim2=3;
	case 2	    % sqrt(param^2 + (x1-x2)^2 + (y1-y2)^2)
		dim1=2;
		dim2=2;
	case 3		% -sqrt(param^2 - (x1-x2)^2 - (y1-y2)^2)
		dim1=2;
		dim2=2;
	case 4		% sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)
		dim1=3;
		dim2=3;
	case 11		% ( x1 x2 + y1 y2)/(param-z2) 
		dim1=2;
		dim2=3;
	case 12		% ( x1 x2 + y1 y2)/(param-z1) 
		dim1=3;
		dim2=2;
end


if numel(n_max) > 2
	ns = n_max;
else
	if numel(n_max) == 2
		n_min = n_max(1);
		n_max = n_max(2);
	else
		n_min = 20;
    end
    if n_min == n_max %% without multiscale
        ns = n_max;
    else
        decrease_rate = sqrt(2);
        nsteps = ceil((log(n_max) - log(n_min))/log(decrease_rate));
        ns = linspace(log2(n_max), log2(n_min), nsteps+1);
        ns = round(2.^ns);
        ns(2:end-1) = round(ns(2:end-1)/10)*10;
        ns = ns([end:-1:1 1 1 ]);%1
    end
end 
nsteps = numel(ns);

if iscell(placefun1)
	pf1 = placefun1;
	placefun1 = @(N)place_points(pf1{1}, N, pf1{2:end});
	assert(dim1==2, 'place_points used in placefun1 returns 2D points instead of 3D as required by the Cost fun');
end

if iscell(placefun2)
	pf2 = placefun2;
	placefun2 = @(N)place_points(pf2{1}, N, pf2{2:end});
	assert(dim2==2, 'place_points used in placefun2 returns 2D points instead of 3D as required by the Cost fun');
end

Z1_final = [];
Z2_final = [];

n_neighbors_0 = n_neighbors;

ni = 0;

%% Main cycle

while ni < numel(ns)

%% Placing points
	ni = ni + 1;
	N = ns(ni);
	disp([num2str(ni) '/' num2str(numel(ns)), ': N = ' num2str(N) '; n_neighbors = ' num2str(n_neighbors)]);
	
	if dim1 == 2
		[X1, Y1] = placefun1(N);
	else
		[X1, Y1, Z1] = placefun1(N);
	end
	
	if dim2 == 2
		[X2, Y2] = placefun2(N);
	else
		[X2, Y2, Z2] = placefun2(N);
	end
	
	
%% Calculating cost function
	
	if ni == 1 
		A = zeros(N^2, N^2);
		is = 1:N^2;
		
		x1s = X1(:);
		y1s = Y1(:);
		if dim1 == 3
			z1s = Z1(:);
		end
		dx = max(X1(:)) - min(X1(:)); % why only at ni==1?
		dy = max(Y1(:)) - min(Y1(:));
	else
		A = spalloc(N^2, N^2, N^2*n_neighbors);
	end
	
	R = zeros(N);
	X1_C = zeros(N);
	Y1_C = zeros(N);
	Z1_C = zeros(N);
	
	for j = 1:N^2
		if ni > 1
			if dim2 == 2
				dist2_squared = (X2(j) - X2_final).^2 + (Y2(j) - Y2_final).^2;
			else
				dist2_squared = (X2(j) - X2_final).^2 + (Y2(j) - Y2_final).^2 + (Z2(j) - Z2_final).^2;
			end
			[~, last_j] = min(dist2_squared);
			
			if dim1 == 2
				dist1_squared = (X1(:) - X1_final(last_j)).^2 + (Y1(:) - Y1_final(last_j)).^2;
			else
				dist1_squared = (X1(:) - X1_final(last_j)).^2 + (Y1(:) - Y1_final(last_j)).^2 + (Z1(:) - Z1_final(last_j)).^2;
			end
			X1_C(j) = X1_final(last_j);
			Y1_C(j) = Y1_final(last_j);
			if dim1 == 3
				Z1_C(j) = Z1_final(last_j);
			end
			
			
			curR = 2*sqrt(n_neighbors*(dx*dy /  N^2) / pi);
			
            count = 1;            
			repeat = true;
            upper_multiple = 5;
			while repeat
				ind2 = find(dist1_squared <= curR^2);
				n_found = numel(ind2);
				if n_found < n_neighbors
					curR = curR * 1.5;  %2
				elseif n_found > upper_multiple * n_neighbors
					curR = curR / 1.5;  %2
				else
					repeat = false;
                end
                count = count+1;
                if count > 1000 % in case we are stuck
                    upper_multiple = upper_multiple*1.1;
                end
                
			end
			R(j) = curR;
			
			[~, sorted_i] = sort(dist1_squared(ind2));
			is = ind2(sorted_i(1:n_neighbors));
			is = sort(is);
	
%			[~, sorted_i] = sort(dist1_squared);
%			is2 = sorted_i(1:n_neighbors);
%			assert(all(is == is2));

			x1s = X1(is);
			y1s = Y1(is);
			if dim1 == 3
				z1s = Z1(is);
			end
		end
		
		switch cf
			case 10
				new_row = -(x1s * X2(j) + y1s * Y2(j));
			case 0
				new_row = -(x1s * X2(j) + y1s * Y2(j) + z1s * Z2(j));
			case 1
				new_row = -log(1 - param*(x1s * X2(j) + y1s * Y2(j) + z1s * Z2(j)));
				new_row(imag(new_row)~=0) = inf;
			case 2
				new_row = sqrt(param_squared + (x1s - X2(j)).^2 + (y1s - Y2(j)).^2);
			case 3
				new_row = -sqrt(param_squared - (x1s - X2(j)).^2 - (y1s - Y2(j)).^2);
				new_row(imag(new_row)~=0) = inf;
			case 4
				new_row = sqrt( (x1s - X2(j)).^2 + (y1s - Y2(j)).^2 + (z1s - Z2(j)).^2);
			case 11
				new_row = (x1s * X2(j) + y1s * Y2(j)) ./ (param - Z2(j));
			case 12
				new_row = (x1s * X2(j) + y1s * Y2(j)) ./ (param - z1s);
		end
		A(is, j) = new_row*sign + 0.000121315649795;  % to remove zero elements!
	end
	A = A';
	
%% Norming cost function
	
	mx = max(A(:));
	mn = min(A(:));
	
	if ni == 1		
		A = 1 + (10000000*(A-mn)/(mx-mn));
	else
		[i, j, s] = find(A);
		clear A
		s = 1 + (10000000*(s-mn)/(mx-mn));
		A = sparse(i, j, s);
		clear i j s
	end
	
%% Solving sparse assignment problem
	
%	tic
	failed = false;
	if ni == 1 
		b = assignmentProblemAuctionAlgorithm(A)';
	else
		assert(nnz(A) == N^2*n_neighbors);
		%b = sparseAssignmentProblemAuctionAlgorithm(A')';
		b_i = sparseAssignmentProblemAuctionAlgorithm(A)';
		if b_i(1) == -1
			failed = true;
		else
			b = perminv(b_i);
            disp(['Average cost = ' num2str(sum(full(A(sub2ind(size(A),1:size(A,1), b_i))))/size(A,1))])
		end
	end
	clear A
%	toc
	
	if failed
		n_neighbors = n_neighbors + 25;
		ns = ns([1:ni, ni:end]);
		assert(numel(ns) < 5*nsteps);
	else
		n_neighbors = max(n_neighbors_0, n_neighbors - 25);
		
		X1_final = X1;
		Y1_final = Y1;
		X2_final = X2(b);
		Y2_final = Y2(b);
		
		if dim1 == 3
			Z1_final = Z1;
		end
		if dim2 == 3
			Z2_final = Z2(b);
		end
		
%% Checking convergence
		
		if ni > 1
			if dim1 == 2
				Delta1 = (X1(b_i) - X1_C(:)).^2 + (Y1(b_i) - Y1_C(:)).^2;
			else	
				Delta1 = (X1(b_i) - X1_C(:)).^2 + (Y1(b_i) - Y1_C(:)).^2 + (Z1(b_i) - Z1_C(:)).^2;
			end
			max_delta_R = sqrt(max(Delta1./(R(:).^2)));
	%		figure(41516)
	%		subplot(5, 3, min(15,ni-1))
	%		hist(sqrt(Delta1)./R(:), 100);
	%		title(num2str([ni N n_neighbors]));
	%		drawnow
	%		disp(['delta / R = ' num2str(max_delta_R*100) '%'])
			if (max_delta_R==0) && (ns(ni-1) == ns(ni))
                disp('Done');
				return
			end
		end
    end
end
disp(['multiscale_LAP Done in ' hms(toc(t))]);