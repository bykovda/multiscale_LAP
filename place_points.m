function [XS2, YS2, test_fun, I_fun, regions] = place_points(model, N, opts)
% place N^2 points 
size1 = opts(1);
I_fun=[];
regions={};

% ä.á.  place_points(model, N, size1[, parameters])
% e.g.  place_points('text', 100, W, 'AB')

N = round(N);

scale1 = 2;
scale2 = 1/2;
done = false;


if ~ischar(model) %raster
	W1 = size(model,2);
	model2 = false(size(model)+2);
	model2(2:end-1, 2:end-1) = model;
	
	window = @(X, maxind) X .* (X >= 1 & X <= maxind) + 1 * (X < 1) + maxind * (X > maxind);
	window_x = @(X)window(round(X), size(model2,2));
	window_y = @(Y)window(round(Y), size(model2,1));
	testcond1 = @(X,Y)model2(sub2ind(size(model2),Y,X));
	testcond = @(X,Y)testcond1(window_x(W1*X/2*2+size(model2,2)/2), window_y(W1*Y/2*2+size(model2,1)/2));
	
	W = 1+0.65;
	S = sum(model(:))/numel(model);
else
	switch lower(model)
		case {'einstein','einstein0.01', 'gradient'}
            W2 = size1;
            if N > 0
                dr = 'D:\_Projects\2018-04-06 - Lloyd spacing\mfkasim91-stippling-lloyds-a326873';
                addpath(dr);
                addpath([dr '\lib']);
                switch lower(model)
                    case 'einstein'
                        st_filename = ['D:\_Projects\2018-04-06 - Lloyd spacing\Einstein-' num2str(N) '.mat'];
                        filepath = [dr '/..'];
                        filename = 'Einstein.png';
                        bg=0;
                    case 'einstein0.01'
                        st_filename = ['D:\_Projects\2018-04-06 - Lloyd spacing\Einstein0.01-' num2str(N) '.mat'];
                        filepath = [dr '/..'];
                        filename = 'Einstein.png';
                        bg=0.01;
                    case 'gradient'
                        filepath = 'D:\_Projects\2018-02-05 - LAP Mirror Lambertian DD\';
                        st_filename = [filepath, '\Gradient-' num2str(N) '.mat'];
                        filename = 'Gradient.png';
                        bg=0;
                end
                
                colorImg = imread([filepath '\' filename]);
                img = double(rgb2gray(colorImg));
                W = size(img, 1);
                H = size(img, 2);
                Cx = W/2;
                Cy = W/2;

                if exist(st_filename, 'file') == 2
                    load(st_filename, 'Px', 'Py');
                else
                    [Px, Py, Ap] = stipple_image(filepath, filename, N^2, 250, [],[],[],bg); %#ok<ASGLU>
                    %figure;hist(Ap,100)
                    save(st_filename, 'Px', 'Py', 'Ap');
                end
                %Cx = (max(Px) +  min(Px))/2;
                %Cy = (max(Py) +  min(Py))/2;
        
                %W = max(Px) - min(Px);
                %H = max(Py) - min(Py);
                XS2 = W2*(Px-Cx)/max(W, H);
                YS2 = W2*(Cy-Py)/max(W, H);
            else
                XS2 = [];
                YS2 = [];
                W = 560;
                H = 560;
            end
			test_fun = @(x, y) abs(x) <= W2/2 * W/max(W, H) & abs(y) <= W2/2 * W/max(W, H);
			
			done = true;
		case 'chessboard'
			test1d = @(x) ((x >= -4) & (x <= -3)) | ((x >= -2) & (x <= -1)) | ((x >= 0) & (x <= 1)) | ((x >= 2) & (x <= 3));
			
			testcond1 = @(X,Y) (xor(test1d(X), test1d(Y)) | abs(X) >= 4 | abs(Y) >= 4) & abs(X) <= 4.05 & abs(Y) <= 4.05;
			testcond = @(X,Y)testcond1(8.1*X, 8.1*Y);

			W = (8.1*sqrt(2))/8.1;
			S = ((4.05*2)^2-(4*2)^2 + 32) / (8.1^2);
            
            center_xy = -3.5:1:3.5;
            ind_xy = [1 0 1 0 1 0 1 0];
            [IX, IY] = ndgrid(ind_xy, ind_xy);
            I = xor(IX, IY);
            [CX, CY] = ndgrid(center_xy,center_xy);
            
            for i = -6:2:6
                regions{end+1} = @(x,y) (8.1*y >= 8.1*x + i);
            end
            for i = -7:2:7
                regions{end+1} = @(x,y) (8.1*y >= -8.1*x + i);
            end
            
            %{
            for i = 1:numel(I)
                if I(i)
                    regions{end+1} = @(x,y) (8.1*x-CX(i)).^2+(8.1*y-CY(i)).^2 <= 1/2;
                end
            end
            regions{end+1} = @(x,y) 8.1*x >= 4  & 8.1*x <=4.05 & abs(8.1*y) <= 4.05;
            regions{end+1} = @(x,y) 8.1*x <= -4 & 8.1*x >=-4.05 & abs(8.1*y) <= 4.05;
            regions{end+1} = @(x,y) 8.1*y >= 4  & 8.1*y <=4.05 & abs(8.1*x) <= 4.05;
            regions{end+1} = @(x,y) 8.1*y <= -4 & 8.1*y >=-4.05 & abs(8.1*x) <= 4.05;
            %}
            
            %figure;
            %plot(CX(I), CY(I), 'x')
            
			
		case {'normal', 'gauss'}
			sigma = 0.5;
			r = linspace(0,size1,N.^2);
			gr = (sqrt(3)+1)/2;
			phi = linspace(0,1,N.^2)*(N^2/gr^2)*2*pi;
			%r = rand(1, N^2)*size;
			%phi = rand(1, N^2)*2*pi;
			u = sigma*sqrt(-2*log(1-r*(1-exp(-size1^2/(2*sigma^2)))));
			XS2 = u.*cos(phi);
			YS2 = u.*sin(phi);
			XS2 = XS2(:);
			YS2 = YS2(:);

			test_fun = @(X,Y)X.^2+Y.^2<=size1^2;	
			done = true;
			
		case 'square'
			%testcond = @(X,Y)abs(X)<=1/2 && abs(Y) <= 1/2;		
			%W = 1.25;
			%S = 1;		
			xs = linspace(-size1/2, size1/2, N)*(N-1)/N;%2019-02-12  added (n-1)/n
			[XS2, YS2] = ndgrid(xs, xs);
			XS2 = XS2(:);
			YS2 = YS2(:);
			test_fun = @(X,Y)(abs(X)<= size1/2) & (abs(Y)<= size1/2);

			done = true;
            
		case 'squarering'            
			r2 = opts(2)/opts(1); % w_inner / w_outer
            
			testcond = @(X,Y)abs(X)<=1/2 & abs(Y) <= 1/2 & ~((abs(X) <= r2/2) & (abs(Y) <= r2/2));
			W = 2%sqrt(2)*1.1;
            S = 1 - r2^2;
            
		case 'ring'
			r2 = opts(2)/opts(1); % r_inner / r_outer		

			testcond = @(X,Y)X.^2+Y.^2<=1 & X.^2+Y.^2>=r2^2;
			W = 3;
			S = pi * (1-r2^2);
			            
		case 'squarehex'
			testcond = @(X,Y)abs(X)<=1/2 & abs(Y) <= 1/2;		
			W = 2;
			S = 1;
			
		case 'rectangle'
			testcond = @(X,Y)abs(X)<=1/2 & abs(Y*2) <= 1/2;		
			W = 1.25;
			S = 1/2;

			%xs = linspace(-size1/2, size1/2, N)*(N-2)/(N-1);
			%[XS2, YS2] = ndgrid(xs, xs);
			%XS2 = XS2(:);
			%YS2 = YS2(:)/2;
			%test_fun = @(X,Y)(abs(X)<= size1/2) & (abs(Y*2)<= size1/2);
            
			%done = true;
		case 'sinrectangle'
			%testcond = @(X,Y)abs(X)<=1/2 & abs(Y*2) <= 1/2;		
			%W = 1.25;
			%S = 1/2;

			xs = linspace(0, 1, N);
            xs2 = xs;
            
            for i = 1:100
                %xs2 = -4*sin(5*pi*xs2/2).^2 / (15 * pi) - 5 .* xs./ 3;
                xs2 = ((4 + 15 * pi).*xs - 2 + 2 .* cos(5 * pi .* xs2)) / (15 * pi);                
            end            
			[XS2, YS2] = ndgrid(xs2, xs);
			XS2 = XS2(:)*size1 - size1/2;
			YS2 = (YS2(:)-1/2)*size1/2;
            
            test_fun = @(X,Y)(abs(X)<= size1/2) & (abs(Y*2)<= size1/2);
			done = true;
			
        case 'sinsquare'
			%testcond = @(X,Y)abs(X)<=1/2 & abs(Y*2) <= 1/2;		
			%W = 1.25;
			%S = 1/2;

			xs = linspace(0, 1, N);
            xs2 = xs;
            
            for i = 1:100
                %xs2 = -4*sin(5*pi*xs2/2).^2 / (15 * pi) - 5 .* xs./ 3;
                %xs2 = ((4 + 15 * pi).*xs - 2 + 2 .* cos(5 * pi .* xs2)) / (15 * pi);                
                %xs2 = sqrt(pi)/2 * exp(xs2.^2) .* ( erf(xs2) - (-1 + 2 .*xs) .* erf(6));
                xs2 = pi*(4.* xs - 2) + sin(xs2);
            end
			[XS2, YS2] = ndgrid(xs2, xs2);
			XS2 = XS2(:)*size1 - size1/2;
			YS2 = YS2(:)*size1 - size1/2;
            
            test_fun = @(X,Y)(abs(X)<= size1/2) & (abs(Y*2)<= size1/2);
            I_fun_1d = @(x)1+bg+(1-bg).*sin(x-pi/2);
            I_fun = @(x,y)I_fun_1d(x).*I_fun_1d(y).*test_fun(x,y);
			done = true;
			

        case 'sinsquare2'
			%testcond = @(X,Y)abs(X)<=1/2 & abs(Y*2) <= 1/2;		
			%W = 1.25;
			%S = 1/2;
            
            
            [xs1, ys1] = place_points('squarehex', N, 1);            
            xs1 = xs1 + 1/2;
            ys1 = ys1 + 1/2;

            xs2 = xs1;
            ys2 = ys1;
            
            for i = 1:500
                %xs2 = -4*sin(5*pi*xs2/2).^2 / (15 * pi) - 5 .* xs./ 3;
                %xs2 = ((4 + 15 * pi).*xs - 2 + 2 .* cos(5 * pi .* xs2)) / (15 * pi);                
                %xs2 = sqrt(pi)/2 * exp(xs2.^2) .* ( erf(xs2) - (-1 + 2 .*xs) .* erf(6));
                %xs2 = 0.9*xs2+0.1*(pi*(4.* xs1 - 2) + sin(xs2));
                %ys2 = 0.9*ys2+0.1*(pi*(4.* ys1 - 2) + sin(ys2));
                
                bg = 0.17321;
                
                %bg = 0.70711;
                
                xs2 = 0.9*xs2+0.1*(pi*(4.* xs1 - 2) + (1-bg)./ (1+bg) .* sin(xs2));
                ys2 = 0.9*ys2+0.1*(pi*(4.* ys1 - 2) + (1-bg)./ (1+bg) .* sin(ys2));
            end
            
            XS2 = xs2(:) * size1 / (4 * pi);
			YS2 = ys2(:) * size1 / (4 * pi);
            
            test_fun = @(X,Y)(abs(X)<= size1/2) & (abs(Y)<= size1/2);
            I_fun_1d = @(x)1+bg+(1-bg).*sin(x*4*pi/size1-pi/2);
            I_fun = @(x,y)I_fun_1d(x).*I_fun_1d(y).*test_fun(x,y);
			done = true;
            
		case 'lower-rectangle'
			%testcond = @(X,Y)abs(X)<=1/2 & abs(Y*2) <= 1/2;		
			%W = 1.25;
			%S = 1/2;

			xs = linspace(-size1/2, size1/2, N);
			[XS2, YS2] = ndgrid(xs, xs);
			XS2 = XS2(:);
			YS2 = YS2(:)/2 - size1/4;
			test_fun = @(X,Y)(abs(X)<= size1/2) & (abs(Y+size1/4)<= size1/4);
			done = true;
			
		case {'sunflower-circle', 'sunflower', 'circle-sunflower'}
			[XS2, YS2] = sunflower_points_on_circle(N^2);
			maxr = 1.000001*sqrt(max(max(XS2.^2+YS2.^2)));
			XS2 = XS2(:) * size1 / maxr;
			YS2 = YS2(:) * size1 / maxr;		
			test_fun = @(X,Y)X.^2+Y.^2<=size1^2;
			done = true;
			
		case 'circle'
			testcond = @(X,Y)X.^2+Y.^2<=1;		

			W = 2.25;
			S = pi;
			
		case 'circle-cutted'
			testcond = @(X,Y)X.^2+Y.^2<=1 & X - 0.25*Y.^2 <= sind(50);

			W = 2.25;
			S = pi;
			
		case 'whirl'
			 [XS2, YS2, test_fun] = place_points('circle', N, size1);
			 phi = atan2(YS2, XS2)+pi;
			 rho = sqrt(YS2.^2 + XS2.^2);
			 %phi = 0.8*sqrt(phi*(2*pi));
			 bg=0.5;
			 phi = (sqrt(bg^2*(2*pi-phi) + phi)*sqrt(2*pi)-bg*2*pi)/(1-bg);
			 XS2 = rho .* cos(phi);
			 YS2 = rho .* sin(phi);
			 return;

		case 'star'
			testa = @(x,y,phi)x*sin(phi) + y*cos(phi) < 1;
			testcond = @(x,y)testa(x, y, 0) + testa(x, y, 2*pi/5) + testa(x, y, 2*2*pi/5) + testa(x, y, 3*2*pi/5) + testa(x, y, 4*2*pi/5)>=4;

			S = 5 * sqrt(10 - 2*sqrt(5));
			W = 8;		


		case 'batman'
			UnitStep = @(x)(x>0)*1.0;
			w = @(x)3*sqrt(1 - (x/7).^2);
			l = @(x)(6/7)*sqrt(10) + (3 + x)/2 - (3/7)*sqrt(10)*sqrt(4 - (x + 1).^2);
			h = @(x)(1/2)*(3*(abs(x - 1/2) + abs(x + 1/2) + 6) - 11*(abs(x - 3/4) + abs(x + 3/4)));
			r = @(x)(6/7)*sqrt(10) + (3 - x)/2 - (3/7)*sqrt(10)*sqrt(4 - (x - 1).^2);
			u = @(x)w(x) + (l(x) - w(x)).*UnitStep(x + 3) + (h(x) - l(x)).*UnitStep(x + 1) + (r(x) - h(x)).* UnitStep(x - 1) + (w(x) - r(x)).*UnitStep(x - 3);
			d = @(x)(1/2)*(3*sqrt(1 - (x/7).^2) + sqrt(1 - (abs(abs(x) - 2) - 1).^2) + abs(x/2) - ((3*sqrt(33) - 7)/112)*x.^2 - 3).*((x + 4)./abs(x + 4) - (x - 4)./abs(x - 4)) - 3*sqrt(1 - (x/7).^2);
			testcond = @(x,y)7*y<u(7*x) & 7*y>d(7*x);

			W = 1.25;
			S = 0.2468;

		case 'triangle'		
			testcond = @(X,Y)Y >=- 1/(2*sqrt(3)) & Y <= 1/sqrt(3) - X*sqrt(3) & Y <= 1/sqrt(3) + X*sqrt(3);		

			W = 1.25;
			S = sqrt(3)/4;

		case 'triangle+square'		
			tri_cond = @(X,Y)Y >=- 1/(2*sqrt(3)) & Y <= 1/sqrt(3) - X*sqrt(3) & Y <= 1/sqrt(3) + X*sqrt(3);		
			sq_cond = @(X,Y)abs(X)<=1/2 & abs(Y) <= 1/2;		

			testcond = @(X,Y)tri_cond(X+0.75, Y+1/(4*sqrt(3))) | sq_cond(1.1*(X-0.75),1.1*Y);
			W = 3;
			S = 1/1.1^2+sqrt(3)/4;

		case 'cross'
			fsq= @(X,Y,dx,dy) (X-dx)<=1/2 & (X-dx)>=-1/2 & (Y-dy)<=1/2 & (Y-dy)>=-1/2 ;
            %testcond = @(X,Y)fsq(X,Y,0,0) | fsq(X,Y,0.9,0) | fsq(X,Y,-0.9,0) | fsq(X,Y,0,0.9) | fsq(X,Y,0,-0.9);
%			W = 3;
%			S = 4.6;

            testcond = @(X,Y)fsq(3*X,3*Y,0,0) | fsq(3*X,3*Y,0.999999,0) | fsq(3*X,3*Y,-0.99999,0) | fsq(3*X,3*Y,0,0.99999) | fsq(3*X,3*Y,0,-0.99999);
			W = 3/3*1.25;
			S = 4.6/9;
        otherwise
            error(['place_points:: Unknown figure "' model '"'])
	end
end

if ~done
	N2 = round( N * sqrt(W^2 / S) );

	mesh = 'hexagonal rotated';
	%mesh = 'square';

	switch mesh
		case {'square', 'square rotated'}
			%dx = (rand()-1/2)*W/(N2-1);
			%dy = (rand()-1/2)*W/(N2-1);
			dx=0.00012;
			dy=0.052156;
			xst = linspace(-W/2, W/2, N2) + dx;
			yst = linspace(-W/2, W/2, N2) + dy;
			[XSt, YSt] = meshgrid(xst * size1, yst * size1);
		case {'hexagonal', 'hex', 'hexagonal rotated', 'hex rotated'}
			q = sqrt(2)/sqrt(sqrt(3));
			a = q * W/(N2-1); %triangle side
			h = a * sqrt(3)/2;

	%		dx = (rand()-1/2)*a;
	%		dy = (rand()-1/2)*h;
			dx=0.000012;
			dy=0.0512156;


			xst = (-W/2):a:(W/2);
			yst = (-W/2):h:(W/2);
			
			
			xst = xst - (max(xst) + a/2 + min(xst))/2 + dx;
			yst = yst - (max(yst) + min(yst))/2 + dy;
			[XSt, YSt] = meshgrid(xst * size1, yst * size1);
			XSt(1:2:end,:) = XSt(1:2:end,:) + size1 * a/2;
	end

	switch mesh
		case {'square rotated', 'hexagonal rotated', 'hex rotated'}
			%phi = 2*pi*rand();
			phi = pi/2+pi/exp(1);
			XSt1 = XSt * cos(phi) + YSt * sin(phi);
			YSt = YSt * cos(phi) - XSt * sin(phi);
			XSt = XSt1;
	end

	test_fun = @(x,y)testcond(x/size1,y/size1);
    regions2 = regions;
    for i = 1:numel(regions2)
        regions{i} = @(x,y)regions2{i}(x/size1,y/size1);
    end
	
	if N == 0 % only test function is needed
		XS2 = [];
		YS2 = [];
    else

        fn = @(s)sum(sum(testcond(XSt*s/size1,YSt*s/size1)));
        assert(fn(scale1)< N^2 && fn(scale2) > N^2)

        %{
        x1 = XSt*scale1/size1;
        y1 = YSt*scale1/size1;
        ind = testcond(x1, y1);
        figure;
        subplot 211
        hold all
        plot(x1(ind),y1(ind), '.k');
        plot(x1(~ind),y1(~ind), '.r');
        daspect([1 1 1])

        x1 = XSt*scale2/size1;
        y1 = YSt*scale2/size1;
        ind = testcond(x1, y1);
        subplot 212
        hold all
        plot(x1(ind),y1(ind), '.k');
        plot(x1(~ind),y1(~ind), '.r');
        daspect([1 1 1])
        %}

        f_r = 0;
        while f_r ~= N^2	
            scale = (scale1 + scale2)/2;
            f_r = fn(scale);
            if  f_r > N^2
                scale2 = scale;		
            else
                scale1 = scale;
            end
            if abs(scale2- scale1)<1e-8 && f_r ~= N^2 
                error('Cannot distribute points')
            end
        end

        ind = testcond(XSt*scale/size1,YSt*scale/size1);

        XS2 = XSt(ind(:))*scale;
        YS2 = YSt(ind(:))*scale;

    end
end

if N > 0
    [~, ind] = sort(YS2+0.0001*XS2, 'descend');
    XS2 = XS2(ind);
    YS2 = YS2(ind);
end

if isempty(I_fun)
    I_fun = @(x,y)test_fun(x,y)*1.0;
end
