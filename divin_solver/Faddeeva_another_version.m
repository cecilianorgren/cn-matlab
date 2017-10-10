function f = faddeeva(zmat)
% # function f = faddeeva(z)
% # Calculates Faddeeva function. 
% #
% # Copyright (C) 2004 Victor Munoz 
% #
% # Based on code in Matpack (version 1.7.3).
% #
% # Copyright (C) 1991-2002 by Berndt M. Gammel
% 
% #Version 1.0

DBL_MAX = 1.797693134862316e+308;
DBL_EPSILON= 2.220446049250313e-16;
M_2_SQRTPI= 2/sqrt(pi);

%     # The maximum value of rmaxreal equals the root of the largest number 
%     # rmax which can still be implemented on the computer in double precision
%     # floating-point arithmetic
rmaxreal = sqrt(DBL_MAX);

%     # rmaxexp  = ln(rmax) - ln(2)
rmaxexp  = log(DBL_MAX) - log(2.0);

%    # the largest possible argument of a double precision goniometric function
rmaxgoni = 1.0 / DBL_EPSILON;
h2 = 0;
u2 = 0;
v2 = 0;
qlambda = 0;

f=zeros(size(zmat));

for fila=1:size(zmat,1),
  for columna=1:size(zmat,2),

    z=zmat(fila,columna);

    xi = real(z);
    yi = imag(z);
    xabs = abs(xi);
    yabs = abs(yi);
    x = xabs / 6.3;
    y = yabs / 4.4;

%    # the following statement protects qrho = (x^2 + y^2) against overflow
    if ((xabs > rmaxreal) | (yabs > rmaxreal)),
      f(fila,columna)=NaN;
      message('Faddeeva, absolute value of argument so large w(z) overflows')
    end


    qrho = x .* x + y .* y;
    xabsq = xabs .* xabs;
    xquad = xabsq - yabs .* yabs;
    yquad = xabs .* 2 .* yabs;
    
    a = qrho < 0.085264;

    if (a),
	
%	# If (qrho < 0.085264) then the Faddeeva-function is evaluated 
%        # using a power-series (Abramowitz/Stegun, equation (7.1.5), p.297). 
%        # n is the minimum number of terms needed to obtain the required 
%        # accuracy.
	
      qrho = (1 - y * 0.85) * sqrt(qrho);
      n = round(qrho * 72 + 6);
%# C++      j = (n << 1) + 1; 
%# e1 << e2 = e1*2^(e2);  
      j = n * 2 + 1;
      xsum = 1.0 / j;
      ysum = 0.0;
      for i=n:-1:1,
	j = j-2;
	xaux = (xsum * xquad - ysum * yquad) / i;
	ysum = (xsum * yquad + ysum * xquad) / i;
	xsum = xaux + 1.0 / j;
      end
      u1 = (xsum * yabs + ysum * xabs) * -M_2_SQRTPI + 1.0;
      v1 = (xsum * xabs - ysum * yabs) * M_2_SQRTPI;
      daux = exp(-xquad);
      u2 = daux * cos(yquad);
      v2 = -daux * sin(yquad);
      
      u = u1 * u2 - v1 * v2;
      v = u1 * v2 + v1 * u2;
    else
	
%       #  If (qrho > 1.0) then w(z) is evaluated using the Laplace continued 
% 	#  fraction.  nu is the minimum number of terms needed to obtain the
% 	#  required accuracy. 
%         #  if ((qrho > 0.085264) && (qrho < 1.0)) then w(z) is evaluated
%         #  by a truncated Taylor expansion, where the Laplace continued
% 	#  fraction is used to calculate the derivatives of w(z). 
%         #  kapn is the minimum number of terms in the Taylor expansion needed
%         #  to obtain the required accuracy. 
%         #  nu is the minimum number of terms of the continued fraction needed
%         #  to calculate the derivatives with the required accuracy. 
	
      if (qrho > 1.0),
	h = 0.0;
	kapn = 0;
	qrho = sqrt(qrho);
%# C++: (int)
	nu = fix(1442 / (qrho * 26 + 77) + 3);
      else
	qrho = (1 - y) * sqrt(1 - qrho);
	h = qrho * 1.88;
	h2 = h * 2;
	kapn = round(qrho * 34 + 7);
	nu   = round(qrho * 26 + 16);
      end
      
      b = h > 0.0;
      
      if (b),
	qlambda = h2.^kapn;
      end
      
      rx =  0.0;
      ry =  0.0;
      sx =  0.0;
      sy = 0.0;
      for n=nu:-1:0,
	np1 = n + 1;
	tx = yabs + h + np1 * rx;
	ty = xabs - np1 * ry;
	c = 0.5 / (tx * tx + ty * ty);
	rx = c * tx;
	ry = c * ty;
	if (b & (n <= kapn)),
	  tx = qlambda + sx;
	  sx = rx * tx - ry * sy;
	  sy = ry * tx + rx * sy;
	  qlambda = qlambda/h2;
	end
      end
      
      if (h == 0.0),
	u = rx * M_2_SQRTPI;
	v = ry * M_2_SQRTPI;
      else
	u = sx * M_2_SQRTPI;
	v = sy * M_2_SQRTPI;
      end
      
      if (yabs == 0.0),
	u = exp(-(xabs * xabs));
      end
    end
    
 %   #  evaluation of w(z) in the other quadrants 

    if (yi < 0.0),
      
      if (a),
	u2 = u2*2;
	v2 = v2*2;
      else
	xquad = -xquad;

      
  %          # the following statement protects 2*exp(-z**2) against overflow
	if ((yquad > rmaxgoni) | (xquad > rmaxexp)),
	  message('Faddeeva, absolute value of argument so large w(z) overflows');
          f(fila,columna)=NaN;
	end
      
	w1 = exp(xquad) * 2;
	u2 =  w1 * cos(yquad);
	v2 = -w1 * sin(yquad);
      end
    
      u = u2 - u;
      v = v2 - v;
      if (xi > 0.0),
	v = -v;
      end
    elseif (xi < 0.0),
      v = -v;
    end
    f(fila,columna)=u+sqrt(-1)*v;
  end
end