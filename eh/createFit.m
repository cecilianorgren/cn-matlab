function [fitresult, gof, fun_out, fun_deriv_out] = createFit(phi, dntrap)
%CREATEFIT(PHI,DNTRAP)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : phi
%      Y Output: dntrap
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 17-Aug-2018 17:40:05


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( phi, dntrap );
xData = tocolumn(xData);
yData = tocolumn(yData);
% Set up fittype and options.
fit_type = 'poly8';
switch fit_type
  case 'poly8'
    ft = fittype( 'poly8' );

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft );

    fun_eval_str = sprintf(...
      'fun_out = @(x) %g*x.^8 + %g*x.^7 + %g*x.^6 + %g*x.^5 + %g*x.^4 + %g*x.^3 + %g*x.^2 + %g*x + %g;',...
      fitresult.p1,fitresult.p2,fitresult.p3,fitresult.p4,fitresult.p5,fitresult.p6,fitresult.p7,fitresult.p8,fitresult.p9);
    eval(fun_eval_str)

    fun_eval_str_deriv = sprintf(...
      'fun_deriv_out = @(x) 8*%g*x.^(8-1) + 7*%g*x.^(7-1) + 6*%g*x.^(6-1) + 5*%g*x.^(5-1) + 4*%g*x.^(4-1) + 3*%g*x.^(3-1) + 2*%g*x.^(2-1) + 1*%g*x.^(1-1) + 0*%g;',...
      fitresult.p1,fitresult.p2,fitresult.p3,fitresult.p4,fitresult.p5,fitresult.p6,fitresult.p7,fitresult.p8,fitresult.p9);
    eval(fun_eval_str_deriv)
  case 'poly7'
    ft = fittype( 'poly7' );

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft );

    fun_eval_str = sprintf(...
      'fun_out = @(x) %g*x.^7 + %g*x.^6 + %g*x.^5 + %g*x.^4 + %g*x.^3 + %g*x.^2 + %g*x + %g;',...
      fitresult.p1,fitresult.p2,fitresult.p3,fitresult.p4,fitresult.p5,fitresult.p6,fitresult.p7,fitresult.p8);
    eval(fun_eval_str)

    fun_eval_str_deriv = sprintf(...
      'fun_deriv_out = @(x) 7*%g*x.^(7-1) + 6*%g*x.^(6-1) + 5*%g*x.^(5-1) + 4*%g*x.^(4-1) + 3*%g*x.^(3-1) + 2*%g*x.^(2-1) + 1*%g*x.^(1-1) + 0*%g;',...
      fitresult.p1,fitresult.p2,fitresult.p3,fitresult.p4,fitresult.p5,fitresult.p6,fitresult.p7,fitresult.p8);
    eval(fun_eval_str_deriv)
  case 'poly6'
    ft = fittype( 'poly6' );

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft );

    fun_eval_str = sprintf(...
      'fun_out = @(x) %g*x.^6 + %g*x.^5 + %g*x.^4 + %g*x.^3 + %g*x.^2 + %g*x + %g;',...
      fitresult.p1,fitresult.p2,fitresult.p3,fitresult.p4,fitresult.p5,fitresult.p6,fitresult.p7);
    eval(fun_eval_str)

    fun_eval_str_deriv = sprintf(...
      'fun_deriv_out = @(x) 6*%g*x.^(6-1) + 5*%g*x.^(5-1) + 4*%g*x.^(4-1) + 3*%g*x.^(3-1) + 2*%g*x.^(2-1) + 1*%g*x.^(1-1) + 0*%g;',...
      fitresult.p1,fitresult.p2,fitresult.p3,fitresult.p4,fitresult.p5,fitresult.p6,fitresult.p7);
    eval(fun_eval_str_deriv)
  case 'poly5'
    ft = fittype( 'poly5' );

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft );

    fun_eval_str = sprintf(...
      'fun_out = @(x) %g*x.^5 + %g*x.^4 + %g*x.^3 + %g*x.^2 + %g*x + %g;',...
      fitresult.p1,fitresult.p2,fitresult.p3,fitresult.p4,fitresult.p5,fitresult.p6);
    eval(fun_eval_str)

    fun_eval_str_deriv = sprintf(...
      'fun_deriv_out = @(x) 5*%g*x.^(5-1) + 4*%g*x.^(4-1) + 3*%g*x.^(3-1) + 2*%g*x.^(2-1) + 1*%g*x.^(1-1) + 0*%g;',...
      fitresult.p1,fitresult.p2,fitresult.p3,fitresult.p4,fitresult.p5,fitresult.p6);
    eval(fun_eval_str_deriv)
  case 'poly4'
    ft = fittype( 'poly4' );

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft );

    fun_eval_str = sprintf(...
      'fun_out = @(x) %g*x.^4 + %g*x.^3 + %g*x.^2 + %g*x + %g;',...
      fitresult.p1,fitresult.p2,fitresult.p3,fitresult.p4,fitresult.p5);
    eval(fun_eval_str)

    fun_eval_str_deriv = sprintf(...
      'fun_deriv_out = @(x) 4*%g*x.^(4-1) + 3*%g*x.^(3-1) + 2*%g*x.^(2-1) + 1*%g*x.^(1-1) + 0*%g;',...
      fitresult.p1,fitresult.p2,fitresult.p3,fitresult.p4,fitresult.p5);
    eval(fun_eval_str_deriv)
  case 'poly2'
    ft = fittype( 'poly2' );

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft );

    fun_eval_str = sprintf(...
      'fun_out = @(x) %g*x.^2 + %g*x + %g;',...
      fitresult.p1,fitresult.p2,fitresult.p3);
    eval(fun_eval_str)

    fun_eval_str_deriv = sprintf(...
      'fun_deriv_out = @(x) 2*%g*x.^(2-1) + 1*%g*x.^(1-1) + 0*%g;',...
      fitresult.p1,fitresult.p2,fitresult.p3);
    eval(fun_eval_str_deriv)
end
% Plot fit with data.
if 0
  %figure( 'Name', 'untitled fit 1' );
  figure(13);
  h = plot( fitresult, xData, yData );
  hca = gca;
  legend( h, 'dntrap vs. phi', 'untitled fit 1', 'Location', 'NorthEast' );
  % Label axes
  hca.XLabel.String = 'phi';
  hca.YLabel.String = 'dntrap';  
  grid on
end
end
