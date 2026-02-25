% Running sc.peace_pos to obtain the time intervals when leea measures the
% 180 deg pitch angle.

c_eval('angleHEEA? = sc.peace_pos(diB?,?,''heea'');',3:4)

c_eval('angle10? = find(angleHEEA?(:,2) < 10);',3:4);
c_eval('angle20? = find(angleHEEA?(:,2) < 20);',3:4); % this is ok

c_eval('jumpAngle10? = find(abs(diff(angle10?))>1);',3:4)
c_eval('jumpAngle20? = find(abs(diff(angle20?))>1);',3:4)

%tintIndices1 = reshape(diB3(jumpAngle103,1),);
%tintIndices2 = 