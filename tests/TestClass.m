classdef TestClass
  % My test class
  properties 
    index
    data   
  end
  
  methods
    function obj = TestClass(data) % constructor    
      obj.index = 1:numel(data);
      obj.data = data(:);
    end    
    function varargout = subsref(obj,idx)
    %SUBSREF handle indexing

    switch idx(1).type           
      case '()'                
        tmpIndex = idx(1).subs{1};
        obj.index = obj.index(tmpIndex);
        obj.data = obj.data(tmpIndex); 
        
        if numel(idx) > 1 
          obj = builtin('subsref',obj,idx(2:end));
        end     

        [varargout{1:nargout}] = obj;                
    end
    end
    function varargout = test_varargout(obj)
      nargout
      if nargout == 1
        varargout{1} = 1;
      elseif nargout == 2
        varargout{1} = 1;
        varargout{2} = 2;
      end
    end        
  end  
end