classdef NewClass
% NewClass

  properties 
    index = [];  
    data = [];
  end
  
  methods
    function obj = NewClass(data)
      dataSize = size(data);
      obj.index = 1:dataSize(1);
      obj.data = data;            
    end
    
    function [varargout] = subsref(obj,idx)      
      switch idx(1).type        
       case '.'
         [varargout{1:nargout}] = builtin('subsref',obj,idx);
       case '()'
         subIndex = builtin('subsref',obj.index,idx(1));
         obj.index = subIndex;
         obj.data = obj.data(subIndex,:); 
         [varargout{1:nargout}] = obj;
      end
    end  
    
    function s = size(obj)
      s = size(obj.data);
    end
  end
end