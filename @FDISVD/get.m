function [ Sigma, V ] = get( obj, varargin )
  if ~isempty(varargin)
      r = varargin{1};
  else
      r = numel(obj.S);
  end
  try
    Sigma = obj.S(1:r);
    V = obj.V(:,1:r);
  catch ME
    fprintf("%s", ME.message);
  end

end