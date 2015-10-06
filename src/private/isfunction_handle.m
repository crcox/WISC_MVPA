function b = isfunction_handle(x)
  if isa(x, 'function_handle');
    b = true;
  else
    b = false;
  end
end
