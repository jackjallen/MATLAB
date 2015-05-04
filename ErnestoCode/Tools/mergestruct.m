function A = mergestruct( A , B )

  if ~isstruct(A) || ~isstruct(B)
    error('mergestruct  struct_A  and  struct_B. Both inputs must be struct.');
  end
  if numel(A) ~= 1 || numel(B) ~= 1
    error('both struct must be 1x1.');
  end
  
  fn = fieldnames(B);
  for f=1:numel(fn)
    A.(fn{f}) = B.(fn{f});
  end

end
