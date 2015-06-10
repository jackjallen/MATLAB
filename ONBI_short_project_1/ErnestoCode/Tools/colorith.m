function c = colorith( i )

  c = rand(1,3);
  switch i
    case 1, c = [1 0 0];
    case 2, c = [0 1 0];
    case 3, c = [0 0 1];
    case 4, c = [1 0 1];
    case 5, c = [0 1 1];
    case 6, c = [0 0 0];
  end

end