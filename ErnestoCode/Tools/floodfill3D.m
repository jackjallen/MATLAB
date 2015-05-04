function I= floodfill3D( I , i,j,k, conn , plane )

  if conn==4 || conn==8
    switch upper(plane)
      case {'XY','YX','Z'}
        im = I(:,:,k);
        if ~im(i,j), im= ~im; end
        im = bwlabel( im , conn );
        im = im==im(i,j);
        I  = logical( I*0 );
        I(:,:,k) = im;
      case {'XZ','ZX','Y'}
        im = squeeze( I(:,j,:) );
        if ~im(i,k), im= ~im; end
        im = bwlabel( im , conn );
        im = im==im(i,k);
        I  = logical( I*0 );
        I(:,j,:) = im;
      case {'YZ','ZY','X'}
        im = squeeze( I(i,:,:) );
        if ~im(j,k), im= ~im; end
        im = bwlabel( im , conn );
        im = im==im(j,k);
        I  = logical( I*0 );
        I(i,:,:) = im;
    end
  elseif conn==6 || conn==18 || conn==26
    if ~I(i,j,k), I= ~I; end
    I= bwlabeln( I , conn);
    I= I==I(i,j,k);
  end

end
