function cdup( c )

  if nargin < 1, c = 1; end
  if ischar(c), c = str2double(c); end

  D = com.mathworks.mde.desk.MLDesktop.getInstance;

  D = D.getMainFrame.getComponent(0).getComponent(1).getComponent(1).getComponent(0).getComponent(0);

  for i = (1:D.getComponentCount)-1
    if strcmp( D.getComponent(i).getName , 'Current Directory History' )
      break
    end
  end
  D = D.getComponent(i);

  cd( D.getItemAt(c) )

end
