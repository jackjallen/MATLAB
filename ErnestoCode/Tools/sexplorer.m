function sexplorer(S)
  Sname = inputname(1);
  if ~any( strcmp( import , 'javax.swing.*' ) )
    import('javax.swing.*');
  end

  h = figure('Units', 'pixels', 'Position',[5    32   300   810],...
             'NumberTitle', 'off',...
             'NextPlot','new'                                ,...
             'IntegerHandle','off',...
             'name',['Structure: ' Sname],...
             'Toolbar', 'none',...
             'Menu','none');
           
  texto = uicontrol('Style','text','Parent',h,'Units','normalized',...
                    'Position',[0.02 .9 .96 .1],'string','',...
                    'FontSize',11,'HorizontalAlignment','left','FontWeight','bold');
           
  root = uitreenode( 'v0' , Sname , Sname , [] , 0 );

  tree = uitree( 'v0' , 'Parent', h , 'ExpandFcn', @ExpandNodeFcn, 'Root', root );
  drawnow;
  set(tree , 'Units', 'normalized'   ,...
             'position', [0 0 1 0.9] ,...
             'NodeSelectedCallback'  , @(h,e)NodeSelectedFcn(h,e)  );
  tree.expand(root);
%   ExpandAll(tree,root);

  try
  figurewindowstate(h,'ontop');
  end

  function NodeSelectedFcn( tree , ev )
    N   = ev.getCurrentNode.getValue;

    cl  = evalin('base', ['class(' N ')'] );
    sz  = num2str( evalin('base', ['size(' N ')'] ) );
    str = { [ N ' :   < ' cl ' >' ] ; [ 'size :  [ ' sz ' ] ' ] };
    set(texto,'String',str);

%     openvar( N );
    clc;
%     evalin( 'base', ['if numel(' N ')<100, disp(' N '); end']);
    evalin( 'base', ['disp(' N ')'] );
  end

  function ExpandAll( tree , N )
    tree.expand(N);
    drawnow;
    for nn= 1:N.getChildCount
      c= N.getChildAt(nn-1);
      ExpandAll( tree , c );
    end
  end

  function isleaf = IsLeaf( D )
    isleaf= [];
    try
        isleaf = fieldnames(D);
    end
    if isempty( isleaf )
      isleaf= 1;
    else
      isleaf= 0;
    end
    
    isleaf= 1;
    switch class(D)
      case 'cell'
        isleaf= 0;
      case 'struct'
        isleaf= 0;
    end
  end

  function nodes = ExpandNodeFcn( tree , N )
    D    = evalin('base', N );
    done = 0;

    if numel(D) > 1
      NN= findstr( N, '.');
      if ~isempty( NN )
        NN= N(NN(end)+1:end);
      else
        NN= N;
      end

      sz = size(D);
      if numel(sz) == 2
        if sz(1) == 1
          sz = sz(2);
        elseif sz(2) == 1
          sz = sz(1);
        end
      end
        
      for nn = 1:numel(D)
        if iscell( D )
          ind_str = {};
          [ind_str{1:numel(sz)}] = ind2sub(sz,nn);
          ind_str = sprintf('%d,',ind_str{:});
          ind_str(end) = [];
          
          
          nodes(nn) =  uitreenode('v0',[N '{', ind_str , '}'], [NN '{', ind_str ,'} <' class( D{nn} ) '>' ], '', 0);
        else
          nodes(nn) =  uitreenode('v0',[N '(', num2str(nn),')'], [NN '(', num2str(nn),')'], '', 0);
        end
        done = 1;
      end
    else
      fnames= [];
      try fnames = fieldnames(D);  end
      for nn=1:numel(fnames)
        nodes(nn) = uitreenode( 'v0' , [N '.' fnames{nn}], [ fnames{nn} ' <' class( D.(fnames{nn}) ) '>' ], '' , IsLeaf(D.(fnames{nn})) );
        done = 1;
      end
    end
    if ~done,      nodes = [];    end  
  end
end






% function sexplorer(S)
%   Sname= inputname(1);
%   if ~any( strcmp( import , 'javax.swing.*' ) )
%     import javax.swing.*;
%   end
% 
%   h = figure('Units', 'pixels', 'Position',[5 80 300 890],...
%              'NumberTitle', 'off',...
%              'name',['Structure: ' Sname],...
%              'Toolbar', 'none',...
%              'Menu','none');
%            
%   texto= uicontrol('Style','text','Parent',h,'Units','normalized',...
%                    'Position',[0.02 .9 .98 .1],'string','',...
%                    'FontSize',11,'HorizontalAlignment','left','FontWeight','bold');
%            
%   root = uitreenode(Sname, Sname, [], 0);
%   SetChildren( root );
% 
%   tree = uitree(h, 'Root', root );
%   set(tree , 'Units', 'normalized'   ,...
%              'position', [0 0 1 0.9] ,...
%              'NodeSelectedCallback'  , @NodeSelectedFcn  );
%            
%   ExpandAll(tree,root);
%   figureontop(h,'on');
% 
%   function NodeSelectedFcn( tree , ev )
%     N   = ev.getCurrentNode.getValue;
%     cl  = evalin('base', ['class(' N ')'] );
%     sz  = num2str( evalin('base', ['size(' N ')'] ) );
%     str = { [ N ' :   < ' cl ' >' ] ; [ 'size :  [ ' sz ' ] ' ] };
%     set(texto,'String',str);
% 
%     clc;
%     evalin( 'base', ['if numel(' N ')<100, disp(' N '); end']);
%   end
% 
%   function SetChildren( N )
%     D = evalin('base', N.getValue );
% 
%     if numel(D) > 1
%       for nn = 1:numel(D)
%         nN =  uitreenode( [N.getValue '(', num2str(nn),')'], [N.getValue '(', num2str(nn),')'], '', 0);
%         N.add( nN );
%         SetChildren( nN );
%       end
%     else
%       fnames= [];
%       try fnames = fieldnames(D); end
%       for nn=1:numel(fnames)
%         isleaf= [];
%         try isleaf = fieldnames(D.(fnames{nn}));  end
%         if isempty( isleaf ), isleaf= 1;
%         else                  isleaf= 0;   end
% 
%         nN = uitreenode( [N.getValue '.' fnames{nn}], fnames{nn}, '' , isleaf );
%         N.add( nN );
%         if ~isleaf
%           SetChildren( nN );
%         end
%       end
%     end
%   end
% 
%   function ExpandAll( tree , N )
%     tree.expand(N);
%     drawnow;
%     for nn= 1:N.getChildCount
%       c= N.getChildAt(nn-1);
%       ExpandAll( tree , c );
%     end
%   end
% 
% end
% 
% 

