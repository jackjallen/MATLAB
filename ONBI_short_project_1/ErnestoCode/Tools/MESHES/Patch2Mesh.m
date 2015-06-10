function M = Patch2Mesh( P )
  
  if isstruct(P)
    try,    M.tri = P.faces;                 end
    try,    M.xyz = P.vertices;              end
    try,    M.xyzDATA = P.facevertexcdata;   end
  elseif ishandle( P )
    try, P = findall( P , 'type','patch' ); end
    if numel(P)>1
      disp('tengo que implementarlo');
    elseif numel(P)==1
      switch lower(get(P,'type'))
        case 'patch'
          M.xyz= get(P,'vertices');
          M.tri= get(P,'faces');
        
      end
    end
  end
  
end







% 
% 
%             AlphaDataMapping: 'scaled'
%                        CData: []
%                    FaceAlpha: 0.5000
%                    FaceColor: [0.8000 0.4000 0.5000]
%                        Faces: [5628x3 double]
%          FaceVertexAlphaData: []
%              FaceVertexCData: []
%                    LineWidth: 0.5000
%                       Marker: 'none'
%              MarkerEdgeColor: 'auto'
%              MarkerFaceColor: 'none'
%                   MarkerSize: 6
%                VertexNormals: [2826x3 double]
%                     Vertices: [2826x3 double]
%                        XData: [3x5628 double]
%                        YData: [3x5628 double]
%                        ZData: [3x5628 double]
