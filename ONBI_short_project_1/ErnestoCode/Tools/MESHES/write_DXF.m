function write_DXF( mesh , filename )

  fields= fieldnames( mesh );

  f= fopen(filename,'w');
  
    %TABLES SECTION
    fprintf(f,'0\nSECTION\n2\nTABLES\n');
      %LTYPE table
      fprintf(f,'0\nTABLE\n2\nLTYPE\n0\nLTYPE\n2\nCONTINUOUS\n70\n0\n3\nSolid\n72\n65\n73\n0\n40\n0\n0\nENDTAB\n');
      %LAYERS table    
      fprintf(f,'0\nTABLE\n2\nLAYER\n');
        %LAYER 0
        fprintf(f,'0\nLAYER\n2\n0\n70\n0\n62\n7\n6\nCONTINUOUS\n');
        
        for ff=1:size(fields,1)
          field= fields{ff};
          if ( strncmp( field, 'lmk',3) )
            fprintf(f,'0\nLAYER\n2\n%s\n70\n0\n62\n%d\n6\nCONTINUOUS\n', field,ff);
          end
        end
    
        %LAYER XYZNORMALS table
        fprintf(f,'0\nLAYER\n2\nXYZNORMALS\n70\n0\n62\n4\n6\nCONTINUOUS\n');
        %LAYER TRINORMALS table
        fprintf(f,'0\nLAYER\n2\nTRINORMALS\n70\n0\n62\n6\n6\nCONTINUOUS\n');
        %LAYER NODES table
        fprintf(f,'0\nLAYER\n2\nNODES\n70\n0\n62\n6\n6\nCONTINUOUS\n');
        %LAYER CELLS table
        fprintf(f,'0\nLAYER\n2\nCELLS\n70\n0\n62\n3\n6\nCONTINUOUS\n');
      %END LAYERS table
      fprintf(f,'0\nENDTAB\n');
      %STYLE table
      fprintf(f,'0\nTABLE\n2\nSTYLE\n0\nSTYLE\n2\nSTANDARD\n70\n0\n40\n0\n41\n1\n50\n0\n71\n0\n42\n1\n3\nromans.shx\n0\nENDTAB\n');
    %END TABLES SECTION
    fprintf(f,'0\nENDSEC\n');
    
    %BLOCKS SECTION
    fprintf(f,'0\nSECTION\n2\nBLOCKS\n');
      %LMK block
      fprintf(f,'0\nBLOCK\n8\n0\n2\nLMK\n70\n2\n10\n0\n20\n0\n30\n0\n');
        fprintf(f,'0\nATTDEF\n8\n0\n62\n0\n10\n0\n20\n0\n30\n0\n40\n1\n1\n0\n11\n0\n21\n0\n31\n0\n3\nLMK_ID\n2\nLMK_ID\n70\n0\n');
        fprintf(f,'0\nCIRCLE\n8\n0\n62\n0\n10\n0\n20\n0\n30\n0\n40\n0.95\n');
        fprintf(f,'0\nCIRCLE\n8\n0\n62\n0\n10\n0\n20\n0\n30\n0\n40\n1.05\n');
      fprintf(f,'0\nENDBLK\n');
      %CELL block
      fprintf(f,'0\nBLOCK\n8\n0\n2\nCELL\n70\n2\n10\n0\n20\n0\n30\n0\n');
        fprintf(f,'0\nATTDEF\n8\n0\n62\n0\n10\n0\n20\n0\n30\n0\n40\n1\n1\n0\n11\n0\n21\n0\n31\n0\n3\nCELL_ID\n2\nCELL_ID\n70\n0\n');
        fprintf(f,'0\nLINE\n8\n0\n62\n0\n10\n1\n20\n1\n30\n0\n11\n-1\n21\n1\n31\n0\n');
        fprintf(f,'0\nLINE\n8\n0\n62\n0\n10\n-1\n20\n1\n30\n0\n11\n-1\n21\n-1\n31\n0\n');
        fprintf(f,'0\nLINE\n8\n0\n62\n0\n10\n-1\n20\n-1\n30\n0\n11\n1\n21\n-1\n31\n0\n');
        fprintf(f,'0\nLINE\n8\n0\n62\n0\n10\n1\n20\n-1\n30\n0\n11\n1\n21\n1\n31\n0\n');
      fprintf(f,'0\nENDBLK\n');
      %NODE block
      fprintf(f,'0\nBLOCK\n8\n0\n2\nNODE\n70\n2\n10\n0\n20\n0\n30\n0\n');
        fprintf(f,'0\nATTDEF\n8\n0\n62\n0\n10\n0\n20\n0\n30\n0\n40\n1\n1\n0\n11\n0\n21\n0\n31\n0\n3\nNODE_ID\n2\nNODE_ID\n70\n0\n');
        fprintf(f,'0\nLINE\n8\n0\n62\n0\n10\n1\n20\n0\n30\n0\n11\n0\n21\n1\n31\n0\n');
        fprintf(f,'0\nLINE\n8\n0\n62\n0\n10\n0\n20\n1\n30\n0\n11\n-1\n21\n0\n31\n0\n');
        fprintf(f,'0\nLINE\n8\n0\n62\n0\n10\n-1\n20\n0\n30\n0\n11\n0\n21\n-1\n31\n0\n');
        fprintf(f,'0\nLINE\n8\n0\n62\n0\n10\n0\n20\n-1\n30\n0\n11\n1\n21\n0\n31\n0\n');
      fprintf(f,'0\nENDBLK\n');
      %NORMAL block
      fprintf(f,'0\nBLOCK\n8\n0\n2\nNORMAL\n70\n0\n10\n0\n20\n0\n30\n0\n');
        fprintf(f,'0\nPOLYLINE\n8\n0\n62\n0\n66\n1\n10\n0\n20\n0\n30\n0\n210\n1\n220\n0\n230\n0\n');
          fprintf(f,'0\nVERTEX\n8\n0\n62\n0\n10\n0\n20\n0\n30\n0\n');
          fprintf(f,'0\nVERTEX\n8\n0\n62\n0\n10\n0\n20\n0.7\n30\n0\n40\n0.15\n');
          fprintf(f,'0\nVERTEX\n8\n0\n62\n0\n10\n0\n20\n1\n30\n0\n');
        fprintf(f,'0\nSEQEND\n');
        fprintf(f,'0\nPOLYLINE\n8\n0\n62\n0\n66\n1\n10\n0\n20\n0\n30\n0\n210\n0\n220\n1\n230\n0\n');
          fprintf(f,'0\nVERTEX\n8\n0\n62\n0\n10\n0\n20\n0\n30\n0\n');
          fprintf(f,'0\nVERTEX\n8\n0\n62\n0\n10\n0\n20\n0.7\n30\n0\n40\n0.15\n');
          fprintf(f,'0\nVERTEX\n8\n0\n62\n0\n10\n0\n20\n1\n30\n0\n');
        fprintf(f,'0\nSEQEND\n');
        fprintf(f,'0\nPOLYLINE\n8\n0\n62\n0\n66\n1\n10\n0\n20\n0\n30\n0\n210\n-0.7071067811865475\n220\n0.7071067811865475\n230\n0\n');
          fprintf(f,'0\nVERTEX\n8\n0\n62\n0\n10\n0\n20\n0\n30\n0\n');
          fprintf(f,'0\nVERTEX\n8\n0\n62\n0\n10\n0\n20\n0.7\n30\n0\n40\n0.15\n');
          fprintf(f,'0\nVERTEX\n8\n0\n62\n0\n10\n0\n20\n1\n30\n0\n');
        fprintf(f,'0\nSEQEND\n');
        fprintf(f,'0\nPOLYLINE\n8\n0\n62\n0\n66\n1\n10\n0\n20\n0\n30\n0\n210\n-0.7071067811865475\n220\n-0.7071067811865475\n230\n0\n');
          fprintf(f,'0\nVERTEX\n8\n0\n62\n0\n10\n0\n20\n0\n30\n0\n');
          fprintf(f,'0\nVERTEX\n8\n0\n62\n0\n10\n0\n20\n0.7\n30\n0\n40\n0.15\n');
          fprintf(f,'0\nVERTEX\n8\n0\n62\n0\n10\n0\n20\n1\n30\n0\n');
        fprintf(f,'0\nSEQEND\n');
      fprintf(f,'0\nENDBLK\n');
    %END BLOCKS SECTION
    fprintf(f,'0\nENDSEC\n');
    
    %ENTITIES SECTION
    fprintf(f,'0\nSECTION\n2\nENTITIES\n');
    
      %FACES
      xyz= [ mesh.xyz(mesh.tri(:,1),:) mesh.xyz(mesh.tri(:,2),:) mesh.xyz(mesh.tri(:,3),:)];
      fprintf(f,'0\n3DFACE\n8\n0\n10\n%f\n20\n%f\n30\n%f\n11\n%f\n21\n%f\n31\n%f\n12\n%f\n22\n%f\n32\n%f\n13\n%f\n23\n%f\n33\n%f\n', ...
          xyz(:,[1:9 7:9])' );  %'
    
      %LANDMARKS
      for ff=1:size(fields,1)
        field= fields{ff};
        if ( strncmp( field, 'lmk',3) )
          lmks= getfield(mesh,field);
          for i=1:size( lmks , 1)
            fprintf(f,'0\nINSERT\n8\n%s\n66\n1\n2\nLMK\n10\n%f\n20\n%f\n30\n%f\n0\nATTRIB\n8\n0\n62\n0\n10\n%f\n20\n%f\n30\n%f\n40\n1.0\n1\n%d\n72\n4\n11\n%f\n21\n%f\n31\n%f\n2\nLMK_ID\n70\n0\n0\nSEQEND\n',...
               field, lmks(i,:), lmks(i,:), i, lmks(i,:) );
          end
        end
      end
        
      %NODE_ID
      xyz= mesh.xyz;
      xyz(:,4)= [1:size(xyz,1)]';
      fprintf(f,'0\nINSERT\n8\nNODES\n66\n1\n2\nNODE\n10\n%f\n20\n%f\n30\n%f\n0\nATTRIB\n8\n0\n62\n0\n10\n%f\n20\n%f\n30\n%f\n40\n1.0\n1\n%d\n72\n4\n11\n%f\n21\n%f\n31\n%f\n2\nNODE_ID\n70\n0\n0\nSEQEND\n',...
        xyz(:,[1:3 1:3 4 1:3])');       %'
    
      %CELL_ID
      xyz= [ mesh.xyz(mesh.tri(:,1),:) mesh.xyz(mesh.tri(:,2),:) mesh.xyz(mesh.tri(:,3),:)];
      xyz= [ sum(xyz(:,[1 4 7]),2)/3  sum(xyz(:,[2 5 8]),2)/3  sum(xyz(:,[3 6 9]),2)/3 ];
      xyz(:,4)= [1:size(xyz,1)]';
      fprintf(f,'0\nINSERT\n8\nCELLS\n66\n1\n2\nCELL\n10\n%f\n20\n%f\n30\n%f\n0\nATTRIB\n8\n0\n62\n0\n10\n%f\n20\n%f\n30\n%f\n40\n1.0\n1\n%d\n72\n4\n11\n%f\n21\n%f\n31\n%f\n2\nCELL_ID\n70\n0\n0\nSEQEND\n',...
        xyz(:,[1:3 1:3 4 1:3])');  %'
    
      %CELLS_NORMALS
      xyz= [ mesh.xyz(mesh.tri(:,1),:) mesh.xyz(mesh.tri(:,2),:) mesh.xyz(mesh.tri(:,3),:)];
      xyz= [ sum(xyz(:,[1 4 7]),2)/3  sum(xyz(:,[2 5 8]),2)/3  sum(xyz(:,[3 6 9]),2)/3 ];
      if isfield(mesh, 'triNORMALS')
        normals= mesh.triNORMALS;
      else
        normals= [ mesh.xyz(mesh.tri(:,1),:) mesh.xyz(mesh.tri(:,2),:) mesh.xyz(mesh.tri(:,3),:)];
        normals= [ normals(:,4:6)-normals(:,1:3) normals(:,7:9)-normals(:,1:3) ];
        normals= [ normals(:,2).*normals(:,6)-normals(:,3).*normals(:,5) ...
                   normals(:,3).*normals(:,4)-normals(:,6).*normals(:,1) ...
                   normals(:,1).*normals(:,5)-normals(:,2).*normals(:,4) ];    
      end
      normals= normals./repmat( sqrt( sum( normals.^2,2)),1,3);
    
      for i=1:size(normals,1)
        if ( abs(normals(i,1)) < 1/64 ) & (abs(normals(i,2)) < 1/64)
        	Ax = cross([0 1 0],normals(i,:));
        else
        	Ax = cross([0 0 1],normals(i,:));
        end
        Ax = Ax/sqrt( sum(Ax.^2,2) );
        Ay = cross(normals(i,:),Ax);
        Ay = Ay/sqrt( sum(Ay.^2,2) );
    
        fprintf(f,'0\nINSERT\n8\nTRINORMALS\n2\nNORMAL\n10\n%f\n20\n%f\n30\n%f\n50\n0\n210\n%f\n220\n%f\n230\n%f\n', ...
          xyz(i,:)*[Ax ; Ay ; normals(i,:)]' , normals(i,:) );  %'
      end    
    
    
      %NODE_NORMALS
      if isfield(mesh, 'xyzNORMALS')
        xyz= mesh.xyz;
        normals= mesh.xyzNORMALS;
        normals= normals./repmat( sqrt( sum( normals.^2,2)),1,3);
    
        for i=1:size(normals,1)
          if ( abs(normals(i,1)) < 1/64 ) & (abs(normals(i,2)) < 1/64)
          	Ax = cross([0 1 0],normals(i,:));
          else
          	Ax = cross([0 0 1],normals(i,:));
          end
          Ax = Ax/sqrt( sum(Ax.^2,2) );
          Ay = cross(normals(i,:),Ax);
          Ay = Ay/sqrt( sum(Ay.^2,2) );
      
          fprintf(f,'0\nINSERT\n8\nXYZNORMALS\n2\nNORMAL\n10\n%f\n20\n%f\n30\n%f\n50\n0\n210\n%f\n220\n%f\n230\n%f\n', ...
            xyz(i,:)*[Ax ; Ay ; normals(i,:)]' , normals(i,:) );  %'
        end    
      end
    
    
    %END ENTITIES SECTION
    fprintf(f,'0\nENDSEC\n');
  
  %END FILE
  fprintf(f,'0\nEOF\n');
  fclose(f);
  
  
  
  
  
