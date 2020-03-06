! Two-dimensional periodic boundary conditions

       array(:,0) = array(:,ny)
       array(:,ny+1) = array(:,1)
       array(nx+1,:) = array(1,:)
       array(0,:) = array(nx,:)
