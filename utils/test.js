var fourierExpansion = function (n,N, t, hl) {

   l = (hl * 1.7) / 2;
   var v = [];
  

  const result = animationData(n,N,t);
   let hx = result[0];
   let hy = result[1];
   let hz = result[2];
  
   let rx = result[3];
   let ry = result[4];
   let rz = result[5];

   for(var i=0; i<N+1;i++) {
    //rz[i]=rz[i]-.35;
    //rz[i]=rz[i]-.35;
   };

   kappa = result[6];
  
   let cen  =  [0,0,0];
 
  var mid = [0, 0, 0];

 
  // creating tetrahedron vertices
  

  for (var i = 1; i < N +1; i++) {
    

    mid[0] = rx[i];
    mid[1] = ry[i];
    mid[2] = rz[i];

    v[6 * i - 5] = mid[0] - hx[i] * l;
    v[6 * i - 4] = mid[1] - hy[i] * l;
    v[6 * i - 3] = mid[2] - hz[i] * l;
    v[6 * i - 2] = mid[0] + hx[i] * l;
    v[6 * i - 1] = mid[1] + hy[i] * l;
    v[6 * i]     = mid[2] + hz[i] * l;
  }
  for (var i = 0; i < 2*N+1; i++) {
    for (var j = 1; j < 4; j++) {
      v[3*i+j] -= cen[j-1]/n;
      if(j==3){
        v[3*i+j] -=.25;
      }
    };
  };
 
  return [v, rx];
};





function updateBandMesh(n,N,t, hingeLength) {
          
    let temp = fourierExpansion(n,N, t, hingeLength);
    var v = temp[0];
 

 

    // N = 105;
    var pos = [];
    var pathMidline = [];
    var indices = [];
    var indicesBack = [];
    var colors = [];
    var colorsBack = [];
    
    var sign_band = (-1)**n;
    for (var i = 1; i < N + 1; i++) {
      var v1 = [v[6 * i - 5], v[6 * i - 3], -v[6 * i - 4]];
      var v2 = [v[6 * i - 2], v[6 * i], -v[6 * i - 1]];

 
      // surface 

      pos.push(v1[0], v1[1], v1[2]);
      pos.push(v2[0], v2[1], v2[2]);

       

      if (i > 0) {
        indices.push(2 * i - 1, 2 * i, 2 * i - 2);
        indices.push(2 * i, 2 * i - 1, 2 * i + 1);

        indicesBack.push(2 * i, 2 * i - 1, 2 * i - 2);
        indicesBack.push(2 * i - 1, 2 * i, 2 * i + 1);

      };
    };

    //   console.log(pathMidline)
    // Update the tube mesh with the new path
    // Remove the existing tube mesh


    i = 1;
    var v2 = [v[6 * i - 5], v[6 * i - 3], -v[6 * i - 4]];
    var v1 = [v[6 * i - 2], v[6 * i], -v[6 * i - 1]];

 
    pos.push(v1[0], v1[1], v1[2]);
    pos.push(v2[0], v2[1], v2[2]);

   

    var normals = [];
    BABYLON.VertexData.ComputeNormals(pos, indices, normals);

    var vertexData = new BABYLON.VertexData();
    vertexData.positions = pos;
    vertexData.indices = indices;
    vertexData.normals = normals;
    vertexData.colors = colors;

    vertexData.applyToMesh(band);



    var normalsBack = [];
    BABYLON.VertexData.ComputeNormals(pos, indicesBack, normalsBack);
    vertexDataBack = new BABYLON.VertexData();
    vertexDataBack.positions = pos;
    vertexDataBack.indices = indicesBack;
    vertexDataBack.normals = normalsBack;
 
    vertexDataBack.applyToMesh(bandBack);

    ///  surface created 

     
    };