#include "SpectralModularity.h"

/*
  Use of the stack for memory allocation. This 
 should be faster for large networks, but will need to reset NSIZE 
 large enough for your network size, and then re-Make. 

 */


SpectralModularity::SpectralModularity() { 

  this->A       = 0;
  this->Bgi     = 0;
  this->NR_Bgi  = 0;
  this->NC_Bgi  = 0;
  this->M       = 0;
  this->usedBgi = false;
  this->PRINT   = false;

  this->specQ   = 0;
  this->NORM    = 0; 
  
  this->tol     = 0.00001;//the tolerance value, 10^-5; eigenvalues below this threshold are not used
  this->MINCn   = 1;//The minimum cluster size

  this->MAXK    = 0;//Counter storing the maximum community number so far
  
}

SpectralModularity::SpectralModularity( network *gg, edgelist *el, double *A, int N, int M ) {

  this->gg      = gg;
  this->A       = A;
  this->Bgi     = 0;
  this->NR_Bgi  = N;//set initial rows for Bgi
  this->NC_Bgi  = N;//set initial cols for Bgi
  this->M       = M;//number of edges
  this->usedBgi = false;
  this->PRINT   = false;
  
  this->specQ   = 0;
  this->NORM    = 0; 

  this->tol     = 0.00001;//the tolerance value, 10^-5; eigenvalues below this threshold are not used
  this->MINCn   = 1;//The minimum cluster size

  this->MAXK    = 0;//Counter storing the maximum community number so far
   
  assignSpace();
  
  setupMatrices();
   
}

int SpectralModularity::calculateSpectralModularity(){

  int k,Ng,KK;

  Ng = NR_Bgi;
  KK = Ng * Ng;
  
  //--- Calculate eigenvectors, and values, from Bgi...
  calculateEigenVectors();

  
  if( PRINT ){ cout << "> max EigenValue is " << betai << endl; }
  
  //--- set up the index vectors, si and SI, for the initial split
  maximiseIndexVectors();


  //--- Calculate the Spectral Modularity
  double deltaQ_old = 0.0;
  double deltaQ_new = 0.0;
  
  deltaModularity(deltaQ_old);
  if( PRINT ){ cout << "> Spectral Q: " << deltaQ_old << endl; }

  double diff = deltaQ_old;

  //--- Fine tuning stage to maximum deltaModularity for the initial split
  visited[Ng];
  
  //reset visited
  for(k=0; k<Ng; k++){ visited[k]=0; }
  while( diff > tol ){

    modifySplit( Ng );
    
    deltaModularity( deltaQ_new );
    if( PRINT ){ cout << "> Modified Q: " << deltaQ_new << endl; }    

    if( deltaQ_new < 0 ) break;
    
    diff = fabs( deltaQ_new - deltaQ_old ); 

    deltaQ_old = deltaQ_new;    

  }

  //--- Keep recorded of maximum fine-tuned Modularity value.
  specQ += deltaQ_old;

  if( PRINT ){ cout << "si[0] " << si[0] << endl; }
  
  if( PRINT ){ cout << "> node list " << endl; }
  for(k=0; k<Ng; k++){

    if(si[k] > 0){
      keys_p[k]  = gg->V[k].id;
      keys_n[k]  = dummy;
      gg->V[k].K = 1;
    } else {
      keys_p[k]  = dummy;
      keys_n[k]  = gg->V[k].id;
      gg->V[k].K = 2;
    }

  }

  //--- update the maximum community number
  MAXK = 2;
  
  if( PRINT ){ gg->printVertices(); }

  
  //--- Recursively split the group of positive eigenvector nodes
  split( Bgi_temp, Ng, keys_p, "splitP" );
  
  //--- Recursively split the group of negative eigenvector nodes
  split( Bgi_temp, Ng, keys_n, "splitN" );

  
  if( PRINT ){ cout << "done." << endl; }
  
  
  return 0;
}

SpectralModularity::~SpectralModularity(){ freeSpace(); }

void SpectralModularity::freeSpace(){

  if( Bgi     !=0 && usedBgi  == true ){ free(Bgi);   }   
 
}


void SpectralModularity::assignSpace(){

  int i,k,Ng,KK; 

  Ng = NR_Bgi;
  KK = Ng * Ng;
  
  Bgi       = (double*)malloc(KK*sizeof(double));
  usedBgi   = true;

  //make a copy of Bgi to pass to splitP/N
  Bgi_temp[KK];
  
  for(k=0; k<KK; k++){
    Bgi[k]      = 0.0;
    Bgi_temp[k] = 0.0;
  }

  keys_p[Ng];
  keys_n[Ng];
  
  for(k=0; k<Ng; k++){
    keys_p[k] = 0;
    keys_n[k] = 0;
  }
  
}

/*
 Utility method used by Geodesic and RandomWalk algorithms
 to set-up the Modularity and Laplacian matrices.
 */
void SpectralModularity::setupMatrices(){

  int i,j,k,KK,Ng;
    
  Ng = NR_Bgi;
  KK = Ng * Ng;
  
  //---norm
  NORM = 1.0/(2.0*(double)M);

  if( PRINT ){ cout << "N: " << NR_Bgi << ", M: " << M << ", NORM: " << NORM << endl; }
  double sum=0.0;

  //--- The Modularity matrix, Bgi  
  for(k=0; k<KK; k++){
    i = floor(k/Ng);
    j = k % Ng;

    double val = A[(i*Ng)+j] - (gg->V[i].degree * gg->V[j].degree * NORM);
    
    Bgi[(i*Ng)+j]      = val;
    Bgi_temp[(i*Ng)+j] = val;

    sum += val;
    
  }

  if( PRINT ){ cout << "Sum(Bgi): " << sum << endl; }
  
}

/*
 Calculate the eigenvalues, betai, and eigenvectors, u, for 
 the current Modularity matrix Bgi.
 */
void SpectralModularity::calculateEigenVectors(){

  int i,j,k,KK,Ng,indx_max;
  
  Ng = NR_Bgi;
  KK = Ng*Ng;

  betai = 0.0; //hold the largest eigenvalue
  u[Ng];       //hold the eigenvector corresponding to largest eigenvalue
  

  //make a copy of Bgi to pass to gsl_matrix_view
  //--- in nr3.h
  MatDoub temp = MatDoub(Ng,Ng);

  for( k=0; k<KK; k++ ){

    //--- the linerar indexing for symmetric matrix
    i = floor( k/Ng );
    j = k % Ng;

    temp[i][j] = 0.0;
    temp[i][j] = Bgi[k];

  }

  //--- in eigen_sym.h
  Symmeig h(temp, true);
  
  //We just need the leading eigenvalue and vector.
  indx_max = 0;
  betai    = h.d[indx_max];
  for( i=0; i<Ng; i++ ){
    if( h.d[i] > betai ){
      betai = h.d[i];
      indx_max = i;
    }
  }
    
  for(k=0; k<Ng; k++){
    u[k] = 0.0;
    u[k] = h.z[k][indx_max];
  }
  //---
   
}


/*
 Update the index vectors, si and SI, for each node in the 
 current split such that:

 si(i) =  1 if eigenvector_max(i) > 0
       = -1 if eigenvector_max(i) < 0

 SI(i,0) = 1        
 SI(i,1) = 0 if eigenvector_max(i) > 0
         = 0
         = 1 if eigenvector_max(i) < 0
 */
void SpectralModularity::maximiseIndexVectors(){

  int i,j,k,KK,Ng;
  
  Ng = NR_Bgi;
  KK = 2*Ng;

  //set size of index vectors
  si[Ng];
  SI[KK]; 

  for(k=0; k<Ng; k++){

    if( u[k] < 0 ){
      si[k] = -1;
    } else {
      si[k] = 1;
    }

  }

  for(k=0; k<KK; k++){
    i = floor(k/2);
    j = k % 2;

    if( u[i] < 0 ){
      if( j == 0 ){ SI[(i*2)+j] = 0; }
      if( j == 1 ){ SI[(i*2)+j] = 1; }
    } else {
      if( j == 0 ){ SI[(i*2)+j] = 1; }
      if( j == 1 ){ SI[(i*2)+j] = 0; }
    }

  }


}

/*
 The change in Modularity used for the Spectral method.
 deltaQ = Sum_k { Sum_ij { si_ki * Bgi_ij * si_jk } }
*/
void SpectralModularity::deltaModularity( double &mod ){

  mod         = 0;

  int i,j,J,k,KK,Ng;
  double sum, sum1, sum2;

  Ng = NR_Bgi;
  KK = Ng*Ng;

  sum = 0.0; sum1 = 0.0; sum2 = 0.0;

  double SIt[(2*Ng)];

  for(i=0; i<Ng; i++){
    SIt[(i*2)+0] = 0;
    SIt[(i*2)+1] = 0;
    sum1         = 0;
    sum2         = 0;
    for(j=0; j<Ng; j++){
      sum1 += Bgi[(i*Ng)+j] * SI[(j*2)+0];
      sum2 += Bgi[(i*Ng)+j] * SI[(j*2)+1];
    }
    SIt[(i*2)+0] = sum1;
    SIt[(i*2)+1] = sum2;
  }    
  
   for(k=0; k<(2*Ng); k++){
     i = floor(k/2);
     j = k % 2;
    
     sum += SI[(i*2)+j] * SIt[(i*2)+j];
    
   }  
 
  mod = NORM * sum;
 
}

/*
 Utility method used by the Spectral method fine-tune an initial 
 given community split. 
 */
void SpectralModularity::modifySplit( int countmax ){

  int i,j,k, KK, Ng, count;
  double qmax, qold;

  count = 0;    Ng   = NR_Bgi; KK = (2*Ng);
  qmax  = 0.0;  qold = 0.0;
    
  double Gsi[Ng];
  double GSI[KK];

  for(k=0; k<Ng; k++){ Gsi[k] = si[k]; }
  
  for(k=0; k<KK; k++){
    i = floor(k/2);
    j = k % 2;
    GSI[(i*2)+j] = SI[(i*2)+j];
  }

  maxModularity( qmax );

  while( count < countmax ){

    if( qmax > qold ){

      for(k=0; k<Ng; k++){ Gsi[k] = si[k]; }
      
      for(k=0; k<KK; k++){
	i = floor(k/2);
	j = k % 2;
	GSI[(i*2)+j] = SI[(i*2)+j];
      }

    }

    qold = qmax;
    qmax = 0.0;

    maxModularity(qmax);

    count++;    

  }
  
  for(k=0; k<Ng; k++){ si[k] = Gsi[k]; }

  for(k=0; k<KK; k++){
    i = floor(k/2);
    j = k % 2;    
    SI[(i*2)+j] = GSI[(i*2)+j];
  }

}

/*
 Utility method used by the Spectral method to find 
 which node, when moved gives the maximum change in the 
 Modularity value.
 */
void SpectralModularity::maxModularity(double &qmax){

  int k,Ng,ind_max;
  double Q;
  
  Ng = NR_Bgi;

  double qstored[Ng];
  Q = 0.0;
     
  for(k=0; k<Ng; k++){
    
    qstored[k] = 0.0;
  
    if( visited[k] < 1 ){

      Q  = 0.0;
      
      deltaModularityMax( k, Q );      
      
      qstored[k] = Q;
      	      
    }

  }

  qmax    =   0;//qstored(0);
  ind_max =  -1;//0; 
  for(k=0; k<Ng; k++){
    
    if( qstored[k] > qmax ){
      qmax    = qstored[k];
      ind_max = k; 
    }

  }
  
  if( ind_max != -1 ){
    visited[ind_max] = 1;
    if( si[ind_max] == 1 ){
      si[ind_max] = -1;
      SI[(ind_max*2)+0] = 0;
      SI[(ind_max*2)+1] = 1;
    } else {
      si[ind_max] = 1;
      SI[(ind_max*2)+0] = 1;
      SI[(ind_max*2)+1] = 0;
    }
  } 
  

}


/*
 The change in Modularity used during the fine-tuning
 method; where node si_i is moved from one community to 
 the other: if si^old_i = +-1 => si^new_i = -+1
 deltaQ = deltaQ^new - deltaQ^old
        = Sum_ij { Big_ij * si^new_i * si^new_j }  
        - Sum_ij { Big_ij * si^old_i * si^old_j }
        = Sum_(i!=k,j!=k) { Bgi_ij * si^new_i * si^new_j
	                    + Sum_(j!=k) Big_kj * si^new_k * si^new_j
                            + Sum_(i!=k) Big_ik * si^new_i * si^new_k 
                            + Big_kk } 
        - Sum_(i!=k,j!=k) { Big_ij si^old_i * si^old_j
                            - Sum_(j!=k) Big_kj * si^old_k * si^old_j
	                    - Sum_(i!=k) Big_ik * si^old_i * si^old_k 
			    - Big_kk }
        = Sum_(j!=k) { Big_kj * ( si^new_k - si^old_k ) * si^old_j }
	+ Sum_(i!=k) { Big_ik * si^old_i * ( si^new_k - si^old_k ) }
	=  2 * ( si^new_k - si^old_k ) * Sum_(i!=k) { Big_ik * si^old_i }
	= -4 * si^old_k * Sum_(i!=k) { Big_ik * si^old_i }
*/
void SpectralModularity::deltaModularityMax( int K, double &mod ){

  int k,Ng;
  double sumi;

  Ng   = NR_Bgi;  
  mod  = 0;
  sumi = 0.0;

  for(k=0; k<Ng; k++){

    if( k != K )
      sumi += Bgi[(k*Ng)+K] * si[k];

  }

  mod = -4.0 * si[K] * sumi;


}

/*
 Calculate the split of nodes belonging to the last group of nodes
 with positive eigenvector values.
 */
void SpectralModularity::split( double *Bgiii, int NR_Bgiii, int *keys, const char *sign ){


  int k,i,j,p,Ng,KK,KKin;
  
  if( PRINT ){ cout << "> In " << sign << " method... " << endl; }

  //--- Starting from the group Modularity matrix Bg,
  //--- resize matrices: Bgi, keysi_p, keysi_n, u and betai.
  KKin  = NR_Bgiii * NR_Bgiii;
  Ng    = 0;

  for(k=0; k<NR_Bgiii; k++){
    if( keys[k] != dummy ) Ng++;    
  }

  if( PRINT ){ cout << "> Ng = " << Ng << ". " << endl; }
  
  //creates a new Ng x Ng matrix Bgii from Bgiii
  KK = Ng * Ng;  
  double Bgii[KK];
  
  for(k=0, p=0; k<KKin; k++){

    i = floor(k/NR_Bgiii);
    j = k % NR_Bgiii;

    if( keys[i] != dummy && keys[j] != dummy ){
      Bgii[p] = 0.0;
      Bgii[p] = Bgiii[(i*NR_Bgiii)+j]; p++; }
  }        

  //--- filter out dummy values from keys
  int keysi_p[Ng];
  int keysi_n[Ng];
  
  if( strcmp (sign,"splitP") == 0 ){
    for(k=0, i=0; k<NR_Bgiii; k++){
      if(keys[k] != dummy){ keysi_p[i++] = keys[k]; }
    }
  } else {
    for(k=0, i=0; k<NR_Bgiii; k++){
      if(keys[k] != dummy){ keysi_n[i++] = keys[k]; }    
    }
  }
  //---
  
  //--- Calculate the Modularity matrix Bgi for the new matrix Bgii
  calculateB(Bgii, Ng);

  
  //--- Calculate eigenvectors, and values, from Bgi...
  calculateEigenVectors();

  if( PRINT ){ cout << "> max EigenValue is " << betai << endl; }

  
  if( betai > tol ){
	
    //--- set up the index vectors, si and SI, for the initial split
    maximiseIndexVectors();

    double deltaQ_old = 0.0;
    double deltaQ_new = 0.0;

    int cp = 0;
    int cn = 0;
    
    //--- Calculate the Spectral Modularity
    deltaModularity(deltaQ_old);
    if( PRINT ){ cout << "> Spectral Q: " << deltaQ_old << endl; }

    double diff = fabs(deltaQ_old);
    int count   = 0;

    //--- Fine tuning stage to maximum deltaModularity for the initial split
    visited[Ng];

    //reset visited
    for(k=0; k<Ng; k++){ visited[k]=0; }
    while( diff > tol ){

      modifySplit( Ng );

      deltaModularity(deltaQ_new);
      if( PRINT ){ cout << "> Modified Q: " << deltaQ_new << endl; }

      diff = fabs( deltaQ_new - deltaQ_old ); 
    
      deltaQ_old = deltaQ_new;

    }
    
    //--- Keep recorded of maximum fine-tuned Modularity value.
    specQ += deltaQ_old;
    for(k=0; k<Ng; k++){
      if(si[k] > 0) cp++;
      else          cn++;
    }

    
    //Minimum cluster size... we can reset this either in the header or using setMinCn();
    if( cp < MINCn || cn < MINCn ) { if( PRINT ){ cout << "> Stop splitting. " << endl; }

      return;

    }

    if( strcmp (sign,"splitP") == 0 ){

    //--- get the maximum community number
    int Ncomp = MAXK  + 1;
    MAXK      = Ncomp;
    
    if( PRINT ){ cout << "si[0] " << si[0] << endl; }
    if( PRINT ){ cout << "> node list " << endl; }

    for(k=0; k<Ng; k++){

      if( si[k] > 0 ){
	keysi_p[k] = (int)keysi_p[k];
	keysi_n[k] = dummy;
	gg->V[ (int)keysi_p[k] ].K = Ncomp;
	if( PRINT ){ cout << "> Node: " << gg->V[ (int)keysi_p[k]].label << " c = " << gg->V[ (int)keysi_p[k]].K << endl; }
      } else {
	keysi_n[k] = (int)keysi_p[k];
	keysi_p[k] = dummy;
	if( PRINT ){ cout << "> Node: " << gg->V[ (int)keysi_n[k]].label << " c = " << gg->V[ (int)keysi_n[k]].K << endl; }
      }

    }

    } else {

      //--- get the maximum community number
      int Ncomn = MAXK  + 1;
      MAXK      = Ncomn;
    
      if( PRINT ){ cout << "si[0] " << si[0] << endl; }    
      if( PRINT ){ cout << "> node list " << endl; }
      
      for(k=0; k<Ng; k++){
	
	if( si[k] < 0 ){
	  keysi_n[k] = keysi_n[k];
	  keysi_p[k] = dummy;
	  gg->V[ (int)keysi_n[k] ].K = Ncomn;
	  if( PRINT ){ cout << "> Node: " << gg->V[ (int)keysi_n[k]].label << " c = " << gg->V[ (int)keysi_n[k]].K << endl; }	  
	} else {
	  keysi_p[k] = keysi_n[k];
	  keysi_n[k] = dummy;
	  if( PRINT ){ cout << "> Node: " << gg->V[ (int)keysi_p[k]].label << " c = " << gg->V[ (int)keysi_p[k]].K << endl; }
	}
	
      }
	
      
    }
    
    //--- Recursively split the group of positive eigenvector nodes
    split( Bgii, Ng, keysi_p, "splitP" );

    //--- Recursively split the group of negative eigenvector nodes    
    split( Bgii, Ng, keysi_n, "splitN" );
      
  } else {
    if( PRINT ){ cout << "> Stop splitting. " << endl; }

    
    return; 
  }
  
}


/*
 Kronecker-delta function
 */
int SpectralModularity::delta( int i, int j){ return (i == j ) ? 1 : 0; }

/*
 Calculte the Modularity matrix when split 
 into more than two communities, see [2]
 in method declarations above.
 */
 void SpectralModularity::calculateB(double *B, int NR_B){

  int i,j,k,p,KK,Ng;

  Ng = NR_B;
  KK = Ng * Ng;
  
  if( Bgi != 0 && usedBgi == true ){

    free(Bgi);

    Bgi = (double*)malloc(KK*sizeof(double));

    NR_Bgi = Ng;
    NC_Bgi = Ng;
    
  } else {

    Bgi = (double*)malloc(KK*sizeof(double));

    NR_Bgi = Ng;
    NC_Bgi = Ng;
    
    usedBgi = true;

  }

  for(i=0; i<Ng; i++){
    for(j=0; j<Ng; j++){
      double sum = 0.0;
      for(k=0; k<Ng; k++){ sum += B[(i*Ng)+k]; }

      Bgi[(i*Ng)+j] = 0.0;
      Bgi[(i*Ng)+j] = B[(i*Ng)+j] -1.0 * delta(i,j) * sum;
    }
  }
  
 
 }

void SpectralModularity::setMinCn( int NEWCn ){

  if( NEWCn > 0 && NEWCn <= gg->getN() ){
    MINCn = NEWCn;
    if( PRINT ){ cout << "> Min Cn = " << MINCn << endl; }
  }
  
}


void SpectralModularity::settol( double NEWtol ){

  if( NEWtol >= 0 ){
    tol = NEWtol;
    if( PRINT ){ cout << "> tol = " << tol << endl; }
  }
  
}

void SpectralModularity::setPrint( bool status ){

  PRINT = status;
  
  
}
