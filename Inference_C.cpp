#include <Rcpp.h> 
#include <iostream>
#include <vector>
#include <string>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericMatrix Inference_C( std::vector<double> p1,
                                 std::vector<double> p2,
                                 std::vector<int> id1,
                                 std::vector<int> org_id,
                                 double pi0 ){
  
  // Step 0. Settings
  int m = p1.size();
  double m0 = pi0*p1.size();
  
  // DEG and FDR
  int i0 = -1; int j0 = -1;
  int DEG = 0;
  double prev_FDR = 0;  double FDR = 0;
  int Num_pts[m]; int output [m][2];
  
  // Initialize Num_pts & output variables :
  for(int i = 0; i < m; i++){
    Num_pts[i] = 1;
    for(int j = 0; j < 2; j++){
      output[i][j] = -2; // i, j
    }    
  }
  
  // Step 1. generate L_ids & R_ids.
  // std::vector<int> L_ids [m]; 
  std::vector<int> R_ids [m];
  for(int i = 0; i < m; i++){
    for(int j = 0; j < i; j++){
      if(id1[j]>id1[i]){ R_ids[i].push_back(j); }else{ Num_pts[i]++; }
    }
  }
  
  // Step 2.
  for(int i = 0; i < m; i++){
    
    // Rectangular area.
    DEG = Num_pts[i];
    FDR = (p2[i]*p1[i])*m0/DEG;
    
    i0 = output[(DEG-1)][0];  j0 = output[(DEG-1)][1];
    if(i0==-2){ i0 = i; j0 = -1;} // In case that there is no previous id.
    
    if(j0 == -1){
      prev_FDR = (p2[i0]*p1[i0])*m0/Num_pts[i0];
    }else{
      prev_FDR = ((p2[i0]*p1[i0])+(p2[R_ids[i0][j0]]*(1-p1[i0])))*m0/(Num_pts[i0]+j0+1);
    }
    
    if(prev_FDR >= FDR){ // better than update the output.
      output[(DEG-1)][0] = i;
      output[(DEG-1)][1] = -1;
    }
    
    
    // L-shape area.
    for(int j = 0; j < R_ids[i].size(); j++){
      
      DEG = (Num_pts[i]+j+1);
      FDR = ((p2[i]*p1[i])+(p2[R_ids[i][j]]*(1-p1[i])))*m0/DEG;
      
      i0 = output[(DEG-1)][0];  j0 = output[(DEG-1)][1];
      if(i0==-2){ i0 = i; j0 = j;} // In case that there is no previous id.
      
      if(j0 == -1){
        prev_FDR = (p2[i0]*p1[i0])*m0/Num_pts[i0];
      }else{
        prev_FDR = ((p2[i0]*p1[i0])+(p2[R_ids[i0][j0]]*(1-p1[i0])))*m0/(Num_pts[i0]+j0+1);
      }
      
      if(prev_FDR >= FDR){
        output[(DEG-1)][0] = i;
        output[(DEG-1)][1] = j;
      }
    }
    
  }
  
  // Step 3. output : DEG, FDR, i, j
  Rcpp::NumericMatrix NM_output(m,4);
  for(int s = 0; s < m; s++){
    
    NM_output(s,0) = (s+1); // DEG
    i0 = output[s][0]; j0 = output[s][1];
    
    if(j0 == -1){
      
      NM_output(s,1) = (p2[i0]*p1[i0])*m0/Num_pts[i0]; // FDR
      NM_output(s,2) = org_id[i0]; // org_id[i]
      NM_output(s,3) = -1; // org_id[j]
      
    }else{
      
      NM_output(s,1) = ((p2[i0]*p1[i0])+(p2[R_ids[i0][j0]]*(1-p1[i0])))*m0/(Num_pts[i0]+j0+1); // FDR
      NM_output(s,2) = org_id[i0]; // org_id[i]
      NM_output(s,3) = org_id[R_ids[i0][j0]]; // org_id[j]
      
    }
    
  }
  
  return NM_output;
}


// // [[Rcpp::export]]
// Rcpp::NumericMatrix Inference_C_int( std::vector<double> p1,
//                                      std::vector<double> p2,
//                                      std::vector<int> id1,
//                                      std::vector<int> org_id,
//                                      double pi0 ){
//   
//   // Step 0. Settings
//   int m = p1.size();
//   double m0 = pi0*p1.size();
//   double int_N = 100000000.0; // double to int, using int(v*int_N)
//   // 2147483647
//   // DEG and iFDR
//   int DEG = 0; unsigned int iFDR = 0;
//   int Num_pts[m]; int output [m][3];
//   
//   // Initialize Num_pts & output variables :
//   for(int i = 0; i < m; i++){
//     Num_pts[i] = 1;
//     output[i][0] = int(int_N*10); // iFDR
//     for(int j = 1; j < 3; j++){
//       output[i][j] = -2; // i, j, k
//     }
//   }
//   
//   // Step 1. generate L_ids & R_ids.
//   // std::vector<int> L_ids [m];
//   std::vector<int> R_ids [m];
//   
//   for(int i = 0; i < m; i++){
//     for(int j = 0; j < i; j++){
//       if(id1[j]>id1[i]){ R_ids[i].push_back(j);}else{ Num_pts[i]++; }
//     }
//   }
//   
//   // Step 2.
//   for(int i = 0; i < m; i++){
//     
//     // Rectangular area.
//     DEG = Num_pts[i];
//     iFDR = int(int_N*(p2[i]*p1[i])*m0/DEG);
//     
//     if(output[(DEG-1)][0]>iFDR){
//       output[(DEG-1)][0] = iFDR;
//       output[(DEG-1)][1] = org_id[i];
//       output[(DEG-1)][2] = -1;
//     }
//     
//     // L-shape area.
//     for(int j = 0; j < R_ids[i].size(); j++){
//       
//       DEG = (Num_pts[i]+j+1);
//       iFDR = int(int_N*((p2[i]*p1[i])+(p2[R_ids[i][j]]*(1-p1[i])))*m0/DEG);
//       
//       if(output[(DEG-1)][0]>iFDR){
//         output[(DEG-1)][0] = iFDR;
//         output[(DEG-1)][1] = org_id[i];
//         output[(DEG-1)][2] = org_id[R_ids[i][j]];
//       }
//     }
//   }
//   
//   // Step 3. datatype change from vector to NM.
//   Rcpp::NumericMatrix NM_output(m,3);
//   for(int s = 0; s < m; s++){
//     NM_output(s,0) = output[s][0]/int_N;
//     for(int t = 1; t < 3; t++){
//       NM_output(s,t) = output[s][t];
//     }
//   }
//   return NM_output;
// }


