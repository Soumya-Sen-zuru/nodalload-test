#include<stdio.h>
#include<iostream>
#include<vector>
#include<cmath>
#include<iomanip>
using namespace std;
double area_load = 10.0;                     // Partial load intensity
 
double nodesX[] = {1.0, 6.0, 5.3, 1.5};      // Plate element node coordinates
double nodesY[] = {0.5, 0.0, 5.0, 4.0};
double load_nodesX[] = {1.25, 3.5, 3.25, 1.5};    // Partial load coordinates
double load_nodesY[] = {1.5, 1.2, 3.0, 2.5};



double * intersection(double ax, double ay, double bx, double by, double cx , double cy,double dx, double dy)
{ static double pointofintersection[1];
  // subtracting co-ordinates of point A from
    // b and c, to make A as origin
    bx -= ax;
    by -= ay;
    cx -= ax;
    cy -= ay;
 	dx -= ax;
    dy -= ay;
    // Determining cross Product
    double cross_product1 = bx * cy - by * cx;
  	 double cross_product2 = bx * dy - by * dx;
    
    if(((cross_product1 > 0.0 )&&(cross_product2 < 0.0))||((cross_product1 < 0.0) &&(cross_product2 > 0.0))) 
    {
    	
    	 double dx = bx - ax;
   		 double dy = by - ay;
		double m1 = dy / dx;
    
    double c1 = ay - m1 * ax;	
    dx = dx - cx;
    dy = dy - cy;


    double m2 = dy / dx;	
    double c2 = cy - m2 * cx;	
    pointofintersection[0] = (c2 - c1) / (m1 - m2);

         pointofintersection[1]= m1 * pointofintersection[0] + c1;	
    
	}
        
  
  
    else if (cross_product1 == 0)
    {
		   pointofintersection[0]=cx ;
    pointofintersection[1]=cy ;
  } 
  else if (cross_product2 == 0)
    {
		   pointofintersection[0]=dx ;
    pointofintersection[1]=dy ;
  } 
  else
  pointofintersection[0]= pointofintersection[1]=0;
     
    return  pointofintersection; 
}

void loadonplate()
{
int j=0,i=0;
while(i<4)
{while(j<4)
{double *p=intersection(nodesX[i],nodesY[i],nodesX[i+1],nodesY[i+1],load_nodesX[j],load_nodesY[j],load_nodesX[j+1],load_nodesY[j+1]);
if(*p!=0.0|| *(p+1)!=0)
load_nodesX[j]=*p;
load_nodesY[j]=*(p+1);
j++;
}
i++;
}
}
struct shape_func {
   double n0, n1, n2, n3;
   double N0, N1, N2, N3;
   double j;
} typedef shape_func;
 
shape_func ShapeFunc(double &a, double &b) {
   shape_func N;
   const int n = *(&nodesX + 1) - nodesX;      // Total sides of polygon
 
   double x, y;
   double A, B, C;
   double u, v;
 
   N.n0 = 0.25*(1.0-a)*(1.0-b);        // Shape functions in (s,t) coordinates
   N.n1 = 0.25*(1.0+a)*(1.0-b);
   N.n2 = 0.25*(1.0+a)*(1.0+b);
   N.n3 = 0.25*(1.0-a)*(1.0+b);
   x = N.n0*load_nodesX[0] + N.n1*load_nodesX[1] + N.n2*load_nodesX[2] + N.n3*load_nodesX[3];
   y = N.n0*load_nodesY[0] + N.n1*load_nodesY[1] + N.n2*load_nodesY[2] + N.n3*load_nodesY[3];
 
   double ax = 0.25*(nodesX[0] + nodesX[1] + nodesX[2] + nodesX[3]);
   double bx = 0.25*(-nodesX[0] - nodesX[1] + nodesX[2] + nodesX[3]);
   double cx = 0.25*(-nodesX[0] + nodesX[1] + nodesX[2] - nodesX[3]);
   double dx = 0.25*(nodesX[0] - nodesX[1] + nodesX[2] - nodesX[3]);
 
   double ay = 0.25*(nodesY[0] + nodesY[1] + nodesY[2] + nodesY[3]);
   double by = 0.25*(-nodesY[0] - nodesY[1] + nodesY[2] + nodesY[3]);
   double cy = 0.25*(-nodesY[0] + nodesY[1] + nodesY[2] - nodesY[3]);
   double dy = 0.25*(nodesY[0] - nodesY[1] + nodesY[2] - nodesY[3]);
 
   // Area load coordinates on isoparametric element
   A = by*dx-bx*dy;
   B = dy*(x-ax)-bx*cy+by*cx-dx*(y-ay);
   C = cy*(x-ax)-cx*(y-ay);
 
   if(dx == 0.0 || dy == 0.0){
       v = (x-ax)/cx;
       u = (y-ay)/by;
   } else {
       v = (-B + sqrt(B*B - 4.0*A*C))/(2.0*A);
       u = (x - ax - bx*v)/(cx + dx*v);
   }
 
   N.N0 = 0.25*(1.0-u)*(1.0-v);       // Shape functions in (u,v) coordinates
   N.N1 = 0.25*(1.0+u)*(1.0-v);
   N.N2 = 0.25*(1.0+u)*(1.0+v);
   N.N3 = 0.25*(1.0-u)*(1.0+v);
   // Jacobian of the Load transformation
   double j_11 = 0.25*((-1.0+b)*load_nodesX[0] + (1.0-b)*load_nodesX[1] + (1.0+b)*load_nodesX[2] + (-1.0-b)*load_nodesX[3]);   //dxi/ds
   double j_12 = 0.25*((-1.0+b)*load_nodesY[0] + (1.0-b)*load_nodesY[1] + (1.0+b)*load_nodesY[2] + (-1.0-b)*load_nodesY[3]);   //deta/ds
   double j_21 = 0.25*((-1.0+a)*load_nodesX[0] + (-1.0-a)*load_nodesX[1] + (1.0+a)*load_nodesX[2] + (1.0-a)*load_nodesX[3]);   //dxi/dt
   double j_22 = 0.25*((-1.0+a)*load_nodesY[0] + (-1.0-a)*load_nodesY[1] + (1.0+a)*load_nodesY[2] + (1.0-a)*load_nodesY[3]);   //deta/dt
   N.j = j_11*j_22 - j_12*j_21;
   return N;
}
int main() {
	loadonplate();
   const int m = 3;        // No. of integration points
 
   double sum0 = 0.0;
   double sum1 = 0.0;
   double sum2 = 0.0;
   double sum3 = 0.0;
 
   // Quadrature weights
   double w[] = {5.0/9.0, 8.0/9.0, 5.0/9.0};
   // Quadrature nodes
   double p[] = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};
   // Numerical integration: Gauss-Legendre quadrature with nine (3x3) integration points
   for(int i=0; i<m; i++){
       for(int j=0; j<m; j++){
           shape_func SF = ShapeFunc(p[i], p[j]);
           sum0 += w[i]*w[j]*SF.N0*SF.j*area_load;
           sum1 += w[i]*w[j]*SF.N1*SF.j*area_load;
           sum2 += w[i]*w[j]*SF.N2*SF.j*area_load;
           sum3 += w[i]*w[j]*SF.N3*SF.j*area_load;
       }
   }
   cout << "Nodal Forces:" << endl;
   cout << setprecision(10) << "Node 1: " << sum0 << endl;
   cout << "Node 2: " << sum1 << endl;
   cout << "Node 3: " << sum2 << endl;
   cout << "Node 4: " << sum3 << endl;
   cout << "Check for total Load: " << sum0+sum1+sum2+sum3 << endl;     
}
