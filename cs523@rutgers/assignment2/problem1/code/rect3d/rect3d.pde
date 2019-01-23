int N=300;
float a=50;
float b=50;
float c=50;
float alpha=75;
float CO=50;
float edge;
float epsilon=1e-2;
float[] vertex_phi_value=new float[(N+2)*(N+2)*(N+2)];

void setup() {
  
  size(800,800);
  background(255); 
  // assign phi value to each vertex
  edge=width/(N+2.0);
  for(int i=0;i<N+2;i++){
    for(int j=0;j<N+2;j++){
      for(int k=0;k<N+2;k++){
      vertex_phi_value[index_in_array(i,j,k)]=abs(phi_function(i*edge,j*edge,k*edge));
      }
    }
  }
  


 
 
  for(int i=1;i<=10;i++){
 
  float x=random(0,width);
  float y=random(0,height);
  float z=random(0,height);
  float ox=x,oy=y,oz=z;
  

  float dx;
  float dy;
  float dz;
  float dl;


  if(abs(phi_interpolation(x,y,z))<epsilon){
    print("Found: \n");
    print("(", x, ",",y,",",z,")\n");
  }
  
  else{
    while(abs(phi_interpolation(x,y,z))>epsilon){
      //print("x: ",x,"  y: ",y,"  z: ",z,"  phi: ",phi_interpolation(x,y,z),"actual value: ",phi_function(x,y,z),"\n");
      
      dx=updateX(x,y,z);
      dy=updateY(x,y,z);
      dz=updateZ(x,y,z);
      dl=sqrt(dx*dx+dy*dy+dz*dz);
      dx=dx/dl;
      dy=dy/dl;
      dz=dz/dl;
      x-=1*abs(phi_interpolation(x,y,z))*dx;
      y-=1*abs(phi_interpolation(x,y,z))*dy;
      z-=1*abs(phi_interpolation(x,y,z))*dz;
      
    }
     print("Test ",i,"\n");
    print("nearest: (",x, ",",y,",",z,")\n");
    print("original: (",ox,",",oy,",",oz,")\n");
    print("phi: ",phi_interpolation(x,y,z),"\n");
    print("actual value: ",phi_function(x,y,z),"\n");
    print("\n");
    //print(vertex_phi_value[index_in_array(x_index(x),y_index(y))]);
  }

}

}






float phi_function(float x_, float y_, float z_)
{  
  float v1=abs(x_-(width/2.0-a))*b*c;
  float v2=abs(x_-(width/2.0+a))*b*c;
  float v3=abs(y_-(height/2.0-b))*a*c;
  float v4=abs(y_-(height/2.0+b))*a*c;
  float v5=abs(z_-(height/2.0-c))*a*b;
  float v6=abs(z_-(height/2.0+c))*a*b;

  
  float ph=(v1+v2+v3+v4+v5+v6)/(2*3*a*b*c)-1;
  
  if(ph>0)  return ph;
  else{
  float min1=v1<v2?v1:v2;
  float min2=v3<v4?v3:v4;
  float min3=v5<v6?v5:v6;
  float min=min1<min2?(min1<min3?min1:min3):(min2<min3?min2:min3);

  //print("min/V: ",min/(8*a*b*c));
  return -min/(8*a*b*c);
  
  
  }

}



float phi_interpolation(float x_,float y_, float z_)
{
  int i=x_index(x_);
  int j=y_index(y_); 
  int k=z_index(z_);
  
  float dx=x_-i*edge;
  float dy=y_-j*edge;
  float dz=z_-k*edge;
  
  float p1=vertex_phi_value[index_in_array(i,j,k)];
  float p2=vertex_phi_value[index_in_array(i+1,j,k)];
  float p3=vertex_phi_value[index_in_array(i,j+1,k)];
  float p4=vertex_phi_value[index_in_array(i+1,j+1,k)];
  //print("p1: ",p1,", p2: ",p2,", p3: ",p3,", p4: ",p4,"\n");
  
  float pa=dx*p2+(1-dx)*p1;
  float pb=dx*p4+(1-dx)*p3;
  //print("pa: ",pa,", pb: ",pb,"\n");
  float p=dy*pb+(1-dy)*pa;
  
  float p5=vertex_phi_value[index_in_array(i,j,k+1)];
  float p6=vertex_phi_value[index_in_array(i+1,j,k+1)];
  float p7=vertex_phi_value[index_in_array(i,j+1,k+1)];
  float p8=vertex_phi_value[index_in_array(i+1,j+1,k+1)];
  //print("p5: ",p5,", p6: ",p6,", p7: ",p7,", p8: ",p8,"\n");
  
  float pc=dx*p6+(1-dx)*p5;
  float pd=dx*p8+(1-dx)*p7;
  //print("pc: ",pc,", pd: ",pd,"\n");
  float pn=dy*pb+(1-dy)*pa;
  
  float preturn=dz*pn+(1-dz)*p;
  //print("preturn: ",preturn,"\n\n");
  return preturn;

}

float updateX(float x_, float y_, float z_)
{
  int x=x_index(x_);
  int y=y_index(y_);
  int z=z_index(z_);
  
  float div_x=(vertex_phi_value[index_in_array(x+1,y,z)]
                  -vertex_phi_value[index_in_array(x,y,z)])/edge;  
  return div_x;
}
float updateY(float x_, float y_, float z_)
{
    int x=x_index(x_);
    int y=y_index(y_);
    int z=z_index(z_);
  
    float div_y=(vertex_phi_value[index_in_array(x,y+1,z)]
                  -vertex_phi_value[index_in_array(x,y,z)])/edge;  
    return div_y;
}

float updateZ(float x_, float y_, float z_)
{
    int x=x_index(x_);
    int y=y_index(y_);
    int z=z_index(z_);
  
    float div_z=(vertex_phi_value[index_in_array(x,y,z+1)]
                  -vertex_phi_value[index_in_array(x,y,z)])/edge;  
    return div_z;
}


int x_index(float x_)
{
  return floor(x_/edge);
}

int y_index(float y_)
{
  return floor(y_/edge);
}

int z_index(float z_)
{
  return floor(z_/edge);
}

int index_in_array(int i_,int j_,int k_)
{
  return i_*(N+2)*(N+2)+j_*(N+2)+k_;
}
