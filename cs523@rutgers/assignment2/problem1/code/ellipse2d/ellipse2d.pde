int N=3000;
float a=100;
float b=100;
float alpha=75;
float CO=50;
float edge;
float epsilon=1e-3;
float[] vertex_phi_value=new float[(N+2)*(N+2)];

void setup() {
  
  size(800,800);
  background(255); 
  // assign phi value to each vertex
  edge=width/(N+2.0);
  for(int i=0;i<N+2;i++){
    for(int j=0;j<N+2;j++){
      vertex_phi_value[index_in_array(i,j)]=abs(phi_function(i*edge,j*edge));
    }
  }
  

  rectMode(CENTER);
  float grey;
  for(int i=1;i<=N;i++){
    for(int j=1;j<=N;j++){
      grey=phi_function(i*edge,j*edge);
      //grey=phi(i*edge,j*edge);
      //print(grey);
      //print("\n");
      noStroke();
      fill(grey*CO,grey*CO,grey*CO);
      rect(i*edge,j*edge,edge,edge);
    }
  }
 ellipseMode(CENTER);
 blendMode(BLEND);
 noFill();
 strokeWeight(1);
 stroke(255);
 ellipse(width/2, height/2, 2*a, 2*b);
 
 for(int i=1;i<=10;i++){
  float x=random(0,width);
  float y=random(0,height);
  float ox=x;
  float oy=y;
  
  

  float dx;
  float dy;
  float dl;


  if(abs(phi_interpolation(x,y))<epsilon){
    print("Found: \n");
    print("(", x, ",",y,")\n");
  }
  
  else{
    while(abs(phi_interpolation(x,y))>epsilon){
      //print("x: ",x,"  y: ",y,"  phi: ",phi_interpolation(x,y),"actual value: ",phi_function(x,y),"\n");
      
      stroke(0);
      strokeWeight(3);
      point(x,y);
      dx=updateX(x,y);
      dy=updateY(x,y);
      dl=sqrt(dx*dx+dy*dy);
      dx=dx/dl;
      dy=dy/dl;
      x-=abs(phi_interpolation(x,y))*dx;
      y-=abs(phi_interpolation(x,y))*dy; 
      
    }
    print("Test ",i,"\n");
    print("nearest: (", x, ",",y,")\n");
    print("original: (",ox,",",oy,")\n");
    print("phi: ",phi_interpolation(x,y),"\n");
    print("actual value: ",phi_function(x,y),"\n");
    print("\n");
    stroke(255,0,0);
    strokeWeight(5);
    point(x,y);
    //print(vertex_phi_value[index_in_array(x_index(x),y_index(y))]);
  }
  saveFrame("circle.jpg");
}


}





float phi_function(float x_, float y_)
{
  float ph= sqrt(pow((x_-width/2)/a,2)+pow((y_-height/2)/b,2))-1;
  //float ph= (x_-width/2.0)*(x_-width/2.0)/(a*a)+(y_-height/2.0)*(y_-height/2.0)/(b*b)-1.0;
  return ph;
}

float phi_interpolation(float x_,float y_)
{
  int i=x_index(x_);
  int j=y_index(y_);  
  
  float dx=x_-i*edge;
  float dy=y_-j*edge;
  
  float p1=vertex_phi_value[index_in_array(i,j)];
  float p2=vertex_phi_value[index_in_array(i+1,j)];
  float p3=vertex_phi_value[index_in_array(i,j+1)];
  float p4=vertex_phi_value[index_in_array(i+1,j+1)];
  //print("p1: ",p1,", p2: ",p2,", p3: ",p3,", p4: ",p4,"\n");
  
  float pa=dx*p2+(1-dx)*p1;
  float pb=dx*p4+(1-dx)*p3;
  //print("pa: ",pa,", pb: ",pb,"\n");
  float p=dy*pb+(1-dy)*pa;
  //print("p: ",p,"\n\n");
  return p;

}

float updateX(float x_, float y_)
{
  int x=x_index(x_);
  int y=y_index(y_);
  
  float div_x=(vertex_phi_value[index_in_array(x+1,y)]
                  -vertex_phi_value[index_in_array(x,y)])/edge;  
  return div_x;
}
float updateY(float x_, float y_)
{
    int x=x_index(x_);
    int y=y_index(y_);
  
    float div_y=(vertex_phi_value[index_in_array(x,y+1)]
                  -vertex_phi_value[index_in_array(x,y)])/edge;  
    return div_y;
}

int x_index(float x_)
{
  return floor(x_/edge);
}

int y_index(float y_)
{
  return floor(y_/edge);
}

int index_in_array(int i_,int j_)
{
  return j_*(N+2)+i_;
}
