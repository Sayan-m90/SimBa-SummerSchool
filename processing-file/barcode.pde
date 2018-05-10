import controlP5.*;
import processing.pdf.*;

// Press 'p' to save a pdf
ControlP5 cp5;
String fileName = "pers";
//String fileName = "spiral_4D_2k_dim4_rho3_log";
//String fileName = "spiral_4D_2k_dim4_rho3_zeta0_9_log";
//String fileName = "spiral_4D_2k_dim4_standard_rips_threshold0_625_log";
//boolean logScale = true;
boolean logScale = false;
int scaling = 3;
int dpi = 200 * scaling;
int w = 3 * dpi - 5;
int h = dpi * 7 / 8;
int boxLeft = dpi / 4;
int boxHeight = h * 6 / 7;
float xMin = 0;
float xMax = 151;
int barCount = 0;
float barWeight = 3;
float[] barCounts = {1,1,1,1};
float[][] diagram = new float[3000000][3];
color barColor = color(0,0,128);
float NEGINF = -100000000;
float INF = 100000000;
int minVal = 0;
int maxVal = 151;
int dim = 2;

void setup() {
  smooth();
  size(w+5,4*h);
  //surface.setSize(w+5, h);
  
  //size(400, 200);
  //surface.setResizable(true);
  
  cp5 = new ControlP5(this);
  cp5.addSlider("xMin")
      .setPosition(50, h)
      .setWidth(300)
      .setRange(minVal,maxVal)
      .setNumberOfTickMarks(41)
      .snapToTickMarks(true)
      .showTickMarks(true);
  cp5.addSlider("xMax")
      .setPosition(50, h+50)
      .setWidth(300)
      .setRange(minVal,maxVal)
      .setNumberOfTickMarks(41)
      .snapToTickMarks(true)
      .showTickMarks(true);
//  cp5.addSlider("xMax").setRange(-10,10);

  loadDiagram();
  drawDiagram();
}

void draw() {
  drawDiagram();
}

void drawDiagram() {
  background(255);
 // surface.setSize(w+5, h);

  ///// BOX
  strokeWeight(1);
  noFill();
  stroke(192);
  for (int i = 1; i < dim+1; ++i) {
    float y = yScale(i, 0);
    line(boxLeft,y,w,y);
  }
  stroke(128);
  strokeWeight(2);
  rect(boxLeft,0, w - boxLeft, boxHeight);
  
  ////// Betti labels
  fill(0);
  for (int i = 0; i < dim+1; ++i) {
    textFont(createFont("Times New Roman",12*scaling,true));
    text("Betti", 3*scaling, yScale(i, 0.55));
    textFont(createFont("Times New Roman",8*scaling,true));
    text(""+i, 27*scaling, yScale(i, 0.6));    
  }
  
  ///// x marks
  /*
  textFont(createFont("Times New Roman",(dpi * 2 / 25),true));
  for(int i = (int)xMin; i <= (int)xMax; ++i) {
    float x = xScale(i);
    line(x, boxHeight, x, boxHeight - (10*scaling));
    String label = "" + i;
    text(label, x - (3*scaling) * label.length(), boxHeight + (10*scaling));
  }
  */
  ///// BARS
  stroke(barColor);
  fill(barColor);
  strokeWeight(barWeight);
//  strokeCap(ROUND);
  for (int d  = 0; d < dim + 1; ++d) {
    int currentDimCount = 1;    
    for (int i = 0; i < barCount; ++i) {
      if ((int)diagram[i][0] == d) {
        drawBar(diagram[i][0], currentDimCount/barCounts[(int)diagram[i][0]], diagram[i][1], diagram[i][2]);
        currentDimCount++;
      }
    }
  }

}

void keyPressed() {
  if (key == 'p') {
    beginRecord(PDF, fileName + ".pdf");
    drawDiagram();
  endRecord();
  } 
}

int loadDiagram() {
    String bar;
    BufferedReader reader = createReader(fileName+ ".dgm");
    
    while(true) {
      try {
        bar  = reader.readLine();
      } catch (IOException e) {
        e.printStackTrace();
        bar = null;
      }
      if (bar == null) {
        return barCount;
      }
  
      String[] barDetails = bar.split(" ");
      diagram[barCount][0] = Float.parseFloat(barDetails[0]);
      diagram[barCount][1] = (barDetails[1].equals("-inf")) ? NEGINF : log2(Float.parseFloat(barDetails[1]));
      diagram[barCount][2] = (barDetails[2].equals("inf")) ? INF : log2(Float.parseFloat(barDetails[2]));
      barCounts[(int)diagram[barCount][0]]++;
      barCount++;
    }
}

void drawBar (float dim, float yShift, float birth, float death) {
  //if(death - birth <= 2)
    //return;
  boolean flag = false;
  birth = (birth == NEGINF) ? xMin : birth;
  birth = (birth > xMax) ? xMax : birth;
  if(death == INF)
  {
    death = xMax;
    flag = true;
  }
  //death = (death == INF) ? xMax : death;
  death = (death < xMin) ? xMin : death;
  // Deal with birth < xMin and death > xMax for finite bars

  float y = yScale(dim, (yShift * 0.8 + 0.1));
  float bx = xScale(birth);
  float dx = xScale(death);
  line(bx, y, xScale(death), y);
  float t = 8;
  
  if (birth < xMin) {
    triangle(bx,y, bx + 2*t, y- t, bx + 2*t, y + t); 
  }
  if (flag == true) {
    triangle(dx,y, dx - 2*t, y-t, dx - 2*t, y + t); 
  }
  
}

float yScale(float dim, float u) {
  return (boxHeight /3) * ((dim) + u);
}

float xScale(float u) {
  float t = (u - xMin) / (xMax- xMin);
  return (1-t) * boxLeft + t * w;
}

float log2(float x) {
  if (logScale == true) {
    return log(x)/log(2.0); 
  }
  else {
    return x;
  }
}
