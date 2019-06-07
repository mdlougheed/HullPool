// Java Libraries
import java.util.*;
import java.util.List;


// A holder class for a centroid and area
class CenterArea{
Vec2 centroid;
float area;
  
   CenterArea(){
     centroid = new Vec2();
     area=0;
   }   
 
}

// Sutherland-Hodgman polygon clipping algorightm
// Adapted from http://rosettacode.org/wiki/Sutherland-Hodgman_polygon_clipping
class SutherlandHodgmanClipper{
ArrayList<Vec2> subject;
CenterArea subjectCenterArea = new CenterArea();
 
ArrayList<Vec2> clipper;
CenterArea clipperCenterArea = new CenterArea();
 
ArrayList<Vec2> result;
Vec2 resultAsArray[];
CenterArea resultCenterArea = new CenterArea();

   SutherlandHodgmanClipper(Vec2[] subjPoints, Vec2[] clipPoints){
   subject = new ArrayList<Vec2>(Arrays.asList(subjPoints));
   clipper = new ArrayList<Vec2>(Arrays.asList(clipPoints));
   result  = new ArrayList<Vec2>(subject);
   
     clipPolygon();
     resultAsArray=result.toArray(new Vec2[result.size()]);
    
    
     if(subjPoints.length >= 3)
     computeCentroidToOut(subjPoints, subjPoints.length, subjectCenterArea);
    
    
     if(clipPoints.length >= 3)
     computeCentroidToOut(clipPoints, clipPoints.length, clipperCenterArea);
     
     if(resultAsArray.length >= 3)
     computeCentroidToOut(resultAsArray, resultAsArray.length, resultCenterArea);
   }


// Centroid and area modified from JBox2D PolygonShape.java
  private final void computeCentroidToOut(Vec2[] vs, int count, CenterArea outArea) {
  assert(count >= 3);
  
  Vec2 out = new Vec2();
  out.setZero();
  float area = 0.0f;
  
  if( count == 2){
    out.set(vs[0]).addLocal(vs[1]).mulLocal(.5f);
    return;
  }
  
  // pRef is the reference point for forming triangles.
  // It's location doesn't change the result (except for rounding error).
  Vec2 pRef = new Vec2();
  pRef.setZero();
  
  Vec2 e1 = new Vec2();
  Vec2 e2 = new Vec2();
  
  float inv3 = 1.0f / 3.0f;
  
  for (int i = 0; i < count; ++i){
    // Triangle vertices.
    Vec2 p1 = pRef;
    Vec2 p2 = vs[i];
    Vec2 p3 = i + 1 < count ? vs[i+1] : vs[0];
  
    e1.set(p2).subLocal(p1);
    e2.set(p3).subLocal(p1);
  
    float D = Vec2.cross(e1, e2);
  
    float triangleArea = 0.5f * D;
    area += triangleArea;
  
    // Area weighted centroid
    e1.set(p1).addLocal(p2).addLocal(p3).mulLocal(triangleArea * inv3);
    out.addLocal(e1);
  }
  
  // Area
//  assert(area > Settings.EPSILON);
  if(area >= Settings.EPSILON)
    outArea.area = area;
  else {
    outArea.area = Settings.EPSILON;
    area = Settings.EPSILON;
  }


  // Centroid
  assert(area >= Settings.EPSILON);
  out.mulLocal(1.0f / area);
  outArea.centroid.set(out.x, out.y);
  } 
 
 
 
  private void clipPolygon() {
  int len = clipper.size();
    for (int i = 0; i < len; i++) {
   
      int len2 = result.size();
      ArrayList<Vec2> input = new ArrayList<Vec2>(result);

      result = new ArrayList<Vec2>();
   
      Vec2 A = clipper.get((i + len - 1) % len);
      Vec2 B = clipper.get(i);
   
      for (int j = 0; j < len2; j++) {
        Vec2 P = input.get((j + len2 - 1) % len2);
        Vec2 Q = input.get(j);
   
        if (isInside(A, B, Q)) {
          if (!isInside(A, B, P))
            result.add(intersection(A, B, P, Q));

          result.add(Q);
          } 
          else if (isInside(A, B, P))
            result.add(intersection(A, B, P, Q));
      }
    }
  }
 

  private boolean isInside(Vec2 a, Vec2 b, Vec2 c) {
    return (a.x - c.x) * (b.y - c.y) > (a.y - c.y) * (b.x - c.x);
  }
 

  private Vec2 intersection(Vec2 a, Vec2 b, Vec2 p, Vec2 q) {
  float A1 = b.y - a.y;
  float B1 = a.x - b.x;
  float C1 = A1 * a.x + B1 * a.y;
   
  float A2 = q.y - p.y;
  float B2 = p.x - q.x;
  float C2 = A2 * p.x + B2 * p.y;
   
  float det = A1 * B2 - A2 * B1;
  float x = (B2 * C1 - B1 * C2) / det;
  float y = (A1 * C2 - A2 * C1) / det;
   
  return new Vec2(x, y);
  }
   
}